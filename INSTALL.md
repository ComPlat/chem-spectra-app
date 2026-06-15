# INSTALL

This guidance is tested on Linux Ubuntu 18.04 & Ubuntu 20.04

## 1. Installation

### 1.1. Install Python

Use the file pyproject.toml to determine the version of Python required.

### 1.2. Clone repository and create env

```sh
git clone --recurse-submodules https://github.com/ComPlat/chem-spectra-app.git
cd chem-spectra-app
```

If the repository was cloned without submodules, initialize the test datasets:

```sh
git submodule update --init --recursive
```

The centralized test datasets live in the `test-datasets/` submodule
([chem-spectra-test-files](https://github.com/ComPlat/chem-spectra-test-files)).
Tests resolve dataset files automatically via `tests/dataset_catalog.py`.
To use a different checkout, set `CHEMSPECTRA_TEST_FILES` to its root path.

ALL the **FOLLOWING** commands are assumed to be executed while inside this
`chem-spectra-app` folder.

```
python3 -m venv .venv && .venv/bin/pip install --upgrade pip
source .venv/bin/activate
pip install .
```

### 1.3. Use docker to run a supporting service

```sh
docker pull proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses
mkdir -p chem_spectra/tmp
chmod -R 755 chem_spectra/tmp
docker run --detach --name msconvert_docker \
    --rm -it \
    -e WINEDEBUG=-all \
    -v ./chem_spectra/tmp:/data \
    proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses bash
```

### 1.4. Add `config.py`

Generate a secret key.

```sh
python -c 'import os; print(os.urandom(16))'
>> b'T\x1d\xb3\xfe\xb6q\xef\xbf\x7f\xcaj\xcbZ\x84\x1ee'
```

Create `config.py`.

```sh
mkdir -p ./instance && touch ./instance/config.py
nano ./instance/config.py
```

Add the following content.

```python
# ./instance/config.py
SECRET_KEY = b'T\x1d\xb3\xfe\xb6q\xef\xbf\x7f\xcaj\xcbZ\x84\x1ee'
IP_WHITE_LIST = 'xxx.xxx.xxx.xxx'
URL_DEEPIR = 'http://xxx.xxx.xxx.xxx:3008/infer_ir'
URL_NSHIFTDB = 'https://nmrshiftdb.nmr.uni-koeln.de/NmrshiftdbServlet/nmrshiftdbaction/quickcheck'
LOGS_FILE = './instance/logging.log' #location of logs file
MAX_ZIP_SIZE = 100 #maximum size of a zip file in MB to prevent zip bomb, default is 100 MB
```

### 1.5. Start server

Using only "one" of following commands.

```sh
# run on the production server
gunicorn -w 4 -b 0.0.0.0:3007 server:app --daemon
```

```sh
# for local development only
export FLASK_APP=chem_spectra && export FLASK_DEBUG=true && flask run --host=0.0.0.0 --port=3007
```

### 1.6 Quick test

You should receive `pong` when executing the following command from another machine.

```sh
curl xxx.xxx.xxx.xxx:3007/ping
```

## 2. Usage

Send a post request with an attachment in `.dx` or `.jdx` extension.

This web service will return a zip file containing an image and a modified jcamp file with peak tables.

```
POST xxx.xxx.xxx.xxx:3007/peak_zip_jcamp_n_img
body = { file: target.jdx }
```

### 2.1. Logging

By default, the logging file is at `./instance/logging.log`.

### 2.2. Implement a new log message

Import `logging` package at where you want to write your log message

```
import logging
```

Write log message as

```
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR) //levels: DEBUG, INFO, ERROR, WARNING, CRITICAL
logger.error('message to log')
```

Note: You need to use function as the same as your logger level, which named as lowercased of level's name, to write your log message to the logs file

### 2.3. Automatic startup in crontab

To make sure ChemSpectra is started on reboot you can use this BASH script in your root crontab (if required, adapt Chemotion ELN username and home directory):

```sh
#!/bin/bash
cd /path/to/chem-spectra-app
source .venv/bin/activate
docker run --detach --name msconvert_docker \
    --rm -it \
    -e WINEDEBUG=-all \
    -v ./chem_spectra/tmp:/data proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses \
    bash
gunicorn -w 4 -b 0.0.0.0:3007 server:app --daemon
```

## 3. Run test

```sh
python -m pytest
```

### 3.1. CI configuration (maintainers)

The `test-datasets` submodule points to
[ComPlat/chem-spectra-test-files](https://github.com/ComPlat/chem-spectra-test-files),
which is currently **private**. GitHub Actions cannot clone it without credentials.

Choose one fix:

1. **Recommended for open test data:** make `chem-spectra-test-files` public.
2. **Repository secret:** add `CHEMSPECTRA_TEST_FILES_PAT` on `chem-spectra-app`
   (Settings → Secrets and variables → Actions). Use a fine-grained PAT or classic
   PAT with read access to `ComPlat/chem-spectra-test-files`.
3. **Organization setting:** in ComPlat org settings, allow GitHub Actions
   workflows to access other private repositories (then the default
   `GITHUB_TOKEN` may be sufficient).

The workflow validates that `test-datasets/catalog/datasets.json` exists after
checkout and prints these instructions if the submodule failed to initialize.
