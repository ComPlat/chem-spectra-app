# INSTALL

This guidance is tested on Linux Ubuntu 18.04 & Ubuntu 20.04

## 1. Installation

### 1.1. Install Anaconda & Docker

Please refer to [https://www.anaconda.com/] & [https://docs.docker.com/install/].

You can follow [unofficial Anaconda & Docker installations for Dummies](INSTALL_BASIC.md) for a test installation.

However, it is highly recommended to refer to official websites.

### 1.2. Create env

_Logout & login to load installations._

```
$ conda create --name chem-spectra python=3.8
$ source activate chem-spectra
```

```
$ sudo apt-get install gcc libxrender1 libxext-dev pkg-config g++
```

```
$ git clone https://github.com/ComPlat/chem-spectra-app.git
$ cd chem-spectra-app
$ pip install -r requirements.txt
```

The project is in `~/chem-spectra-app`.


### 1.3. Use msconvert in Docker

```
$ docker pull chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
```

__Make sure you are in the `chem-spectra-app` folder.__
```
$ cd ~/chem-spectra-app
```

```
$ mkdir chem_spectra/tmp
$ sudo chmod -R 755 chem_spectra/tmp
```
To running docker for converting Mass spectrum, run this command:
```
$ docker run --detach --name msconvert_docker \
    --rm -it \
    -e WINEDEBUG=-all \
    -v [YOUR LOCATION TO chem-spectra-app]/chem_spectra/tmp:/data \
    chambm/pwiz-skyline-i-agree-to-the-vendor-licenses bash
```

For example, if your location of `chem-spectra-app` is `/home/ubuntu/chem-spectra-app`
```
$ docker run --detach --name msconvert_docker \
    --rm -it \
    -e WINEDEBUG=-all \
    -v /home/ubuntu/chem-spectra-app/chem_spectra/tmp:/data \
    chambm/pwiz-skyline-i-agree-to-the-vendor-licenses bash
```

### 1.4. Add `config.py`

Generate a secret key.

```
$ python -c 'import os; print(os.urandom(16))'
>> b'T\x1d\xb3\xfe\xb6q\xef\xbf\x7f\xcaj\xcbZ\x84\x1ee'
```

Create `config.py`.

```
$ mkdir -p ./instance && touch ./instance/config.py
$ vim ./instance/config.py
```

Add content.

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

```
# run on the production server
$ gunicorn -w 4 -b 0.0.0.0:3007 server:app --daemon
```


```
# for local development only
$ export FLASK_APP=chem_spectra && export FLASK_DEBUG=true && flask run --host=0.0.0.0 --port=3007
```

### 1.6 Quick test

You should receive `pong` when executing the following command from another machine.

```
$ curl xxx.xxx.xxx.xxx:3007/ping
```

## 2. Usage

Send a post request with an attachment in `.dx` or `.jdx` extension.

This web service will return a zip file containing an image and a modified jcamp file with peak tables.

```
POST xxx.xxx.xxx.xxx:3007/peak_zip_jcamp_n_img
body = { file: target.jdx }
```

### 2.1. Logging
By default, the logging file is at 
```
./instance/logging.log
```

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

sudo -H -u production bash -c "cd /home/production/chem-spectra-app && \
  source /home/production/anaconda3/bin/activate chem-spectra && \
  gunicorn -w 4 -b 0.0.0.0:3007 server:app --daemon"

# Remember to modify path according to your installation
docker run --detach --name msconvert_docker \
    --rm -it \
    -e WINEDEBUG=-all \
    -v /home/production/chem-spectra-app/chem_spectra/tmp:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses \
    bash
```

## 3. Run test

```
$ coverage run -m pytest --disable-pytest-warnings
$ coverage report
```

## 4. Linting

```
$ flake8
```


## DEBUG

### msconvert_docker

```
$ docker exec -it msconvert_docker wine msconvert --help
```
