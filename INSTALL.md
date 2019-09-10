# INSTALL

This guidance is tested on Linux Ubuntu 18.04.

## 1. Installation

### 1.1. Install Anaconda & Docker

Please refer to `https://www.anaconda.com/` & `https://docs.docker.com/install/`.

You can follow [unofficial Anaconda & Docker installations for Dummies](INSTALL_BASIC.md) for a test installation.

However, it is highly recommended to refer to official websites.

### 1.2. Create env

_Logout & login to load installations._

```
$ conda create --name chem-spectra python=3.5
$ source activate chem-spectra
```

```
$ conda install -c rdkit rdkit
```

```
$ sudo apt-get install gcc libxrender1 libxext-dev
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
```

### 1.5. Start server

Using only "one" of following commands.

```
# run on the production server
$ gunicorn -w 4 -b 0.0.0.0:3007 server:app --daemon
```


```
# for local development only
$ export FLASK_APP=chem_spectra && export FLASK_ENV=development && flask run --host=0.0.0.0 --port=3007
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
