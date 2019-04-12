# INSTALL


### 0. prepare

0.1. install Anaconda.

0.2. create env

```
$ conda create --name chem-spectra python=3.5
$ source activate chem-spectra
```

```
$ conda install -c openbabel openbabel
```

```
$ git clone git@git.scc.kit.edu:qj9692/chem-spectra-app.git
$ cd chem-spectra-app
$ python setup.py install
```

0.3 install nmrglue

```
$ pip uninstall nmrglue

$ cd ..
$ git clone git@bitbucket.org:ioc-general/nmrglue.git
$ cd nmrglue
$ git co show-all-data
$ pip install -e .
$ python setup.py install
```

0.4 docker msconvert

```
$ docker pull chambm/pwiz-skyline-i-agree-to-the-vendor-licenses

$ docker run --detach --name msconvert_docker \
    --rm -it \
    -e WINEDEBUG=-all \
    -v /ABSOLUTE_PATH_TO/chem-spectra-app/chem_spectra/tmp:/data \
    chambm/pwiz-skyline-i-agree-to-the-vendor-licenses bash

```


### 1. Add `config.py`

```
$ python -c 'import os; print(os.urandom(16))'

>> b'T\x1d\xb3\xfe\xb6q\xef\xbf\x7f\xcaj\xcbZ\x84\x1ee'
```

```python
# ./instance/config.py

SECRET_KEY = b'T\x1d\xb3\xfe\xb6q\xef\xbf\x7f\xcaj\xcbZ\x84\x1ee'
IP_WHITE_LIST = 'xxx.xxx.xxx.xxx'
```

### 2. Run

```
$ gunicorn -w 4 -b 0.0.0.0:2412 server:app --daemon
```

### 3. Usage

Send a post request with an attachment in `.dx` or `.jdx` extension.

This web service will return a zip file containing an image and a modified jcamp file with peak tables.

```
POST xxx.xxx.xxx.xxx:2412/peak_zip_jcamp_n_img
body = { file: target.jdx }
```
