# INSTALL


### 0. prepare

##### 0.1. install Anaconda.

##### 0.2. create env

```
$ conda create --name chem-spectra python=3.5
$ conda deactivate; conda activate chem-spectra
```

```
$ conda install -c rdkit rdkit
```

```
$ sudo apt-get install gcc libxrender1 libxext-dev
```

```
$ git clone git@git.scc.kit.edu:qj9692/chem-spectra-app.git
$ cd chem-spectra-app
$ python setup.py install # TBD
$ pip install numpy
```

##### 0.3 install nmrglue

```
$ pip uninstall nmrglue

$ cd ..
$ git clone git@bitbucket.org:ioc-general/nmrglue.git
$ cd nmrglue
$ git checkout show-all-data
$ pip install -e .
```

##### 0.4 docker msconvert

```
$ cd ../chem-spectra-app
$ mkdir chem_spectra/tmp
$ sudo chmod -R 755 chem_spectra/tmp

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

$ mkdir -p ./instance && touch ./instance/config.py
```

```python
# ./instance/config.py

SECRET_KEY = b'T\x1d\xb3\xfe\xb6q\xef\xbf\x7f\xcaj\xcbZ\x84\x1ee'
IP_WHITE_LIST = 'xxx.xxx.xxx.xxx'
URL_DEEPIR = 'http://xxx.xxx.xxx.xxx:2512/infer_ir'
URL_NSHIFTDB = 'https://nmrshiftdb.nmr.uni-koeln.de/NmrshiftdbServlet/nmrshiftdbaction/quickcheck'
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

### 4. Run test

```
$ coverage run -m pytest --disable-pytest-warnings
$ coverage report
```

### 5. Linting

```
$ flake8
```


### DEBUG

##### MS

```
$ docker exec -it msconvert_docker wine msconvert --help
```
