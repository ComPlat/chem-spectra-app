# INSTALL

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
$ pip install waitress

$ waitress-serve --port=2412 --call 'chem_spectra:create_app'
```

### 3. Usage

Send a post request with an attachment in `.dx` or `.jdx` extension.

This web service will return a zip file containing an image and a modified jcamp file with peak tables.

```
POST xxx.xxx.xxx.xxx:2412/peak_zip_jcamp_n_img
body = { file: target.jdx }
```
