# chem-spectra-app

This backend web service provides NMR/IR/MS processing for jcamp/RAW/mzML files.

The frontend is provided by chem-spectra-client.

### Install

see [INSTALL.md][INSTALL]


##### Run test

```
$ coverage run -m pytest --disable-pytest-warnings
$ coverage report

$ coverage run -m pytest --disable-pytest-warnings -k ./tests/test_spectra_im.py -k 'test_meta_1H' -vv
```

### Demo & Manual

[demo & step-by-step manual](https://github.com/ComPlat/react-spectra-editor/blob/master/DEMO_MANUAL.md)





[LICENSE]: LICENSE
[INSTALL]: INSTALL.md
