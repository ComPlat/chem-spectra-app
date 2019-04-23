# Chem_spectra

This web service provides automatically peaks-picking for jcamp files.

### Funding

This project has been funded by the  ![DFG](http://www.dfg.de/includes/images/dfg_logo.gif)


### License

see [LICENSE][LICENSE]


### Install

see [INSTALL.md][INSTALL]


### Run development

```
$ export FLASK_ENV=development
$ export FLASK_APP=chem_spectra
$ flask run --host=0.0.0.0 --port=2412
```

### Run Production

```
$ waitress-serve --port=2412 --call 'chem_spectra:create_app'
```


### Run test

```
$ coverage run -m pytest --disable-pytest-warnings
$ coverage report
```




[LICENSE]: LICENSE
[INSTALL]: INSTALL.md
