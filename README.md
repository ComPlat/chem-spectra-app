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


### Architecture
To see archicture, you can view it at [ARCHITECTURE.md](./docs/ARCHITECTURE.md)


## Acknowledgments

This project has been funded by the **[DFG]**.

[![DFG Logo]][DFG]


Funded by the [Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)](https://www.dfg.de/) under the [National Research Data Infrastructure – NFDI4Chem](https://nfdi4chem.de/) – Projektnummer **441958208** since 2020.


[DFG]: https://www.dfg.de/en/
[DFG Logo]: https://www.dfg.de/zentralablage/bilder/service/logos_corporate_design/logo_negativ_267.png
[LICENSE]: LICENSE
[INSTALL]: INSTALL.md
