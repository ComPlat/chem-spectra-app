# ChemSpectra App

This backend web service provides NMR/IR/MS processing for jcamp/RAW/mzML files.

![GitHub release (release name instead of tag name)](https://img.shields.io/github/v/release/ComPlat/chem-spectra-app?include_prereleases&label=version)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)
![Testing](https://github.com/ComPlat/chem-spectra-app/actions/workflows/unit_test.yml/badge.svg)

The frontend is provided by [chem-spectra-client](https://github.com/ComPlat/chem-spectra-client).

## Documentation

1. [Installation](INSTALL.md)
2. [Demo & step-by-step manual](https://github.com/ComPlat/react-spectra-editor/blob/master/DEMO_MANUAL.md)
3. Software's [Archicture](./docs/ARCHITECTURE.md)

### Run test

```bash
coverage run -m pytest --disable-pytest-warnings
coverage report

coverage run -m pytest --disable-pytest-warnings -k ./tests/test_spectra_im.py -k 'test_meta_1H' -vv
```

## Acknowledgments

This project has been funded by the **[DFG](https://www.dfg.de/en/)**.

![DFG Logo](https://chemotion.net/img/logos/DFG_logo.png)

Funded by the [Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)](https://www.dfg.de/) under the [National Research Data Infrastructure – NFDI4Chem](https://nfdi4chem.de/) – Projektnummer **441958208** since 2020.
