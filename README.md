# ChemSpectra App

This backend web service processes spectra -- NMR, IR, UV/VIS, MS, LC/MS, XRD, CV and more -- from JCAMP-DX, Bruker, BagIt, RAW and mzML files.

![GitHub release (release name instead of tag name)](https://img.shields.io/github/v/release/ComPlat/chem-spectra-app?include_prereleases&label=version)
![Python 3.12 | 3.14](https://img.shields.io/badge/python-3.12%20%7C%203.14-blue.svg)
![Testing](https://github.com/ComPlat/chem-spectra-app/actions/workflows/unit_test.yml/badge.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

The frontend is provided by [chem-spectra-client](https://github.com/ComPlat/chem-spectra-client).

## Documentation

1. [Installation](INSTALL.md)
2. [Developer onboarding](docs/onboarding.md)
3. [Backend architecture](docs/architecture.md)
4. [API reference](docs/api-reference.md)
5. [Core runtime flows](docs/core-flows.md)
6. [Diagram maintenance](docs/dev/diagrams.md)
7. [Demo & step-by-step manual](https://github.com/ComPlat/react-spectra-editor/blob/master/DEMO_MANUAL.md)

### Run tests

```bash
just sync        # once: the locked environment, via uv
just check       # ruff, mypy and pytest -- what CI runs
just cov         # coverage report
uv run python -m pytest tests/test_spectra_im.py -k test_meta_1H -vv
```

`just` lists every recipe; `INSTALL.md` has the rest.

## Acknowledgments

<a href="https://www.dfg.de/en/"><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/8/86/DFG-logo-blau.svg/330px-DFG-logo-blau.svg.png" alt="DFG Logo" height="60"></a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href="https://www.helmholtz.de/en/"><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/7/70/Helmholtz-Gemeinschaft_Logo_2019.svg/330px-Helmholtz-Gemeinschaft_Logo_2019.svg.png" alt="Helmholtz Logo" height="60"></a>

Funded by the [Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)](https://www.dfg.de/) under the [National Research Data Infrastructure – NFDI4Chem](https://nfdi4chem.de/) – Projektnummer **441958208** since 2020.
