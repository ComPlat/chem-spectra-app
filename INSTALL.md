# INSTALL

## 1. Requirements

- **Python ≥ 3.12.** CI tests 3.12 and 3.14; the `Dockerfile.p2d` image runs 3.14.
- **[uv](https://docs.astral.sh/uv/)** for the environment and
  **[just](https://just.systems/)** for the recipes (`just` lists them).
- **Docker**, for the msconvert sidecar that converts RAW mass-spectrometry files.

## 2. Local setup

```sh
git clone https://github.com/ComPlat/chem-spectra-app.git
cd chem-spectra-app
just sync        # uv sync --frozen --all-groups: the locked dependencies, dev tools included
```

All commands below run from the repository root.

### 2.1. msconvert

```sh
just msconvert        # starts `msconvert_docker`, mounting ./chem_spectra/tmp as /data
just msconvert-stop
```

Locally the app runs `docker exec … msconvert_docker wine msconvert …` for real.
Run the ProteoWizard build CI tests against: the version is pinned in
`.github/workflows/unit_test.yml`.

### 2.2. `instance/config.py` (optional)

Every key has a default; create the file only to change one.

```python
# instance/config.py
URL_DEEPIR = 'http://xxx.xxx.xxx.xxx:3008/infer_ir'
URL_NSHIFTDB = 'https://nmrshiftdb.nmr.uni-koeln.de/NmrshiftdbServlet/nmrshiftdbaction/quickcheck'
LOGS_FILE = '/path/to/logging.log'  # default: instance/logging.log
MAX_ZIP_SIZE = 200                  # MB, uncompressed; guards against zip bombs
```

Any key can also be set as a `FLASK_`-prefixed environment variable, which
overrides the file: for example `FLASK_MAX_ZIP_SIZE=500`.

## 3. Run

```sh
# development, on the port chem-spectra-client proxies to
FLASK_APP=chem_spectra FLASK_DEBUG=true uv run flask run --host=0.0.0.0 --port=3007

# production-style
just serve       # gunicorn -w 4 on :4000
```

Check that it answers:

```sh
curl localhost:3007/ping     # pong
```

## 4. Usage

Send a JCAMP-DX file (`.jdx`, `.dx`, `.jcamp`), a Bruker zip, a BagIt archive,
or a RAW / mzML file:

```sh
curl -F file=@target.jdx localhost:3007/zip_jcamp_n_img -o spectrum.zip
```

The answer is `spectrum.zip`, holding the composed JCAMP with its peak tables and
the preview image. A file the service will not process comes back as **422**
with `{"error": "<reason>"}`; see `docs/api-reference.md`.

## 5. Deployment

**Production** is the Chemotion stack, deployed with
[ChemOrc](https://github.com/Chemotion/ChemOrc). Its compose file runs the
`ptrxyz/chemotion:spectra-<version>` image, built by `chemotion-bakery` and not
from this repository's Dockerfile. That image sits beside the
`ptrxyz/chemotion:msconvert-<version>` sidecar, on the stack's network, with a
shared `spectra` volume.

**Developer deployments** of a branch use this repository's
`docker-compose.p2d.yml`. It builds the app from `Dockerfile.p2d` and runs it
beside the same msconvert sidecar, as a Dokploy *Docker Compose* application:
- route a domain to the `app` service on port 4000 in the Domains tab;
- set environment variables in the UI. They reach the container through
  `.env`, which the compose file loads.

**Both images** wire msconvert the same way:
- `/bin/docker` is a shim that POSTs the msconvert command to the sidecar
  (`MSC_HOST=msconvert`, `MSC_PORT=4000`);
- the two containers exchange files through the shared volume, mounted at
  `/shared`.

`just image` builds the developer image locally.

## 6. Tests and checks

```sh
just test        # pytest
just check       # ruff, mypy and pytest: what CI runs
just cov         # coverage report
```

The RAW conversion test needs `just msconvert` running.

## 7. Logging

Logs go to `instance/logging.log`, or to `LOGS_FILE` when set. Log through the
package logger:

```python
import logging

logger = logging.getLogger(__name__)
logger.warning('message')
```

`chem_spectra.*` module loggers propagate to the one file handler the app
factory installs.
