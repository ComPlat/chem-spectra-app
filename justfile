# Task runner for chem-spectra-app. `uv` manages Python and dependencies from
# uv.lock; run `just` to list recipes.
#
# The suite needs a msconvert container for the RAW fixtures -- see `msconvert`
# below. Without it, tests/test_ms.py::test_ms_raw_converter_composer fails.

# list available recipes
default:
    @just --list

# install/refresh .venv from the lockfile, all dependency groups
sync:
    uv sync --frozen --all-groups

# re-resolve and rewrite uv.lock (after editing pyproject.toml)
lock:
    uv lock

# the test suite
test:
    uv run python -m pytest -q

# coverage over the app
cov:
    uv run python -m coverage run --source=chem_spectra -m pytest -q
    uv run python -m coverage report

# modernization lint (ruff pyupgrade rules, target py312)
lint:
    uv run ruff check chem_spectra tests

# apply auto-fixable lint findings
lint-fix:
    uv run ruff check chem_spectra tests --fix

# static type checking; the exemption list lives in pyproject.toml
mypy:
    uv run mypy chem_spectra

# ruff + mypy -- the static gate
static: lint mypy

# the local gate: lint + types + tests
check: lint mypy test

# start the msconvert container the RAW fixtures need. Pinned to the build CI
# uses (.github/workflows/unit_test.yml): `latest` moved to 4.0 on 2026-10-09,
# which reports success without writing the mzML.
msconvert:
    mkdir -p chem_spectra/tmp && chmod -R 777 chem_spectra/tmp
    docker run --detach --name msconvert_docker --rm -it \
        -e WINEDEBUG=-all \
        -v {{justfile_directory()}}/chem_spectra/tmp:/data \
        proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:3.0.26121-ed8dc8a@sha256:52d830174eec595fc9fb17628effbbacec8b931fc13620a87fb44c77263daaef bash

# stop it again
msconvert-stop:
    -docker stop msconvert_docker

# build the production image
image:
    docker build -f Dockerfile.p2d -t chem-spectra-app .

# run the app locally the way the image does
serve:
    uv run gunicorn --timeout 600 -w 4 -b 0.0.0.0:4000 server:app
