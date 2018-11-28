import os
import tempfile

import pytest
from chem_spectra import create_app


@pytest.fixture
def app():
    app = create_app({
        'TESTING': True,
        'IP_WHITE_LIST': '127.0.0.1'
    })

    yield app


@pytest.fixture
def client(app):
    return app.test_client()
