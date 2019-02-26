import os
import tempfile

import pytest
from chem_spectra import create_app

from fixtures.mock_predict import ResponsePredictNmr


@pytest.fixture
def app(mocker):
    mocker.patch(
        'requests.post',
        return_value=ResponsePredictNmr()
    )

    app = create_app({
        'TESTING': True,
        'IP_WHITE_LIST': '127.0.0.1'
    })

    yield app


@pytest.fixture
def client(app):
    return app.test_client()
