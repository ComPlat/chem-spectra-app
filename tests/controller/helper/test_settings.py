from array import array
import pytest
from chem_spectra import create_app
from chem_spectra.controller.helper.settings import get_ip_white_list


@pytest.fixture
def app():
    return create_app()


def test_get_white_list_ip(app):
    with app.app_context():
        trusted_servers = get_ip_white_list()
        assert type(trusted_servers) is list
        assert len(trusted_servers) >= 0
