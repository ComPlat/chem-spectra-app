"""The application factory's four contracts, each with a probe that fails first.

None of these were covered. Every one was measured on the previous code before
being fixed, and the measurement is recorded in the test that pins it.
"""

import logging
import os

import pytest

from chem_spectra import create_app, DEFAULT_LOG_FILE


@pytest.fixture
def no_flask_env(monkeypatch):
    """Clear the ambient FLASK_* environment for tests on the production path.

    Those tests now read the developer's real environment -- which is the whole
    point of the change below -- so without this they would be exactly as
    leaky as the behaviour this branch removed from the *test* path.
    """
    for key in [k for k in os.environ if k.startswith('FLASK_')]:
        monkeypatch.delenv(key)


# - - - one log handler, not one per app - - -

def test_building_several_apps_leaves_one_handler():
    """Measured before the fix: 3 handlers after 3 calls, so 3 copies of every
    record. `create_app` called `addHandler` unconditionally, and the test
    suite builds an app per fixture."""
    logger = logging.getLogger('chem_spectra')
    before = list(logger.handlers)
    try:
        for _ in range(3):
            create_app({'IP_WHITE_LIST': ''})
        installed = [h for h in logger.handlers
                     if getattr(h, '_chem_spectra_file_handler', False)]
        assert len(installed) == 1
    finally:
        for handler in list(logger.handlers):
            if handler not in before:
                logger.removeHandler(handler)
                handler.close()
        for handler in before:
            if handler not in logger.handlers:
                logger.addHandler(handler)


def test_the_package_logger_is_the_one_flask_exposes():
    """Why one handler covers everything: `app.name` is this package, so
    `app.logger` *is* `logging.getLogger('chem_spectra')`, and every
    `chem_spectra.*` module logger propagates up to it.

    If this ever stops holding, `_configure_logging` is configuring a logger
    that no longer sees the application's own records."""
    app = create_app({'IP_WHITE_LIST': ''})
    assert app.name == 'chem_spectra'
    assert app.logger is logging.getLogger('chem_spectra')


# - - - no shipped secret - - -

def test_no_secret_key_is_shipped(no_flask_env, instance_config):
    """The factory used to set `SECRET_KEY='dev'`, copied from the Flask
    tutorial. Nothing in this app uses sessions, flashing or signed cookies --
    grep finds no other mention of it -- so the constant bought nothing and
    would have signed cookies with a publicly known value the moment someone
    added a session.

    Absent, Flask raises at the point of use with a message that says exactly
    what is missing, which is the loud failure the constant was hiding.

    Production is unaffected: the deployed `instance/config.py` sets a real
    random SECRET_KEY of its own. `instance_config` here writes a file that
    does not, so this asserts the factory's own default and not a local one."""
    instance_config("IP_WHITE_LIST = ''\n")
    app = create_app()
    assert app.config['SECRET_KEY'] is None


# - - - the environment overrides the deployment file - - -

@pytest.fixture
def instance_config(request):
    """Write instance/config.py for one test, refusing to clobber a real one."""
    app = create_app({'IP_WHITE_LIST': ''})
    os.makedirs(app.instance_path, exist_ok=True)
    path = os.path.join(app.instance_path, 'config.py')
    if os.path.exists(path):
        pytest.skip('instance/config.py exists; not overwriting a real config')
    def write(body):
        with open(path, 'w') as handle:
            handle.write(body)
        return path
    request.addfinalizer(lambda: os.path.exists(path) and os.remove(path))
    return write


def test_the_environment_overrides_the_instance_file(instance_config, monkeypatch):
    """The order used to be the other way round, and it mattered here: the
    production image bakes instance/config.py in at build time, so a baked
    value could not be overridden per deployment. Measured on the old order --
    with both set, the file won."""
    instance_config("IP_WHITE_LIST = 'from-file'\n")
    monkeypatch.setenv('FLASK_IP_WHITE_LIST', 'from-env')
    assert create_app().config['IP_WHITE_LIST'] == 'from-env'


def test_the_instance_file_still_applies_when_the_environment_is_silent(
        instance_config, no_flask_env):
    """Reversing the order must not stop the file being read at all."""
    instance_config("IP_WHITE_LIST = 'from-file'\n")
    assert create_app().config['IP_WHITE_LIST'] == 'from-file'


def test_a_test_configuration_ignores_the_environment(monkeypatch):
    """An explicit test config wins outright, and the environment is not
    consulted: a stray FLASK_* variable in a developer's shell must not change
    what the suite tests. Previously the environment was read first, so it
    leaked into every key the test config did not set."""
    monkeypatch.setenv('FLASK_LOGS_FILE', '/nonexistent/should-not-be-used.log')
    app = create_app({'IP_WHITE_LIST': '127.0.0.1'})
    assert 'LOGS_FILE' not in app.config
    assert app.config['IP_WHITE_LIST'] == '127.0.0.1'


# - - - the instance directory - - -

def test_the_instance_directory_is_created_and_errors_are_not_swallowed(monkeypatch):
    """`os.makedirs` in a bare `try/except OSError` also swallowed permission
    failures, leaving the app running with no writable instance path and the
    file handler about to fail. `exist_ok=True` handles the one case that was
    meant to be tolerated and lets the rest surface."""
    app = create_app({'IP_WHITE_LIST': ''})
    assert os.path.isdir(app.instance_path)

    def refuse(*args, **kwargs):
        raise PermissionError('read-only filesystem')

    monkeypatch.setattr(os, 'makedirs', refuse)
    with pytest.raises(PermissionError):
        create_app({'IP_WHITE_LIST': ''})


def test_the_default_log_file_matches_the_instance_directory():
    """The default is CWD-relative while instance_path is absolute; in the
    image WORKDIR is /app and instance_path is /app/instance, so they agree.
    They have to, or the handler writes somewhere the factory never created."""
    app = create_app({'IP_WHITE_LIST': ''})
    assert os.path.abspath(DEFAULT_LOG_FILE) == os.path.join(
        app.instance_path, 'logging.log')
