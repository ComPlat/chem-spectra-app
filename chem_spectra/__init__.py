import logging
import os

from flask import Flask, jsonify

# Marks the handler this module installs, so a second create_app() replaces it
# instead of stacking another one beside it. See _configure_logging.
_OURS = '_chem_spectra_file_handler'

DEFAULT_LOG_FILE = './instance/logging.log'


def _configure_logging(app):
    """Attach exactly one file handler, however often create_app() is called.

    The logger is `logging.getLogger('chem_spectra')` -- which is the very
    object Flask exposes as `app.logger`, since `app.name` is this package's
    name. Every `chem_spectra.*` module logger propagates up to it, which is
    why one handler here covers the whole package.

    It is addressed by name rather than as `app.logger` on purpose: reading
    that property runs Flask's `create_logger`, which installs its own stderr
    handler when the logger has none yet. That would be a second sink nobody
    asked for on this branch.

    This used to call `addHandler` unconditionally, so a process that built
    several apps -- the test suite builds one per fixture -- ended up with one
    handler per call and wrote every record that many times over. Measured
    before the fix: 3 handlers after 3 calls.
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    for stale in [h for h in logger.handlers if getattr(h, _OURS, False)]:
        logger.removeHandler(stale)
        stale.close()

    handler = logging.FileHandler(app.config.get('LOGS_FILE') or DEFAULT_LOG_FILE)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    setattr(handler, _OURS, True)
    logger.addHandler(handler)


def create_app(test_config=None):
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        IP_WHITE_LIST=''
    )

    if test_config is None:
        # Deployment config first, environment second, so the environment can
        # override it. The production image bakes instance/config.py in at
        # build time (Dockerfile.p2d ADDs it from the payload server), so the
        # environment is the only per-deployment knob there is -- and with the
        # order reversed, as it was, a baked value could not be overridden at
        # all. Measured: with both set, the file won.
        app.config.from_pyfile('config.py', silent=True)
        app.config.from_prefixed_env()
    else:
        # An explicit test configuration wins outright, and the environment is
        # deliberately not consulted: a stray FLASK_* variable in a developer's
        # shell should not change what the suite tests.
        app.config.from_mapping(test_config)

    os.makedirs(app.instance_path, exist_ok=True)

    _configure_logging(app)

    # ping api
    @app.route('/ping')
    def ping():
        return 'pong'

    # file api
    from chem_spectra.controller.file_api import file_api
    app.register_blueprint(file_api)

    # inference api
    from chem_spectra.controller.inference_api import infer_api
    app.register_blueprint(infer_api)

    # transform api
    from chem_spectra.controller.transform_api import trans_api
    app.register_blueprint(trans_api)

    # spectra layout api
    from chem_spectra.controller.spectra_layout_api import spectra_layout_api
    app.register_blueprint(spectra_layout_api)

    # A conversion the client asked for that this data cannot support is a bad
    # request, not a server fault, and the reason is worth returning -- the
    # alternative is silently producing a ruined spectrum.
    from chem_spectra.lib.converter.jcamp.technique import UnconvertibleSpectrum

    @app.errorhandler(UnconvertibleSpectrum)
    def _unconvertible_spectrum(err):
        return jsonify(error=str(err)), 422

    return app
