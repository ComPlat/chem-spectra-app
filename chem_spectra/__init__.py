import os

from flask import Flask


def create_app(test_config=None):
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        IP_WHITE_LIST=''
    )

    if test_config is None:
        app.config.from_pyfile('config.py', silent=True)
    else:
        app.config.from_mapping(test_config)

    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

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

    return app
