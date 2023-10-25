from flask import jsonify, Blueprint
from flask_jwt_extended import jwt_required, create_access_token, get_jwt_identity, create_refresh_token


refresh_token_api = Blueprint('refresh_token_api', __name__)

@refresh_token_api.route("/api/v1/chemspectra/login", methods=["POST"])
def login():
    access_token = create_access_token(identity="user")
    refresh_token = create_refresh_token(identity="user")
    return jsonify(access_token=access_token, refresh_token=refresh_token)


@refresh_token_api.route("/api/v1/chemspectra/refresh", methods=["POST"])
@jwt_required(refresh=True)
def refresh():
    identity = get_jwt_identity()
    access_token = create_access_token(identity=identity)
    return jsonify(access_token=access_token)