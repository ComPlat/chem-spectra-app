from flask import (
    current_app, g
)


def get_ip_white_list():
    if 'ip_white_list' not in g:
        g.ip_white_list = current_app.config['IP_WHITE_LIST'].split(';')

    return g.ip_white_list
