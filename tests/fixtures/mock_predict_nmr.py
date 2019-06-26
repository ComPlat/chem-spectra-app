import json


request_predict_nmr = {
    "layout": "13C",
    "peaks": [
        {
            "x": 117,
            "y": 117.5,
        },
        {
            "x": 123,
            "y": 123.5,
        },
        {
            "x": 139,
            "y": 139.5,
        }
    ],
    "shift": {
        "ref": {
            "name": "Dichloromethane-d2 (quin)",
            "value": 54,
            "label": "CD$2Cl$2",
        },
        "peak": "false",
        "enable": "true",
    },
    "molfile": "62-53-3\n\n\n  7  7  0  0  0  0  0  0  0  0999 V2000\n   -0.1116    1.3830    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8260    0.9705    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8260    0.1454    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.1116   -0.2670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6028    0.1454    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6028    0.9705    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.1116   -1.0920    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0  0  0  0\n  2  3  1  0  0  0  0\n  3  4  2  0  0  0  0\n  4  5  1  0  0  0  0\n  5  6  2  0  0  0  0\n  6  1  1  0  0  0  0\n  4  7  1  0  0  0  0\nM  END"  # noqa:
}


response_predict_nmr = {
    "result": [
        {
            "id": 1,
            "shifts": [
                {
                    "atom": 1,
                    "diff": 22.712247548698286,
                    "prediction": 11.667065438606755,
                    "real": -11.045182110091531,
                    "status": "reject",
                },
                {
                    "atom": 2,
                    "diff": 30.449941791740933,
                    "prediction": 19.722594964537528,
                    "real": -10.727346827203405,
                    "status": "reject",
                },
            ],
            "statistics": {
                "accept": 5,
                "missing": 1,
                "reject": 6,
                "total": 15,
                "warning": 3,
            },
            "type": "nmr;13C;1d",
        }
    ]
}


class RequestPredictNmr:
    def __init__(self):
        pass

    def json(self):
        return json.dumps(request_predict_nmr)


class ResponsePredictNmr:
    def __init__(self):
        pass

    def json(self):
        return response_predict_nmr
