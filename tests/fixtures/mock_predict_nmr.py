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
            "type": "nmr;13C;1d",
            "statistics": {
                "accept": 0,
                "warning": 0,
                "reject": 5,
                "missing": 1,
                "total": 6
            },
            "shifts": [
                {
                   "atom": 1,
                   "prediction": 135.5500030517578,
                   "real": 117.0,
                   "diff": 18.550003051757812,
                   "status": "reject"
                },
                {
                   "atom": 2,
                   "prediction": 139.1999969482422,
                   "real": 123.0,
                   "diff": 16.199996948242188,
                   "status": "reject"
                },
                {
                   "atom": 6,
                   "prediction": 139.1999969482422,
                   "real": 123.0,
                   "diff": 16.199996948242188,
                   "status": "reject"
                },
                {
                   "atom": 3,
                   "prediction": 121.5999984741211,
                   "real": 103.0,
                   "diff": 18.599998474121094,
                   "status": "reject"
                },
                {
                   "atom": 5,
                   "prediction": 121.5999984741211,
                   "real": 103.0,
                   "diff": 18.599998474121094,
                   "status": "reject"
                },
                {
                   "atom": 4,
                   "prediction": 147.13500213623047,
                   "real": 0.0,
                   "diff": 0.0,
                   "status": "missing"
                }
            ]
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
        return json.dumps(response_predict_nmr)
