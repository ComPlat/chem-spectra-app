# import json


response_predict_ir = {
    "outline": {
        "code": 200,
        "text": "Load from files.",
    },
    "output": {
        "result": [
            {
                'type': 'ir',
                'fgs': [
                    {
                        'sma': 'c-,:C(-,:C)=O',
                        'confidence': 85.11,
                        'status': 'accept',
                    },
                    {
                        'sma': 'C-,:C(=O)-,:O-,:C',
                        'confidence': 93.2,
                        'status': 'accept',
                    },
                    {
                        'sma': 'c-,:[Cl]',
                        'confidence': 87.30,
                        'status': 'reject',
                    },
                ],
            },
        ],
    }
}


# class RequestPredictNmr:
#     def __init__(self):
#         pass

#     def json(self):
#         return json.dumps(request_predict_ir)


class ResponsePredictIr:
    def __init__(self):
        pass

    def json(self):
        return response_predict_ir
