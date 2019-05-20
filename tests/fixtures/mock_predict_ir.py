# import json


response_predict_ir = {
    "outline": {
        "code": 200,
        "text": "Load from files.",
    },
    "output": {
        "result": [
            {
                "c-,:[N&+](=O)-,:[O&-]": {
                    "confidence": 99.19,
                    "status": "accept",
                },
                "c-,:C(=O)-,:O-,:C": {
                    "confidence": 95.2,
                    "status": "accept",
                },
                "c-,:[Br]": {
                    "confidence": 88.46,
                    "status": "accept"
                }
            }
        ]
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
