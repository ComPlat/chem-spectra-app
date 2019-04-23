import json


response_predict_ir = {
   "status": True,
   "results": {
      "C-,:C1(-,:C)-,:S-,:C-,:C-,:S-,:1": {
         "valid": False,
         "accuracy": 0,
         "exist": False
      },
      "C-,:C(-,:C)=O": {
         "valid": True,
         "accuracy": 75,
         "exist": True
      }
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
        return json.dumps(response_predict_ir)
