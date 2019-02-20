import requests

url_nmrshiftdb = 'http://nmr-sdbtest.nmr.uni-koeln.de/NmrshiftdbServlet/nmrshiftdbaction/quickcheck'

headers = {
    'Content-Type': 'application/json'
}


def build_data(typ, shifts, molecule):
    return {
        'inputs':[
            {
                'id': 1,
                'type': typ,
                'shifts': shifts
            },
        ],
        'moltxt': molecule
    }


def extract_x(peaks):
    total = []
    for p in peaks:
        total.append(str(p['x']))

    return ';'.join(total)


def predict_by_peaks(payload):
    layout = payload.get('layout')
    shifts = extract_x(payload.get('peaks'))
    molecule = payload.get('molecule')

    if not molecule:
        return False

    if layout == '1H':
        typ = 'nmr;1H;1d'
        data = build_data(typ, shifts, molecule)
        rsp = requests.post(url_nmrshiftdb, headers=headers, json=data)
        return rsp
    elif layout == '13C':
        typ = 'nmr;13C;1d'
        data = build_data(typ, shifts, molecule)
        rsp = requests.post(url_nmrshiftdb, headers=headers, json=data)
        return rsp

    return False
