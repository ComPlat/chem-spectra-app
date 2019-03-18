import requests

url_nmrshiftdb = 'http://nmrshiftdb.nmr.uni-koeln.de/NmrshiftdbServlet/nmrshiftdbaction/quickcheck'

headers = {
    'Content-Type': 'application/json'
}


def build_data(typ, shifts, molfile):
    return {
        'inputs':[
            {
                'id': 1,
                'type': typ,
                'shifts': shifts
            },
        ],
        'moltxt': molfile
    }


def extract_x(peaks):
    total = []
    for p in peaks:
        total.append(str(p['x']))

    return ';'.join(total)


def predict_by_peaks(layout, molfile, peaks):
    shifts = extract_x(peaks)

    if layout == '1H':
        typ = 'nmr;1H;1d'
        data = build_data(typ, shifts, molfile)
        rsp = requests.post(url_nmrshiftdb, headers=headers, json=data)
        return rsp
    elif layout == '13C':
        typ = 'nmr;13C;1d'
        data = build_data(typ, shifts, molfile)
        rsp = requests.post(url_nmrshiftdb, headers=headers, json=data)
        return rsp

    return False
