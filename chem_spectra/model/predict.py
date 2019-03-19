import requests

url_nmrshiftdb = 'http://nmrshiftdb.nmr.uni-koeln.de/NmrshiftdbServlet/nmrshiftdbaction/quickcheck'

headers = {
    'Content-Type': 'application/json'
}


def build_data(typ, peak_xs, molfile, solvent):
    return {
        'inputs':[
            {
                'id': 1,
                'type': typ,
                'shifts': peak_xs,
                'solvent': solvent,
            },
        ],
        'moltxt': molfile
    }


def extract_x(peaks):
    total = []
    for p in peaks:
        total.append(str(p['x']))

    return ';'.join(total)


def predict_by_peaks(layout, molfile, peaks, shift):
    peak_xs = extract_x(peaks)
    solvent = shift.get('ref', {}) .get('nsdb')

    if layout == '1H':
        typ = 'nmr;1H;1d'
        data = build_data(typ, peak_xs, molfile, solvent)
        rsp = requests.post(url_nmrshiftdb, headers=headers, json=data)
        return rsp
    elif layout == '13C':
        typ = 'nmr;13C;1d'
        data = build_data(typ, peak_xs, molfile, solvent)
        rsp = requests.post(url_nmrshiftdb, headers=headers, json=data)
        return rsp

    return False
