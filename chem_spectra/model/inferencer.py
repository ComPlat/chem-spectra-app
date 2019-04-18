import requests

url_nmrshiftdb = 'http://nmrshiftdb.nmr.uni-koeln.de/NmrshiftdbServlet/nmrshiftdbaction/quickcheck'

headers = {
    'Content-Type': 'application/json'
}


class InferencerModel:
    def __init__(self, layout, molfile, peaks, shift):
        self.layout = layout
        self.molfile = molfile
        self.peaks = peaks
        self.shift = shift


    def predict_by_peaks(self):
        peak_xs = self.__extract_x()
        solvent = self.shift.get('ref', {}) .get('nsdb')

        if self.layout == '1H':
            typ = 'nmr;1H;1d'
            data = self.__build_data(typ, peak_xs, solvent)
            rsp = requests.post(url_nmrshiftdb, headers=headers, json=data)
            return rsp
        elif self.layout == '13C':
            typ = 'nmr;13C;1d'
            data = self.__build_data(typ, peak_xs, solvent)
            rsp = requests.post(url_nmrshiftdb, headers=headers, json=data)
            return rsp

        return False


    def __extract_x(self):
        total = []
        for p in self.peaks:
            total.append(str(p['x']))

        return ';'.join(total)


    def __build_data(self, typ, peak_xs, solvent):
        return {
            'inputs':[
                {
                    'id': 1,
                    'type': typ,
                    'shifts': peak_xs,
                    'solvent': solvent,
                },
            ],
            'moltxt': self.molfile
        }
