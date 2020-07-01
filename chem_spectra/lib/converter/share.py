import json


def parse_params(params):
    default_itg = { 'stack': [], 'refArea': 1, 'refFactor': 1, 'shift': 0 }
    default_mpy = { 'stack': [], 'smExtext': False, 'shift': 0 }
    if not params:
        return {
            'select_x': None,
            'ref_name': None,
            'ref_value': None,
            'peaks_str': None,
            'delta': 0.0,
            'mass': 0,
            'scan': None,
            'thres': None,
            'clear': False,
            'integration': default_itg,
            'multiplicity': default_mpy,
            'fname': '',
        }

    select_x = params.get('select_x', None)
    ref_name = params.get('ref_name', None)
    ref_value = params.get('ref_value', None)
    peaks_str = params.get('peaks_str', None)
    delta = 0.0
    mass = params.get('mass', 0)
    scan = params.get('scan', None)
    thres = params.get('thres', None)
    clear = params.get('clear', False)
    integration = params.get('integration')
    integration = json.loads(integration) if integration else default_itg
    multiplicity = params.get('multiplicity')
    multiplicity = json.loads(multiplicity) if multiplicity else default_mpy
    ext = params.get('ext', '')
    fname = params.get('fname', '').split('.')
    fname = fname[:-2] if (len(fname) > 2 and (fname[-2] in ['edit', 'peak'])) else fname[:-1]
    fname = '.'.join(fname)

    try:
        if select_x and float(select_x) != 0.0 and ref_name != '- - -':
            delta = float(ref_value) - float(select_x)
    except:  # noqa
        pass

    return {
        'select_x': select_x,
        'ref_name': ref_name,
        'ref_value': ref_value,
        'peaks_str': peaks_str,
        'delta': delta,
        'mass': mass,
        'scan': scan,
        'thres': thres,
        'clear': clear,
        'integration': integration,
        'multiplicity': multiplicity,
        'ext': ext,
        'fname': fname,
    }


def parse_solvent(base):
    if base.ncl == '13C':
        ref_name = (
            base.params['ref_name'] or
            base.dic.get('$CSSOLVENTNAME', [''])[0]
        )
        # if ref_name and ref_name != '- - -':
        if ref_name: # skip when the solvent is exist.
            return

        sn = base.dic.get('.SOLVENTNAME', [''])
        sr = base.dic.get('.SHIFTREFERENCE', [''])
        sn = sn if isinstance(sn, str) else sn[0]
        sr = sr if isinstance(sr, str) else sr[0]
        orig_solv = (sn + sr).lower()

        if 'acetone' in orig_solv:
            base.dic['$CSSOLVENTNAME'] = ['Acetone-d6 (sep)']
            base.dic['$CSSOLVENTVALUE'] = ['29.920']
            base.dic['$CSSOLVENTX'] = ['0']
            base.solv_peaks = [(27.0, 33.0), (203.7, 209.7)]
        elif 'dmso' in orig_solv:
            base.dic['$CSSOLVENTNAME'] = ['DMSO-d6']
            base.dic['$CSSOLVENTVALUE'] = ['39.51']
            base.dic['$CSSOLVENTX'] = ['0']
            base.solv_peaks = [(36.0, 43.0)]
        elif 'methanol-d4' in orig_solv or 'meod' in orig_solv:
            base.dic['$CSSOLVENTNAME'] = ['Methanol-d4 (sep)']
            base.dic['$CSSOLVENTVALUE'] = ['49.15']
            base.dic['$CSSOLVENTX'] = ['0']
            base.solv_peaks = [(44, 54)]
        elif 'dichloromethane-d2' in orig_solv:
            base.dic['$CSSOLVENTNAME'] = ['Dichloromethane-d2 (quin)']
            base.dic['$CSSOLVENTVALUE'] = ['54.0']
            base.dic['$CSSOLVENTX'] = ['0']
            base.solv_peaks = [(51.0, 57.0)]
        elif 'acetonitrile-d3' in orig_solv:
            base.dic['$CSSOLVENTNAME'] = ['Acetonitrile-d3 (sep)']
            base.dic['$CSSOLVENTVALUE'] = ['1.39']
            base.dic['$CSSOLVENTX'] = ['0']
            base.solv_peaks = [(-1.0, 4.0)]
        elif 'benzene' in orig_solv:
            base.dic['$CSSOLVENTNAME'] = ['Benzene (t)']
            base.dic['$CSSOLVENTVALUE'] = ['128.390']
            base.dic['$CSSOLVENTX'] = ['0']
            base.solv_peaks = [(125.4, 131.4)]
        elif 'chloroform-d' in orig_solv or 'cdcl3' in orig_solv:
            base.dic['$CSSOLVENTNAME'] = ['Chloroform-d (t)']
            base.dic['$CSSOLVENTVALUE'] = ['77.00']
            base.dic['$CSSOLVENTX'] = ['0']
            base.solv_peaks = [(74.0, 80.0)]
