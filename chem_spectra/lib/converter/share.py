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
