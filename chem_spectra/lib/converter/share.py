import json


def parse_params(params):
    default_itg = { 'stack': [], 'refArea': 1, 'refFactor': 1, 'shift': 0 }
    default_mpy = { 'stack': [], 'smExtext': False, 'shift': 0 }
    default_wavelength = { 'name': 'CuKalpha', 'value': 0.15406, 'label': 'Cu K-alpha', 'unit': 'nm'}
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
            'waveLength': default_wavelength,
            'list_max_min_peaks': None,
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
    waveLength = params.get('waveLength')
    waveLength = json.loads(waveLength) if waveLength else default_wavelength
    listMaxMinPeaks = params.get('list_max_min_peaks')
    listMaxMinPeaks = json.loads(listMaxMinPeaks) if listMaxMinPeaks else None
    

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
        'waveLength': waveLength,
        'list_max_min_peaks': listMaxMinPeaks,
    }


def parse_solvent(base):
    if base.ncl in ['1H', '13C']:
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

        if base.ncl == '13C':
            if 'acetone' in orig_solv:
                base.dic['$CSSOLVENTNAME'] = ['Acetone-d6 (sep)']
                base.dic['$CSSOLVENTVALUE'] = ['29.920']
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(27.0, 33.0), (203.7, 209.7)]
            elif 'dmso' in orig_solv:
                peak = 39.51
                delta = 3
                base.dic['$CSSOLVENTNAME'] = ['DMSO-d6']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'methanol-d4' in orig_solv or 'meod' in orig_solv:
                peak = 49.15
                delta = 5
                base.dic['$CSSOLVENTNAME'] = ['Methanol-d4 (sep)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'dichloromethane-d2' in orig_solv:
                peak = 54.0
                delta = 3
                base.dic['$CSSOLVENTNAME'] = ['Dichloromethane-d2 (quin)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'acetonitrile-d3' in orig_solv:
                peak = 1.39
                delta = 3
                base.dic['$CSSOLVENTNAME'] = ['Acetonitrile-d3 (sep)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'benzene' in orig_solv:
                peak = 128.39
                delta = 3
                base.dic['$CSSOLVENTNAME'] = ['Benzene (t)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'chloroform-d' in orig_solv or 'cdcl3' in orig_solv:
                peak = 77.0
                delta = 3
                base.dic['$CSSOLVENTNAME'] = ['Chloroform-d (t)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
        elif base.ncl == '1H':
            if 'acetonitrile' in orig_solv:
                peak = 1.94
                delta = 0.05
                base.dic['$CSSOLVENTNAME'] = ['Acetonitrile-d3 (quin)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'acetone' in orig_solv:
                peak = 2.05
                delta = 0.05
                base.dic['$CSSOLVENTNAME'] = ['Acetone-d6 (quin)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'benzene' in orig_solv:
                peak = 7.16
                delta = 0.01
                base.dic['$CSSOLVENTNAME'] = ['Benzene (s)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'deuterium' in orig_solv:
                peak = 4.75
                delta = 0.01
                base.dic['$CSSOLVENTNAME'] = ['Deuterium oxide (s)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'dichloromethane' in orig_solv:
                peak = 5.32
                delta = 0.01
                base.dic['$CSSOLVENTNAME'] = ['Dichloromethane-d2 (t)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'dmso' in orig_solv:
                peak = 2.50
                delta = 0.02
                base.dic['$CSSOLVENTNAME'] = ['DMSO-d6 (quin)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'chloroform-d' in orig_solv or 'cdcl3' in orig_solv:
                peak = 7.27
                delta = 0.01
                base.dic['$CSSOLVENTNAME'] = ['Chloroform-d (s)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]


def reduce_pts(xys):
    num_pts_limit = 4000
    filter_ratio = 0.001
    if xys.shape[0] == 0:
        return xys
    filter_y = filter_ratio * xys[:, 1].max()
    while True:
        if xys.shape[0] < num_pts_limit:
            break
        xys = xys[xys[:, 1] > filter_y]
        filter_y = filter_y * 2
    return xys
