import json


def parse_params(params):
    default_itg = {'stack': [], 'refArea': 1, 'refFactor': 1, 'shift': 0}
    default_mpy = {'stack': [], 'smExtext': False, 'shift': 0}
    default_wavelength = {'name': 'CuKalpha', 'value': 0.15406, 'label': 'Cu K-alpha', 'unit': 'nm'}
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
            'cyclicvolta': None,
            'jcamp_idx': 0,
            'axesUnits': None,
            'detector': None,
            'dsc_meta_data': None,
        }

    select_x = params.get('select_x', None)
    ref_name = params.get('ref_name', None)
    ref_value = params.get('ref_value', None)
    peaks_str = params.get('peaks_str', None)
    delta = 0.0
    mass = params.get('mass', 0)
    mass = mass if mass else 0
    scan = params.get('scan', None)
    thres = params.get('thres', None)
    clear = params.get('clear', False)
    clear = clear if clear else False
    integration = params.get('integration')
    integration = json.loads(integration) if integration else default_itg
    multiplicity = params.get('multiplicity')
    multiplicity = json.loads(multiplicity) if multiplicity else default_mpy
    ext = params.get('ext', '')
    ext = ext if ext else ''
    fname = params.get('fname', '').split('.')
    fname = fname[:-2] if (len(fname) > 2 and (fname[-2] in ['edit', 'peak'])) else fname[:-1]
    fname = '.'.join(fname)
    waveLength = params.get('waveLength')
    waveLength = json.loads(waveLength) if waveLength else default_wavelength

    jcamp_idx = params.get('jcamp_idx', 0)
    jcamp_idx = jcamp_idx if jcamp_idx else 0
    axesUnitsJson = params.get('axesUnits')
    axesUnitsDic = json.loads(axesUnitsJson) if axesUnitsJson else None
    axesUnits = None
    if axesUnitsDic != None and 'axes' in axesUnitsDic:
        axes = axesUnitsDic.get('axes', [{'xUnit': '', 'yUnit': ''}])
        try:
            axesUnits = axes[jcamp_idx]
        except:
            pass

    cyclicvolta = params.get('cyclic_volta')
    cyclicvolta = json.loads(cyclicvolta) if cyclicvolta else None
    listMaxMinPeaks = None
    user_data_type_mapping = params.get('data_type_mapping')
    detector = params.get('detector')
    detector = json.loads(detector) if detector else None
    dsc_meta_data = params.get('dsc_meta_data')
    dsc_meta_data = json.loads(dsc_meta_data) if dsc_meta_data else None
    if (cyclicvolta is not None):
        spectraList = cyclicvolta['spectraList']
        if (len(spectraList) > 0):
            spectra = spectraList[jcamp_idx]
            listMaxMinPeaks = spectra['list']

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
        'cyclicvolta': cyclicvolta,
        'jcamp_idx': jcamp_idx,
        'axesUnits': axesUnits,
        'user_data_type_mapping': user_data_type_mapping,
        'detector': detector,
        'dsc_meta_data': dsc_meta_data,
    }


def parse_solvent(base):
    if base.ncl in ['1H', '13C']:
        ref_name = (
            base.params['ref_name'] or
            base.dic.get('$CSSOLVENTNAME', [''])[0]
        )
        # if ref_name and ref_name != '- - -':
        if ref_name:  # skip when the solvent is exist.
            return

        sn = base.dic.get('.SOLVENTNAME', [''])
        sr = base.dic.get('.SHIFTREFERENCE', [''])
        sn = sn if isinstance(sn, str) else sn[0]
        sr = sr if isinstance(sr, str) else sr[0]
        orig_solv = (sn + sr).lower()

        if base.ncl == '13C':
            if 'acetone' in orig_solv:
                base.dic['$CSSOLVENTNAME'] = ['Acetone-d6 (sep)']
                base.dic['$CSSOLVENTVALUE'] = ['29.640']
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(27.0, 33.0), (203.7, 209.7)]
            elif 'dmso' in orig_solv:
                peak = 39.52
                delta = 3
                base.dic['$CSSOLVENTNAME'] = ['DMSO-d6']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'methanol-d4' in orig_solv or 'meod' in orig_solv:
                peak = 49.00
                delta = 5
                base.dic['$CSSOLVENTNAME'] = ['Methanol-d4 (sep)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'dichloromethane-d2' in orig_solv:
                peak = 53.84
                delta = 3
                base.dic['$CSSOLVENTNAME'] = ['Dichloromethane-d2 (quin)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'acetonitrile-d3' in orig_solv:
                peak = 1.32
                delta = 3
                base.dic['$CSSOLVENTNAME'] = ['Acetonitrile-d3 (sep)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'benzene' in orig_solv:
                peak = 128.06
                delta = 3
                base.dic['$CSSOLVENTNAME'] = ['Benzene (t)']
                base.dic['$CSSOLVENTVALUE'] = [str(peak)]
                base.dic['$CSSOLVENTX'] = ['0']
                base.solv_peaks = [(peak - delta, peak + delta)]
            elif 'chloroform-d' in orig_solv or 'cdcl3' in orig_solv:
                peak = 77.16
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
                peak = 4.79
                delta = 0.01
                base.dic['$CSSOLVENTNAME'] = ['D2O (s)']
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
                peak = 7.26
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
