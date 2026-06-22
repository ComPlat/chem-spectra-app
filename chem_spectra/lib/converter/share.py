import json


def _is_placeholder_solvent_name(solvent_name):
    if solvent_name is None:
        return True
    return solvent_name.strip().upper() in ['', '- - -', 'AUTO-OFFSET']


def _solvent_config(ncl, solvent_name):
    if _is_placeholder_solvent_name(solvent_name):
        return None

    solvent_name = solvent_name.lower()

    if ncl == '13C':
        if 'acetone' in solvent_name:
            return ('Acetone-d6 (sep)', '29.640', [(27.0, 33.0), (203.7, 209.7)])
        if 'dmso' in solvent_name:
            peak = 39.52
            delta = 3
            return ('DMSO-d6', str(peak), [(peak - delta, peak + delta)])
        if 'methanol-d4' in solvent_name or 'meod' in solvent_name:
            peak = 49.00
            delta = 5
            return ('Methanol-d4 (sep)', str(peak), [(peak - delta, peak + delta)])
        if 'dichloromethane-d2' in solvent_name:
            peak = 53.84
            delta = 3
            return ('Dichloromethane-d2 (quin)', str(peak), [(peak - delta, peak + delta)])
        if 'acetonitrile-d3' in solvent_name:
            peak = 1.32
            delta = 3
            return ('Acetonitrile-d3 (sep)', str(peak), [(peak - delta, peak + delta)])
        if 'benzene' in solvent_name:
            peak = 128.06
            delta = 3
            return ('Benzene (t)', str(peak), [(peak - delta, peak + delta)])
        if 'chloroform-d' in solvent_name or 'cdcl3' in solvent_name:
            peak = 77.16
            delta = 3
            return ('Chloroform-d (t)', str(peak), [(peak - delta, peak + delta)])

    if ncl == '1H':
        if 'acetonitrile' in solvent_name:
            peak = 1.94
            delta = 0.05
            return ('Acetonitrile-d3 (quin)', str(peak), [(peak - delta, peak + delta)])
        if 'acetone' in solvent_name:
            peak = 2.05
            delta = 0.05
            return ('Acetone-d6 (quin)', str(peak), [(peak - delta, peak + delta)])
        if 'benzene' in solvent_name:
            peak = 7.16
            delta = 0.01
            return ('Benzene (s)', str(peak), [(peak - delta, peak + delta)])
        if 'deuterium' in solvent_name:
            peak = 4.79
            delta = 0.01
            return ('D2O (s)', str(peak), [(peak - delta, peak + delta)])
        if 'dichloromethane' in solvent_name:
            peak = 5.32
            delta = 0.01
            return ('Dichloromethane-d2 (t)', str(peak), [(peak - delta, peak + delta)])
        if 'dmso' in solvent_name:
            peak = 2.50
            delta = 0.02
            return ('DMSO-d6 (quin)', str(peak), [(peak - delta, peak + delta)])
        if 'chloroform-d' in solvent_name or 'cdcl3' in solvent_name:
            peak = 7.26
            delta = 0.01
            return ('Chloroform-d (s)', str(peak), [(peak - delta, peak + delta)])

    return None


def _apply_solvent_config(base, config, overwrite_meta=False):
    if config is None:
        return False

    name, peak_value, solv_peaks = config
    base.solv_peaks = solv_peaks

    existing_name = base.dic.get('$CSSOLVENTNAME', [''])[0]
    existing_value = base.dic.get('$CSSOLVENTVALUE', [''])[0]
    existing_x = base.dic.get('$CSSOLVENTX', [''])[0]

    if overwrite_meta or not existing_name:
        base.dic['$CSSOLVENTNAME'] = [name]
    if overwrite_meta or not existing_value:
        base.dic['$CSSOLVENTVALUE'] = [peak_value]
    if overwrite_meta or not existing_x:
        base.dic['$CSSOLVENTX'] = ['0']

    return True


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
        param_ref_name = base.params['ref_name']
        existing_ref_name = base.dic.get('$CSSOLVENTNAME', [''])[0]

        if _apply_solvent_config(base, _solvent_config(base.ncl, param_ref_name)):
            return
        if _apply_solvent_config(base, _solvent_config(base.ncl, existing_ref_name)):
            return
        if (
            not _is_placeholder_solvent_name(param_ref_name)
        ) or (
            not _is_placeholder_solvent_name(existing_ref_name)
        ):
            return

        sn = base.dic.get('.SOLVENTNAME', [''])
        sr = base.dic.get('.SHIFTREFERENCE', [''])
        sn = sn if isinstance(sn, str) else sn[0]
        sr = sr if isinstance(sr, str) else sr[0]
        orig_solv = (sn + sr).lower()
        _apply_solvent_config(
            base,
            _solvent_config(base.ncl, orig_solv),
            overwrite_meta=True,
        )


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
