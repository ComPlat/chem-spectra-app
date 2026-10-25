from chem_spectra.model.inferencer import InferencerModel as InferModel
from chem_spectra.model.molecule import MoleculeModel

SIMULATED_NUCLEUS = {'1H': 'H', '13C': 'C'}


def decorate_sim_property(jbcv, molfile, isSimulateNMR=False):
    if molfile is None:
        return jbcv
    if jbcv.ncl in ('1H', '13C') and molfile.name:
        # if (not jbcv.simu_peaks) or isSimulateNMR:
        #     deco_jbcv = __simulate_nmr(jbcv, molfile)
        #     return deco_jbcv
        if isSimulateNMR:
            deco_jbcv = __simulate_nmr(jbcv, molfile)
            return deco_jbcv

    return jbcv


def __simulate_nmr(jbcv, molfile):
    layout = jbcv.ncl
    mm = None
    try:
        mm = MoleculeModel(molfile, layout, decorate=True)
    except Exception as error:  # noqa: F841
        return {'invalid_molfile': True, 'origin_jbcv': jbcv}
    jbcv.simu_peaks = InferModel.simulate_nmr(
        mm=mm,
        layout=layout,
    )
    # How many atoms of the nucleus the molecule has. NMRShiftDB predicts
    # fewer when it leaves atoms out (exchangeable NH/OH protons), and the
    # editor otherwise caps the signals a user may send at that smaller
    # count (#339). Only written with a simulation to go with it.
    jbcv.simu_atoms = (
        mm.count_atoms(SIMULATED_NUCLEUS[layout])
        if jbcv.simu_peaks and layout in SIMULATED_NUCLEUS else None
    )
    return jbcv
