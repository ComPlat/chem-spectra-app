from chem_spectra.model.inferencer import InferencerModel as InferModel
from chem_spectra.model.molecule import MoleculeModel


def decorate_sim_property(jbcv, molfile, isSimulateNRM=False):
    if molfile is None:
        return jbcv
    if jbcv.ncl in ('1H', '13C') and molfile.name:
        if (not jbcv.simu_peaks) or isSimulateNRM:
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
    return jbcv
