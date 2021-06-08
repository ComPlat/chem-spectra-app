from chem_spectra.model.inferencer import InferencerModel as InferModel
from chem_spectra.model.molecule import MoleculeModel


def decorate_sim_property(jbcv, molfile, isSimulateNRM=False):
    # if jbcv.ncl in ('1H', '13C') and not jbcv.simu_peaks and molfile.name:
    #     deco_jbcv = __simulate_nmr(jbcv, molfile)
    #     return deco_jbcv
    if jbcv.ncl in ('1H', '13C') and molfile.name:
        if (not jbcv.simu_peaks) or isSimulateNRM:
          deco_jbcv = __simulate_nmr(jbcv, molfile)
          return deco_jbcv

    return jbcv


def __simulate_nmr(jbcv, molfile):
    layout = jbcv.ncl
    mm = MoleculeModel(molfile, layout, decorate=True)
    jbcv.simu_peaks = InferModel.simulate_nmr(
        mm=mm,
        layout=layout,
    )
    return jbcv
