import re

import openbabel as ob


obconv = ob.OBConversion()
obconv.SetInAndOutFormats('mol', 'can')


def molfile2chem(molfile):
    obmol = ob.OBMol()
    mfstr = molfile.stream.read().decode('utf-8')
    obconv.ReadString(obmol, mfstr)

    mass = obmol.GetExactMass()

    smi = obconv.WriteString(obmol)
    smi = re.sub('\s+', '', smi)
    smi

    return smi, mass
