import re
import openbabel as ob


class MoleculeModel:
    def __init__(self, molfile):
        self.molfile = molfile
        self.obconv = ob.OBConversion()
        self.obmol = ob.OBMol()
        self.__load_to_obmol()
        self.can = self.__set_can()
        self.mass = self.__set_mass()


    def __load_to_obmol(self):
        self.obconv.SetInAndOutFormats('mol', 'can')
        self.obconv.ReadString(self.obmol, self.molfile)


    def __set_can(self):
        self.can = self.obconv.WriteString(self.obmol)
        self.can = re.sub('\s+', '', self.can)
        return self.can


    def __set_mass(self):
        self.mass = self.obmol.GetExactMass()
        return self.mass
