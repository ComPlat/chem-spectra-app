from rdkit import Chem
from rdkit.Chem import Descriptors

from chem_spectra.model.helper.share import store_str_in_tmp


class MoleculeModel:
    def __init__(self, molfile):
        self.molfile = molfile.core
        self.mol = self.__set_mol()
        self.smi = self.__set_smi()
        self.mass = self.__set_mass()


    def __set_mol(self):
        tf = store_str_in_tmp(self.molfile, suffix='.mol')
        mol = Chem.MolFromMolFile(tf.name)
        tf.close
        return mol


    def __set_smi(self):
        smi = Chem.MolToSmiles(self.mol, canonical=True)
        return smi


    def __set_mass(self):
        mass = Descriptors.ExactMolWt(self.mol)
        return mass
















#         self.obconv = ob.OBConversion()
#         self.obmol = ob.OBMol()
#         self.__load_to_obmol()
#         self.can = self.__set_can()
#         self.mass = self.__set_mass()


#     def __load_to_obmol(self):
#         self.obconv.SetInAndOutFormats('mol', 'can')
#         self.obconv.ReadString(self.obmol, self.molfile)


#     def __set_can(self):
#         self.can = self.obconv.WriteString(self.obmol)
#         self.can = re.sub('\s+', '', self.can)
#         return self.can


#     def __set_mass(self):
#         self.mass = self.obmol.GetExactMass()
#         return self.mass




# def clear_mapnum(mol):
#     [atom.ClearProp('molAtomMapNumber') for atom in mol.GetAtoms() if atom.HasProp('molAtomMapNumber')]


# def get_unique_fg_smas(smi):
#     results = []
#     mol = Chem.MolFromSmiles(smi)
#     fgs = ifg.identify_functional_groups(mol)

#     for fg in fgs:
#         target = fg.type
#         mol = Chem.MolFromSmarts(target)
#         clear_mapnum(mol)
#         sma = Chem.MolToSmarts(mol)
#         results.append(sma)

#     return list(set(results))
