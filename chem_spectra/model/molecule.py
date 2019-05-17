from rdkit import Chem
from rdkit.Chem import Descriptors

from chem_spectra.lib.shared.buffer import store_str_in_tmp
import chem_spectra.lib.chem.ifg as ifg


class MoleculeModel:
    def __init__(self, molfile):
        is_molfile_str = type(molfile).__name__ == 'str'
        self.molfile = molfile if is_molfile_str else molfile.core
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

    def __clear_mapnum(self, mol):
        [
            atom.ClearProp('molAtomMapNumber')
            for atom in mol.GetAtoms()
            if atom.HasProp('molAtomMapNumber')
        ]

    def fgs(self):
        results = []
        fgs = ifg.identify_functional_groups(self.mol)

        for fg in fgs:
            target = fg.type
            mol = Chem.MolFromSmarts(target)
            self.__clear_mapnum(mol)
            sma = Chem.MolToSmarts(mol)
            results.append(sma)

        return list(set(results))
