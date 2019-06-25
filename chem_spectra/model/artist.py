from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D  # Needed to show molecules
from matplotlib import colors as mcolors

from chem_spectra.lib.shared.buffer import store_str_in_tmp

colors = [
    'lightpink',
    'lime',
    'yellow',
    'lightskyblue',
    'lavender',
    'gold',
    'azure',
    'beige',
]


class ArtistModel:
    def __init__(self, molfile=False, predictions=[]):
        is_molfile_str = type(molfile).__name__ == 'str'
        self.molfile = molfile if is_molfile_str else molfile.core
        self.predictions = predictions
        self.mol = self.__set_mol()
        self.fgs = self.__set_fgs()

    def __set_mol(self):
        tf = store_str_in_tmp(self.molfile, suffix='.mol')
        mol = Chem.MolFromMolFile(tf.name)
        tf.close
        return mol

    def __set_fgs(self):
        return [x['sma'] for x in self.predictions]

    @classmethod
    def draw_ir(cls, molfile=False, predictions=[]):
        instance = cls(molfile=molfile, predictions=predictions)
        svgs = instance.__draw_ir()
        return svgs

    def __draw_ir(self):
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)

        for i, fg in enumerate(self.fgs):
            fg = Chem.MolFromSmarts(fg)
            target_atoms = self.mol.GetSubstructMatch(fg)

            target_bonds = []
            for idx, b in enumerate(self.mol.GetBonds()):
                bb = b.GetBeginAtomIdx()
                be = b.GetEndAtomIdx()
                if bb in target_atoms and be in target_atoms:
                    target_bonds.append(idx)

            farber = colors[i % len(colors)]
            color_atoms = {}
            for t in target_atoms:
                color_atoms[t] = mcolors.to_rgba(farber)
            color_bonds = {}
            for t in target_bonds:
                color_bonds[t] = mcolors.to_rgba(farber)

            drawer.DrawMolecule(
                self.mol,
                highlightAtoms=target_atoms,
                highlightAtomColors=color_atoms,
                highlightBonds=target_bonds,
                highlightBondColors=color_bonds,
            )
        svg = drawer.GetDrawingText().replace('svg:', '')
        return [svg]
