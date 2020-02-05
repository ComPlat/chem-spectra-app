import matplotlib
matplotlib.use('Agg')

from rdkit import Chem  # noqa: E402
from rdkit.Chem.Draw import rdMolDraw2D  # noqa: E402  # To show molecules
from matplotlib import colors as mcolors  # noqa: E402


colors = [
    'yellow',
    'lightskyblue',
    'lightpink',
    'lime',
    'lavender',
    'gold',
    'azure',
    'beige',
]

svg_target = "xmlns:xlink='http://www.w3.org/1999/xlink'\n                  xml:space='preserve'\n"  # noqa: E501
svg_vb = "viewbox='0 0 400 400'"
svg_size = "width='400px' height='400px'"


class ArtistModel:
    def __init__(self, mm=False, predictions=[], layout=False):
        self.predictions = predictions
        self.layout = layout
        self.mol = mm.mol

    @classmethod
    def draw_ir(cls, mm=False, predictions=[], layout=False):
        instance = cls(mm=mm, predictions=predictions, layout=layout)
        svgs = instance.__draw_ir()
        return svgs

    def __draw_ir(self):
        fgs = [x['sma'] for x in self.predictions]
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        atom_counts = {}

        for i, fg in enumerate(fgs):
            fg = Chem.MolFromSmarts(fg)
            matches = self.mol.GetSubstructMatches(fg)

            for target_atoms in matches:
                target_bonds = []
                for idx, b in enumerate(self.mol.GetBonds()):
                    bb = b.GetBeginAtomIdx()
                    be = b.GetEndAtomIdx()
                    if bb in target_atoms and be in target_atoms:
                        target_bonds.append(idx)

                farber = colors[i % len(colors)]
                color_atoms = {}
                radius_atoms = {}
                for t in target_atoms:
                    if atom_counts.get(t):
                        atom_counts[t] += 1
                    else:
                        atom_counts[t] = 1
                    color_atoms[t] = mcolors.to_rgba(farber)
                    radius_atoms[t] = 0.3
                    # adjust for overlaps
                    if atom_counts.get(t) == 2:
                        radius_atoms[t] = 0.2
                    elif atom_counts.get(t) >= 3:
                        radius_atoms[t] = 0.1
                color_bonds = {}
                for t in target_bonds:
                    color_bonds[t] = mcolors.to_rgba(farber)

                drawer.DrawMolecule(
                    self.mol,
                    highlightAtoms=target_atoms,
                    highlightAtomColors=color_atoms,
                    highlightAtomRadii=radius_atoms,
                    highlightBonds=target_bonds,
                    highlightBondColors=color_bonds,
                )
        svg = drawer.GetDrawingText().replace('svg:', '').replace(
                svg_target + svg_size,
                "{} {}".format(svg_target, svg_vb),
            )
        return [svg]

    @classmethod
    def draw_nmr(cls, mm=False, predictions=[], layout=False):
        instance = cls(mm=mm, predictions=predictions, layout=layout)
        svgs = instance.__draw_nmr()
        return svgs

    def __identify_targets(self):
        atom_symbol = 'C'
        if self.layout == '1H':
            atom_symbol = 'H'
        elif self.layout == '19F':
            atom_symbol = 'F'

        targets = []
        for idx in range(self.mol.GetNumAtoms()):
            if self.mol.GetAtomWithIdx(idx).GetSymbol() == atom_symbol:
                targets.append(idx)
        return targets

    def __draw_nmr(self):
        targets = self.__identify_targets()
        colors = {}

        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        opts = drawer.drawOptions()

        for t in targets:
            opts.atomLabels[t] = str(t + 1)
            colors[t] = mcolors.to_rgba('yellow')

        drawer.DrawMolecule(
            self.mol,
            highlightAtoms=targets,
            highlightAtomColors=colors,
            highlightBonds=[]
        )
        drawer.FinishDrawing()

        svg = drawer.GetDrawingText().replace('svg:', '').replace(
                svg_target + svg_size,
                "{} {}".format(svg_target, svg_vb),
            )
        return [svg]
