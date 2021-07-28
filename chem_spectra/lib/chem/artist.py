from collections import defaultdict
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


class ArtistLib:
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
        option = drawer.drawOptions()
        option.padding=0.15
        option.bondLineWidth = 2
        option.annotationFontScale = 0.2

        atom_highlights = defaultdict(list)
        bond_highlights = defaultdict(list)
        atom_highlight_rads = {}
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
                for t in target_atoms:
                    if atom_counts.get(t):
                        atom_counts[t] += 1
                    else:
                        atom_counts[t] = 1
                    atom_highlights[t] = [mcolors.to_rgba(farber)]
                    atom_highlight_rads[t] = 0.3
                    # adjust for overlaps
                    # if atom_counts.get(t) == 2:
                    #     atom_highlight_rads[t] = 0.2
                    # elif atom_counts.get(t) >= 3:
                    #     atom_highlight_rads[t] = 0.1
                for t in target_bonds:
                    bond_highlights[t] = [mcolors.to_rgba(farber)]

        drawer.DrawMoleculeWithHighlights(
            self.mol,
            '',
            dict(atom_highlights),
            dict(bond_highlights),
            atom_highlight_rads,
            {}
        )
        drawer.FinishDrawing()

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
        elif self.layout == '31P':
            atom_symbol = 'P'
        elif self.layout == '15N':
            atom_symbol = 'N'
        elif self.layout == '29Si':
            atom_symbol = 'Si'

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
