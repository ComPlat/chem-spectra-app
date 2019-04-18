from rdkit import Chem
import hashlib

import deep_ir.lib.chem.ifg as ifg


def hashrize(name):
    b_name = str.encode(name)
    hash_name = hashlib.md5(b_name).hexdigest()
    return hash_name


def clear_mapnum(mol):
    [atom.ClearProp('molAtomMapNumber') for atom in mol.GetAtoms() if atom.HasProp('molAtomMapNumber')]


def get_unique_fg_smas(smi):
    results = []
    mol = Chem.MolFromSmiles(smi)
    fgs = ifg.identify_functional_groups(mol)

    for fg in fgs:
        target = fg.type
        mol = Chem.MolFromSmarts(target)
        clear_mapnum(mol)
        sma = Chem.MolToSmarts(mol)
        results.append(sma)

    return list(set(results))
