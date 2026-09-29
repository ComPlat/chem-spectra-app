"""Read a JCAMP file as the blocks it is actually made of.

`ng.jcampdx.read()` returns one flat dict per file: every repeated LDR from
every block accumulates into a single list, and the data arrays from every
block land together under `data['real']`. This app was built on that shape --
`target_idx` indexes both the LDR lists and the data array, and the two line up
only because of how the merge happens to work.

`read_blocks()` returns one dict per block in file order instead, which is what
the file says. This module is the thin layer over it: it answers "which block
holds the measurement" and "what does that block say", so callers stop indexing
a merged list and start naming a block.

Deliberately not an adapter. Nothing here rebuilds the flat dict to keep old
call sites working -- `flat_ldrs()` exists for one caller, the original-metadata
dump, which is specified to write out every LDR in the file.
"""

import nmrglue as ng


class Block:
    """One JCAMP block: its LDRs, and its data if it has any."""

    def __init__(self, raw, index):
        self._raw = raw
        self.index = index
        self._data = None
        self._data_read = False

    @property
    def datatype(self):
        """`##DATA TYPE=`, upper-cased, or '' -- matching how the classifier
        has always compared it."""
        return (self.ldr('DATATYPE') or '').upper()

    @property
    def dataclass(self):
        return (self.ldr('DATACLASS') or '').upper()

    @property
    def parent(self):
        """Index of the enclosing block, or None.

        Not used for choosing the measurement: on `1H.dx` the spectrum block
        reports None while the FID block reports 0, so this does not reliably
        describe LINK nesting. File order does.
        """
        return self._raw.get('_parent')

    def ldr(self, key):
        """First value of `key` in this block, or None."""
        values = self._raw.get(key)
        return values[0] if values else None

    def ldrs(self, key):
        """Every value of `key` in this block, in order."""
        return list(self._raw.get(key) or [])

    def has(self, key):
        return bool(self._raw.get(key))

    @property
    def data(self):
        """This block's data array, read once.

        NTUPLES give `{'real': [...], 'imaginary': [...]}`; `(X++(Y..Y))` gives
        a 1-D array; `(XY..XY)`, XYPOINTS and PEAKTABLE give `(1, N, 2)`; a
        block with no data gives None. nmrglue applies the block's own
        factors, so callers must not apply them again.
        """
        if not self._data_read:
            self._data = ng.jcampdx.getdataarray(self._raw, show_all_data=True)
            self._data_read = True
        return self._data

    def __repr__(self):
        return '<Block {} {!r}>'.format(self.index, self.datatype)


class JcampFile:
    """Every block of one file, in the order the file declares them."""

    def __init__(self, blocks):
        self.blocks = blocks

    def __len__(self):
        return len(self.blocks)

    def __iter__(self):
        return iter(self.blocks)

    def __getitem__(self, index):
        return self.blocks[index]

    @property
    def datatypes(self):
        """Every block's `##DATA TYPE=`, in file order."""
        return [block.datatype for block in self.blocks]

    def first_matching(self, predicate):
        """The first block in file order for which `predicate` holds.

        File order is the rule #291 established for choosing between competing
        datatypes, applied to blocks rather than to a merged list.
        """
        for block in self.blocks:
            if predicate(block):
                return block
        return None

    def blocks_with_datatype(self, predicate):
        return [b for b in self.blocks if predicate(b.datatype)]

    def flat_ldrs(self):
        """Every LDR of every block, concatenated in file order.

        **Only for the original-metadata dump**, which is specified to write
        the whole file back out as `###KEY= v1, v2, ...`. Nineteen golden files
        compare that output byte for byte, including values joined across
        blocks such as `###TITLE= X, X, X`, so the concatenation order is part
        of the contract.

        Any other caller wants a specific block and should name it.
        """
        flat = {}
        for block in self.blocks:
            for key, values in block._raw.items():
                if key.startswith('_'):
                    continue
                flat.setdefault(key, []).extend(values)
        return flat


def block_from_headers(headers=None):
    """A single block standing for a core that is not a JCAMP file.

    `FidBaseConverter`, `NMRiumDataConverter` and `CdfMSConverter` build their
    data themselves and set the headers they need -- `FIRSTX`, `$OFFSET`,
    `.OBSERVEFREQUENCY` -- into a dict of their own. They still reach
    `JcampTechniqueConverter`, which asks the target block for LDRs, so they
    carry one of these: for them that dict *is* the only block.

    Scalars are wrapped so `.ldr()` reads them the same way as a parsed LDR;
    keys nmrglue's Bruker reader nests as sub-dicts are skipped, since they are
    not LDRs.
    """
    raw = {}
    for key, value in (headers or {}).items():
        if isinstance(value, dict):
            continue
        raw[key] = value if isinstance(value, list) else [value]
    block = Block(raw, 0)
    block._data_read = True
    return block


def read_jcamp(path, read_err='ignore'):
    return JcampFile([
        Block(raw, index)
        for index, raw in enumerate(ng.jcampdx.read_blocks(path, read_err=read_err))
    ])
