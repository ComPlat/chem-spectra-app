"""Record what every fixture produces today, so the nmrglue migration can be
measured rather than argued about.

The migration moves the app from a flat `ng.jcampdx.read()` dict to per-block
reads. Almost every observable is meant to stay identical across that move;
the few that change must be justified one by one. This script writes the
observables down first.

Run it on the current pin to make the baseline, then again after each step:

    python -m tests.tools.snapshot_baseline          # write
    python -m tests.tools.snapshot_baseline --check  # compare, exit 1 on drift

**What it deliberately does not record.** The snapshot has to survive the
refactor it is measuring, so it only reads surfaces that outlive it:

- not `target_idx`, which A3 deletes;
- not descriptor field names, which get renamed;
- the chosen block is identified by its datatype string and its x range, never
  by an index.

Everything is reached through `TransformerModel`, the same entry point the
endpoints use.
"""

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

from werkzeug.datastructures import FileStorage

from chem_spectra.controller.helper.file_container import FileContainer
from chem_spectra.model.transformer import TransformerModel

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / 'tests' / 'fixtures'
# Deliberately not under tests/fixtures/. These are harness output, not
# fixtures, and `test_shorten_label.py` asserts that every basename under
# tests/fixtures is short enough to survive legend truncation -- which the
# flattened snapshot names are not. Keeping them out also stops the harness
# from globbing its own composed .jdx files back in as inputs.
SNAPSHOTS = ROOT / 'tests' / 'snapshots'

# Bruker FIDs arrive as zips and take a different branch through the
# transformer, so they carry the molfile the zip path expects.
MOLFILE = FIXTURES / 'source' / 'molfile' / 'svs813f1_B.mol'


def _digest(value):
    return hashlib.sha1(repr(value).encode()).hexdigest()[:16]


def _array(values):
    """A sequence recorded so that a diff names what moved.

    The hash catches any change at all; first/last/len say where to look.
    """
    if values is None:
        return None
    seq = list(values)
    if not seq:
        return {'len': 0}
    return {
        'len': len(seq),
        'first': round(float(seq[0]), 8),
        'last': round(float(seq[-1]), 8),
        'sha1': hashlib.sha1(
            b''.join(f'{float(v):.10g}'.encode() for v in seq)
        ).hexdigest(),
    }


def _peaks(peaks):
    if not peaks:
        return None
    return {'n': len(peaks.get('x') or []), 'sha1': _digest(peaks)}


def _text(lines):
    """Composed JCAMP. Stored whole, since a hash alone cannot be reviewed."""
    if not lines:
        return None
    body = lines.read() if hasattr(lines, 'read') else ''.join(lines)
    if isinstance(body, bytes):
        body = body.decode('utf-8', 'replace')
    return body


def _observables(converter, composer):
    """The per-fixture record.

    `target` is how the chosen block is identified without an index: its
    declared datatype plus the x range it produced. If A3 selects a different
    block, one of the two moves.
    """
    xs = getattr(converter, 'xs', None)
    ys = getattr(converter, 'ys', None)
    record = {
        'typ': getattr(converter, 'typ', None),
        'threshold': getattr(converter, 'threshold', None),
        'title': getattr(converter, 'title', None),
        'target': {
            'datatype': getattr(converter, 'datatype', None),
            'dataclass': getattr(converter, 'dataclass', None),
        },
        'xs': _array(xs),
        'ys': _array(ys),
        'factor': getattr(converter, 'factor', None),
        'boundary': getattr(converter, 'boundary', None),
        'label': getattr(converter, 'label', None),
        'auto_peaks': _peaks(getattr(converter, 'auto_peaks', None)),
        'edit_peaks': _peaks(getattr(converter, 'edit_peaks', None)),
        # Integrals live in two tables, not one attribute named for them.
        'itg_table': getattr(converter, 'itg_table', None),
        'mpy_itg_table': getattr(converter, 'mpy_itg_table', None),
    }
    record.update(_ms_observables(converter))
    if composer is not None:
        jcamp = _text(composer.tf_jcamp())
        record['jcamp'] = {
            'len': len(jcamp) if jcamp else 0,
            'sha1': hashlib.sha1(jcamp.encode()).hexdigest() if jcamp else None,
        }
        # The text itself goes to a sibling file, not into this JSON. Inlined,
        # a 300 KB output is one JSON string, so any change to it shows up in
        # `git diff` as a single unreadable line -- and A1's whole contract is
        # that a later PR either produces an empty diff or explains each
        # entry. As a file, a changed ##BLOCKS reads as one changed line.
        record['_jcamp_text'] = jcamp
    return json.loads(json.dumps(record, default=_jsonable))


def _ms_observables(converter):
    """MS takes a different converter with a different vocabulary.

    `JcampMSConverter` has no xs/ys, factor or peaks; its spectrum lives in
    `.data` as an (N, 2) array per run. Without this the MS fixtures record
    almost nothing, and a migration could lose a whole run unnoticed -- which
    is exactly what A4 changes, since it replaces the cross-block merge.
    """
    if type(converter).__name__ not in ('JcampMSConverter', 'MSConverter'):
        return {}
    data = getattr(converter, 'data', None)
    runs = []
    if data is not None:
        for run in data:
            pairs = run.reshape(-1, 2) if hasattr(run, 'reshape') else run
            runs.append({
                'mz': _array([p[0] for p in pairs]),
                'intensity': _array([p[1] for p in pairs]),
            })
    return {
        'ms': {
            'runs': runs,
            'n_spectra': len(getattr(converter, 'spectra', []) or []),
            'auto_scan': getattr(converter, 'auto_scan', None),
            'exact_mz': getattr(converter, 'exact_mz', None),
            'thres': getattr(converter, 'thres', None),
            'n_datatables': len(getattr(converter, 'datatables', []) or []),
        }
    }


def _jsonable(value):
    """numpy scalars and arrays reach here; nothing else should."""
    if hasattr(value, 'tolist'):
        return value.tolist()
    return str(value)


def fixtures():
    """Every JCAMP the app can be asked to read, plus the Bruker zips.

    The files under `result/`, `edit/` and `auto/` are composed outputs that
    the app reads back in, so they are inputs too.
    """
    found = []
    for path in sorted(FIXTURES.rglob('*')):
        if path.suffix.lower() in {'.jdx', '.dx', '.jcamp'}:
            found.append(path)
    for path in sorted((FIXTURES / 'source' / 'bruker').glob('*.zip')):
        found.append(path)
    return found


def record(path):
    ext = path.suffix.lower().lstrip('.')
    # to_converter reads params['ext'] unconditionally, so it is always set.
    params = {'ext': ext}
    with open(path, 'rb') as handle:
        container = FileContainer(FileStorage(handle))
        molfile = None
        if ext == 'zip' and MOLFILE.exists():
            with open(MOLFILE, 'rb') as molhandle:
                molfile = FileContainer(FileStorage(molhandle))
                return _drive(container, molfile, params)
        return _drive(container, molfile, params)


def _drive(container, molfile, params):
    model = TransformerModel(container, molfile=molfile, params=params)
    converter = model.to_converter()
    composer, _ = model.to_composer()
    # The zip and LC/MS paths hand back lists; the first entry is the one the
    # endpoints render.
    if isinstance(converter, list):
        converter = converter[0] if converter else None
    if isinstance(composer, list):
        composer = composer[0] if composer else None
    if converter is None or converter is False:
        return {'error': 'no converter'}
    # The zip/FID path hands back a `FidBaseConverter`, which carries none of
    # the observables -- those live on the `JcampTechniqueConverter` the
    # composer was built from. That core is the technique converter on every
    # path, so read from it wherever there is one.
    subject = getattr(composer, 'core', None) if composer else None
    return _observables(subject if subject is not None else converter,
                        composer or None)


def snapshot_path(path):
    return SNAPSHOTS / (
        str(path.relative_to(FIXTURES)).replace(os.sep, '__') + '.json'
    )


def composed_path(path):
    """Where the composed JCAMP for a fixture is written, beside its JSON."""
    return snapshot_path(path).with_suffix('.jdx')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--check', action='store_true',
        help='compare against the stored snapshots instead of writing them',
    )
    args = parser.parse_args()

    SNAPSHOTS.mkdir(parents=True, exist_ok=True)
    drift = []
    for path in fixtures():
        try:
            current = record(path)
        except Exception as err:            # a fixture that stops parsing is
            current = {'error': f'{type(err).__name__}: {err}'}   # itself news
        composed = current.pop('_jcamp_text', None)
        target = snapshot_path(path)
        text_target = composed_path(path)
        body = json.dumps(current, indent=2, sort_keys=True) + '\n'
        if args.check:
            if not target.exists():
                drift.append(f'{path.name}: no snapshot recorded')
            elif target.read_text(encoding='utf-8') != body:
                drift.append(f'{path.name}: observables differ')
            stored = text_target.read_text(encoding='utf-8') if text_target.exists() else None
            if composed != stored:
                drift.append(f'{path.name}: composed JCAMP differs')
        else:
            target.write_text(body, encoding='utf-8')
            if composed is None:
                text_target.unlink(missing_ok=True)
            else:
                text_target.write_text(composed, encoding='utf-8')
            print(f'wrote {target.relative_to(ROOT)}')

    if args.check:
        for line in drift:
            print(line)
        print(f'{len(drift)} drift(s) across {len(fixtures())} fixtures')
        return 1 if drift else 0
    return 0


if __name__ == '__main__':
    sys.exit(main())
