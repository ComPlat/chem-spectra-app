"""The msconvert command carries its filters without literal quotes.

The argv goes to msconvert as a list -- through `docker exec`, or through the
/bin/docker shim to the sidecar -- so no shell strips quotes on the way.
ProteoWizard 3.0 stripped them itself. 4.0 does not: it answers
`Unknown filter '"peakPicking'`, exits 1 and writes no mzML, which reached the
app only as "no mzML appeared".
"""

from types import SimpleNamespace

from chem_spectra.lib.converter.ms import MSConverter


def _command():
    converter = MSConverter.__new__(MSConverter)
    converter.hash_str = 'job'
    converter.tf = SimpleNamespace(name='/tmp/upload.RAW')
    return converter._MSConverter__build_cmd_msconvert()


def test_no_argument_carries_literal_quotes():
    assert not [arg for arg in _command() if '"' in arg or "'" in arg]


def test_the_filters_are_peak_picking_and_zero_samples():
    command = _command()
    filters = [command[i + 1] for i, arg in enumerate(command) if arg == '--filter']
    assert filters == ['peakPicking true 1-', 'zeroSamples removeExtra']


def test_the_prefix_the_sidecar_shim_asserts_is_unchanged():
    # fake-docker.py, with MSC_VALIDATE set, requires exactly this prefix
    assert _command()[:7] == [
        'docker', 'exec', '-d', 'msconvert_docker', 'wine', 'msconvert',
        '/data/job/upload.RAW',
    ]
