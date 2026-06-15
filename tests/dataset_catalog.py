"""Resolve centralized test datasets from chem-spectra-test-files."""

from __future__ import annotations

import json
import os
from functools import lru_cache
from pathlib import Path

_APP_ROOT = Path(__file__).resolve().parents[1]
_SUBMODULE_ROOT = _APP_ROOT / "test-datasets"


def test_files_root() -> Path:
    """Return chem-spectra-test-files repository root."""
    env_root = os.environ.get("CHEMSPECTRA_TEST_FILES")
    if env_root:
        return Path(env_root).expanduser().resolve()

    submodule_root = _SUBMODULE_ROOT.resolve()
    catalog_path = submodule_root / "catalog" / "datasets.json"
    if catalog_path.is_file():
        return submodule_root

    raise FileNotFoundError(
        "Centralized test datasets not found. "
        f"Expected submodule at {submodule_root} "
        "(run: git submodule update --init --recursive). "
        "Alternatively, set CHEMSPECTRA_TEST_FILES to the "
        "chem-spectra-test-files repository root."
    )


@lru_cache(maxsize=1)
def _load_catalog() -> dict:
    catalog_path = test_files_root() / "catalog" / "datasets.json"
    if not catalog_path.is_file():
        raise FileNotFoundError(
            f"Missing centralized catalog: {catalog_path}. "
            "Initialize the test-datasets submodule "
            "(git submodule update --init --recursive) or set "
            "CHEMSPECTRA_TEST_FILES to the chem-spectra-test-files root."
        )
    with catalog_path.open(encoding="utf-8") as handle:
        return json.load(handle)


def dataset_path(dataset_id: str) -> Path:
    """Resolve a catalog dataset ID to an absolute file path."""
    catalog = _load_catalog()
    datasets = catalog.get("datasets", {})
    if dataset_id not in datasets:
        raise KeyError(f"Unknown dataset ID: {dataset_id}")

    entry = datasets[dataset_id]
    file_path = (test_files_root() / entry["path"]).resolve()
    if not file_path.is_file():
        raise FileNotFoundError(
            f"Dataset file missing for {dataset_id}: {file_path}"
        )
    return file_path


def dataset_path_str(dataset_id: str) -> str:
    return str(dataset_path(dataset_id))


# Legacy paths under tests/fixtures/source/ mapped to catalog IDs.
LEGACY_SOURCE_PATH_TO_ID: dict[str, str] = {
    "13C-CPD.dx": "NMR-019",
    "13C-DEPT135.dx": "NMR-020",
    "1H.dx": "NMR-021",
    "IR.dx": "IR-004",
    "JPK-948.jdx": "NMR-022",
    "MS.dx": "MS-005",
    "SVS-790A_13C.jdx": "NMR-023",
    "auto/auto_13C-DEPT135.dx": "NMR-024",
    "auto/auto_1H.dx": "NMR-025",
    "auto/auto_IR.dx": "IR-005",
    "edit/edit_13C-DEPT135.dx": "NMR-026",
    "edit/edit_1H.dx": "NMR-027",
    "edit/edit_IR.dx": "IR-006",
    "edit/inherit_1H.dx": "NMR-028",
    "bagit/aif/aif.zip": "AIF-003",
    "bagit/cv/File053_BagIt.zip": "CV-B-001",
    "bagit/dls_acf/dls_acf.zip": "DLS-002",
    "bagit/dls_intensity/dls_intensity.zip": "DLS-003",
    "bagit/dsc/dsc.zip": "DSC-B-002",
    "bagit/emissions/emissions.zip": "EM-002",
    "bruker/1H.zip": "NMR-BZ-019",
    "cyclicvoltammetry/RCV_LSH-R444_full+Fc.jdx": "CV-009",
    "molfile/invalid_molfile.mol": "MOL-001",
    "molfile/svs813f1_B.mol": "MOL-002",
    "ms/MS_ESI.RAW": "MS-R-004",
    "ms/MS_ESI.jdx": "MS-006",
    "ms/ms_v6.dx": "MS-007",
    "ms/svs813f1.jdx": "MS-008",
    "ms/svs813f1.mzML": "MS-M-005",
    "xrd/Test_data_op_001_XRD.jdx": "XRD-004",
    "xrd/Test_data_op_002_XRD.jdx": "XRD-005",
    "CHI-224_10.jdx": "NMR-029",
}


def legacy_source_path(relative_path: str) -> str:
    """Resolve a legacy fixtures/source relative path via the catalog."""
    normalized = relative_path.lstrip("/")
    dataset_id = LEGACY_SOURCE_PATH_TO_ID.get(normalized)
    if dataset_id is None:
        raise KeyError(
            f"No centralized dataset mapping for legacy path: {relative_path}"
        )
    return dataset_path_str(dataset_id)
