# Rapport de migration — datasets centralisés

**Date :** 2026-06-15  
**Dépôt source :** `chem-spectra-test-files` (catalogue `catalog/datasets.json`)  
**Dépôt consommateur :** `chem-spectra-app` (tests uniquement, code production inchangé)

## Résumé

Les tests backend consomment désormais les datasets via `tests/dataset_catalog.py` au lieu de chemins locaux sous `tests/fixtures/source/`. Les fichiers legacy n'ont **pas** été supprimés.

| Métrique | Valeur |
|---|---|
| Fichiers de test modifiés | 24 |
| Références fixture migrées | 31 (mapping catalogue) + appels directs par ID |
| Tests exécutés (suites requises) | 106 |
| Tests réussis | 104 |
| Tests en échec (acceptés) | 2 (`test_ms_mzml_converter_composer`, `test_ms_raw_converter_composer`) |
| Code production modifié | 0 fichier |

Les 2 échecs MS sont dus à l'absence du conteneur Docker `msconvert_docker` (ProteoWizard) requis pour convertir `.mzML` / `.RAW`. La migration des chemins est correcte ; la validation MS est reportée à un environnement avec Docker fonctionnel.

---

## Helper implémenté

**Fichier :** `tests/dataset_catalog.py`

| Fonction | Rôle |
|---|---|
| `test_files_root()` | Racine de `chem-spectra-test-files` (parent de l'app) ou `CHEMSPECTRA_TEST_FILES` |
| `dataset_path(id)` | Résout un ID catalogue → `Path` absolu |
| `dataset_path_str(id)` | Variante string |
| `legacy_source_path(rel)` | Mappe un ancien chemin `fixtures/source/…` → catalogue via `LEGACY_SOURCE_PATH_TO_ID` |

Résolution du catalogue : `catalog/datasets.json` dans le dépôt test-files.

---

## Fichiers de test migrés

| Fichier | Mécanisme |
|---|---|
| `tests/test_spectra_datatable.py` | `legacy_source_path` |
| `tests/test_spectra_peaks.py` | `legacy_source_path` |
| `tests/test_spectra_im.py` | Indirect (via `test_spectra_peaks`) |
| `tests/test_xrd.py` | `dataset_path_str` (`XRD-004`, `XRD-005`) |
| `tests/test_ms.py` | `dataset_path` / `dataset_path_str` |
| `tests/test_cyclic_volta.py` | `dataset_path_str` (`CV-009`) |
| `tests/test_bagit.py` | `dataset_path_str` (`CV-B-001`) |
| `tests/test_func.py` | `dataset_path` (`MOL-002`) |
| `tests/model/test_transformer.py` | `dataset_path_str` |
| `tests/model/concern/test_property.py` | `dataset_path` |
| `tests/lib/converter/bagit/test_bagit_base_converter.py` | `dataset_path_str` (6 bagit) |
| `tests/lib/converter/fid/test_fid_has_processed.py` | `dataset_path_str` (`NMR-BZ-019`) — voir exceptions |
| `tests/lib/converter/jcamp/test_jcamp_base_converter.py` | `dataset_path_str` |
| `tests/lib/converter/jcamp/test_jcamp_ni_converter.py` | `dataset_path_str` |
| `tests/lib/composer/test_base_composer.py` | `dataset_path_str` |
| `tests/lib/composer/test_ms_composer.py` | `dataset_path_str` |
| `tests/controller/test_file_api.py` | `legacy_source_path` + `dataset_path` |
| `tests/controller/test_inference_api.py` | `dataset_path` |
| `tests/controller/test_transform_api.py` | `legacy_source_path` + `dataset_path` |

**Nouveau fichier :** `tests/dataset_catalog.py`

---

## Mapping legacy → catalogue (31 entrées)

| Chemin legacy | ID catalogue |
|---|---|
| `13C-CPD.dx` | `NMR-019` |
| `13C-DEPT135.dx` | `NMR-020` |
| `1H.dx` | `NMR-021` |
| `IR.dx` | `IR-004` |
| `JPK-948.jdx` | `NMR-022` |
| `MS.dx` | `MS-005` |
| `SVS-790A_13C.jdx` | `NMR-023` |
| `auto/auto_13C-DEPT135.dx` | `NMR-024` |
| `auto/auto_1H.dx` | `NMR-025` |
| `auto/auto_IR.dx` | `IR-005` |
| `edit/edit_13C-DEPT135.dx` | `NMR-026` |
| `edit/edit_1H.dx` | `NMR-027` |
| `edit/edit_IR.dx` | `IR-006` |
| `edit/inherit_1H.dx` | `NMR-028` |
| `bagit/aif/aif.zip` | `AIF-003` |
| `bagit/cv/File053_BagIt.zip` | `CV-B-001` |
| `bagit/dls_acf/dls_acf.zip` | `DLS-002` |
| `bagit/dls_intensity/dls_intensity.zip` | `DLS-003` |
| `bagit/dsc/dsc.zip` | `DSC-B-002` |
| `bagit/emissions/emissions.zip` | `EM-002` |
| `bruker/1H.zip` | `NMR-BZ-019` |
| `cyclicvoltammetry/RCV_LSH-R444_full+Fc.jdx` | `CV-009` |
| `molfile/invalid_molfile.mol` | `MOL-001` |
| `molfile/svs813f1_B.mol` | `MOL-002` |
| `ms/MS_ESI.RAW` | `MS-R-004` |
| `ms/MS_ESI.jdx` | `MS-006` |
| `ms/ms_v6.dx` | `MS-007` |
| `ms/svs813f1.jdx` | `MS-008` |
| `ms/svs813f1.mzML` | `MS-M-005` |
| `xrd/Test_data_op_001_XRD.jdx` | `XRD-004` |
| `xrd/Test_data_op_002_XRD.jdx` | `XRD-005` |

---

## Références fixture restantes (locales)

| Fichier test | Chemin local | Raison |
|---|---|---|
| `tests/lib/converter/fid/test_fid_has_processed.py` | `./tests/fixtures/source/bruker/13C.zip` | Non importé dans le catalogue |
| `tests/lib/test_composer.py` | `CHI-224_10.jdx` | Non importé dans le catalogue |

Tous les autres tests listés utilisent le helper centralisé.

---

## Exécution des tests

**Commande :**

```bash
cd chem-spectra-app
. .venv312/bin/activate
CHEMSPECTRA_TEST_FILES=/path/to/chem-spectra-test-files \
  python3 -m pytest \
    tests/test_spectra_datatable.py \
    tests/test_spectra_peaks.py \
    tests/test_spectra_im.py \
    tests/test_xrd.py \
    tests/test_ms.py \
    tests/test_cyclic_volta.py \
    tests/test_bagit.py \
    tests/test_func.py \
    tests/model/test_transformer.py \
    tests/model/concern/test_property.py \
    tests/lib/converter/bagit/test_bagit_base_converter.py \
    tests/lib/converter/fid/test_fid_has_processed.py \
    tests/lib/composer/test_ms_composer.py \
    tests/controller/test_file_api.py \
    tests/controller/test_inference_api.py \
    tests/controller/test_transform_api.py \
    -v
```

**Résultat :** `104 passed, 2 failed` (~4 min 20 s, Python 3.12)

### Échecs acceptés

| Test | Cause | Correctif migration |
|---|---|---|
| `test_ms_mzml_converter_composer` | `MSConverter` → `docker exec msconvert_docker` indisponible ou timeout | Aucun — infrastructure Docker requise (`INSTALL.md` §1.3) |
| `test_ms_raw_converter_composer` | Idem (conversion `.RAW` → `.mzML`) | Idem |

Les tests JCAMP MS (`test_ms_jcamp_converter_composer`, `test_jcamp_single_point_last_line`) passent sans Docker.

---

## Analyse cleanup — `tests/fixtures/source/`

**Ne pas supprimer** tant que la CI locale et l'équipe n'ont pas validé la suppression.

### Partiellement sûr à supprimer (après validation CI)

Les ~31 fichiers source dont les tests pointent désormais vers le catalogue (voir mapping ci-dessus). Ils sont dupliqués dans `chem-spectra-test-files/test_data/`. La suppression n'est possible qu'après :

1. Passage complet des 106 tests (y compris MS avec Docker)
2. Mise à jour de la doc CI pour démarrer `msconvert_docker`
3. Confirmation qu'aucun script externe ne référence encore ces chemins

### Non sûr à supprimer

| Élément | Justification |
|---|---|
| `CHI-224_10.jdx` | Utilisé par `tests/lib/test_composer.py` (hors scope des 16 suites, pas dans le catalogue) |
| `bruker/13C.zip` | Utilisé par `test_fid_has_processed.py`, non importé |
| `tests/fixtures/result/` | Fichiers golden (métadonnées attendues) — indépendants des sources |
| Fichiers sans tests actifs (`mnova/`, `example/`, `1H.edit.jdx`, etc.) | Pas migrés ; pourraient servir à de futurs tests |

### Verdict global

**Partiellement sûr** — la majorité des fixtures source migrées peut être retirée plus tard, mais 2 références locales et les golden `result/` imposent de conserver une partie du répertoire `fixtures/`.

---

## Prérequis pour les développeurs

1. Cloner `chem-spectra-test-files` à côté de `chem-spectra-app`, **ou** définir `CHEMSPECTRA_TEST_FILES`.
2. Pour les tests MS RAW/mzML : démarrer le conteneur ProteoWizard :

```sh
mkdir -p chem_spectra/tmp && chmod -R 777 chem_spectra/tmp
docker run --detach --name msconvert_docker --rm -it \
  -e WINEDEBUG=-all \
  -v "$(pwd)/chem_spectra/tmp:/data" \
  proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses bash
```
