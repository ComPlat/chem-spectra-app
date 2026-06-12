# API Reference

This document lists the HTTP routes exposed by `chem-spectra-app`.

Use `docs/architecture.md` for component boundaries and `docs/core-flows.md` for runtime behavior. This file is the compact route-level reference for handlers, main inputs, and main outputs.

## `file_api`

| Route | Handler | Main input | Main output |
|---|---|---|---|
| `POST /api/v1/chemspectra/file/convert` | `chemspectra_file_convert()` | `file`, optional `molfile`, form params | JSON with base64 `jcamp` and `img`, or BagIt `list_jcamps` |
| `POST /api/v1/chemspectra/file/save` | `chemspectra_file_save()` | `src`, `dst` or `dst_list`, optional `molfile`, form params | `spectrum.zip` containing source, generated files, predictions, optional CSV |
| `POST /api/v1/chemspectra/file/refresh` | `chemspectra_file_refresh()` | `dst` or `dst_list`, optional `molfile`, form params | JSON with base64 `jcamp` and `img`, or `spectrum.zip` for multiple files |
| `POST /api/v1/chemspectra/molfile/convert` | `chemspectra_molfile_convert()` | `molfile` | JSON with `smi`, `mass`, and molecule `svg` |

## `transform_api`

| Route | Handler | Main input | Main output |
|---|---|---|---|
| `POST /zip_jcamp_n_img` | `zip_jcamp_n_img()` | `file`, optional `molfile`, form params | `spectrum.zip` with JCAMP/image/optional CSV, plus `X-Extra-Info-JSON` |
| `POST /zip_jcamp` | `zip_jcamp()` | `file`, optional `molfile`, form params | `spectrum.zip` with JCAMP files |
| `POST /zip_image` | `zip_image()` | `file`, optional `molfile`, form params | `spectrum.zip` with image files |
| `POST /jcamp` | `jcamp()` | `file`, optional `molfile`, form params | `spectrum.jdx` |
| `POST /image` | `image()` | `file`, optional `molfile`, form params | `spectrum.png` |
| `POST /nmrium` | `nmrium()` | NMRium `file` | `spectrum.jdx`, or `404` if conversion fails |
| `POST /combine_images` | `combine_images()` | `files[]`, form params, optional `extras` | `spectrum.zip` with combined image output |

## `inference_api`

| Route | Handler | Main input | Main output |
|---|---|---|---|
| `POST /predict/by_peaks_json` | `chemspectra_predict_by_peaks_json()` | JSON `layout`, `peaks`, `shift`, `molfile` | prediction JSON |
| `POST /api/v1/chemspectra/predict/nmr_peaks_json` | `chemspectra_predict_by_peaks_json()` | JSON `layout`, `peaks`, `shift`, `molfile` | prediction JSON |
| `POST /predict/by_peaks_form` | `chemspectra_predict_by_peaks_form()` | form `layout`, `peaks`, `shift`, `molfile`, optional `spectrum` | prediction JSON |
| `POST /api/v1/chemspectra/predict/nmr_peaks_form` | `chemspectra_predict_by_peaks_form()` | form `layout`, `peaks`, `shift`, `molfile`, optional `spectrum` | prediction JSON |
| `POST /predict/infrared` | `chemspectra_predict_infrared()` | `layout`, `spectrum`, `molfile` | prediction JSON |
| `POST /api/v1/chemspectra/predict/infrared` | `chemspectra_predict_infrared()` | `layout`, `spectrum`, `molfile` | prediction JSON |
| `POST /predict/ms` | `chemspectra_predict_ms()` | `layout`, `spectrum`, `molfile` | prediction JSON |
| `POST /api/v1/chemspectra/predict/ms` | `chemspectra_predict_ms()` | `layout`, `spectrum`, `molfile` | prediction JSON |

## `spectra_layout_api`

| Route | Handler | Main input | Main output |
|---|---|---|---|
| `GET /api/v1/chemspectra/spectra_layouts` | `get_spectra_layouts()` | none | JSON layout mapping from JCAMP data type configuration |
