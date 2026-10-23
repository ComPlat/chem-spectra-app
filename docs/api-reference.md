# API Reference

This document lists the HTTP routes exposed by `chem-spectra-app`.

Use `docs/architecture.md` for component boundaries and `docs/core-flows.md` for runtime behavior. This file is the compact route-level reference for handlers, main inputs, and main outputs.

## Errors

A request the service will not process -- no file, a file it cannot read, data it cannot convert -- is answered with **`422`** and a JSON body naming the reason:

```json
{"error": "this file could not be read as a spectrum. ..."}
```

The ELN shows `error` to the user. The same shape comes back with **`502`** when the msconvert sidecar is unavailable.

The **prediction routes** (`/predict/*`) are different: they always answer **`200`**, with the reason in `outline`, because the ELN keeps a prediction response only when the status is 200 and the editor's forecast panel reads `outline.code`:

```json
{"outline": {"code": 503, "text": "IR prediction is not configured on this server."}}
```

| `outline.code` | meaning | the editor shows |
|---|---|---|
| `200` | the prediction, under `output` | the prediction |
| `400` | the predictor could not assign these peaks to this molecule, or the IR spectrum could not be read | "Peak & Element count mismatch" |
| `422` | nothing to ask about: no or unreadable molfile; no peaks, or peaks that are not a list of objects with a numeric `x`; a `shift` that is not an object; no spectrum, or one that cannot be read (for MS prediction, one that is not a mass spectrum); an NMR layout other than 1H, 13C, 19F | "Unknown state" |
| `502` | the predictor answered with an error, or with something this service cannot read; `text` quotes what came back | "Service is not available" |
| `503` | the predictor is not configured here (`URL_NSHIFTDB`, `URL_DEEPIR` unset), or cannot be reached; also msconvert, when a RAW or mzML upload has to be converted first | "Service is not available" |
| `504` | the predictor did not answer within `PREDICT_TIMEOUT` seconds (default 120) | "Service is not available" |

`text` gives the reason in each case. A predictor that answers with an `outline` of its own, as DeepIR does, keeps its code and reason. The HTTP status is set in one place, `_respond()` in `inference_api.py`, so that it can follow `outline.code` once the ELN reads the body of a non-200 answer.

## `file_api`

| Route | Handler | Main input | Main output |
|---|---|---|---|
| `POST /api/v1/chemspectra/file/convert` | `chemspectra_file_convert()` | `file`, optional `molfile`, form params | JSON with base64 `jcamp` and `img`; or `list_jcamps` for a BagIt archive or a processed Bruker upload (one JCAMP per spectrum, as `/zip_jcamp_n_img` returns them) |
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
| `POST /nmrium` | `nmrium()` | NMRium `file` | `spectrum.jdx` |
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
| `GET /api/v1/chemspectra/spectra_layouts` | `fetch_mapping()` | none | JSON layout mapping from JCAMP data type configuration |
