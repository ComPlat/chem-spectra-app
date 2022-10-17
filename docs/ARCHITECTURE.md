# Architecture for version 0.10.15

## Table of contents
* [Overview](#overview-a-nameoverview)
* [Components](#components)
* [Sequence diagrams](#sequence-diagrams)
    1. [File api](#1-file-api)
    2. [Transform api](#2-transform-api)
    3. [Inference api](#3-inference-api)

## Overview
![Overview of structure](./diagrams/img/OverView.svg)

## Components
![Backend's component](./diagrams/img/BackendComponents.svg)

## Sequence diagrams
### 1. File api

#### *a. Convert file*
```
Endpoint:
/api/v1/chemspectra/file/convert  [POST]
```
![Convert file](./diagrams/img/FileConvert.svg)

#### *b. Convert molfile*
```
Endpoint: 
/api/v1/chemspectra/molfile/convert [POST]
```
![Convert molfile](./diagrams/img/FileConvertMolfile.svg)

#### *c. Refresh file*
```
Endpoint: 
/api/v1/chemspectra/file/refresh [POST]
```
![Refresh file](./diagrams/img/FileRefresh.svg)

#### *d. Save file*
```
Endpoint: 
/api/v1/chemspectra/file/save [POST]
```
![Refresh file](./diagrams/img/FileSave.svg)

### 2. Transform api
#### *a. Transform to get jcamp files and images as zip format*
```
Endpoint:
/zip_jcamp_n_img  [POST]
```
![Transform get jcamps and images as zip](./diagrams/img/TransformGetZipJcampAndImage.svg)

#### *b. Transform to get jcamp files as zip format*
```
Endpoint:
/zip_jcamp  [POST]
```
![Transform get jcamp files as zip](./diagrams/img/TransformGetZipJcamp.svg)

#### *c. Transform to get image as zip format*
```
Endpoint:
/zip_image  [POST]
```
![Transform get images as zip](./diagrams/img/TransformGetZipImage.svg)

#### *d. Transform to get jcamp file*
```
Endpoint:
/jcamp  [POST]
```
![Transform get jcamp file](./diagrams/img/TransformGetJcamp.svg)

#### *e. Transform to get image in PNG format*
```
Endpoint:
/image  [POST]
```
![Transform get image](./diagrams/img/TransformGetImage.svg)

### 3. Inference api

#### *a. Predict NMR signals with peaks data as FORM request*
```
Endpoint:
/predict/by_peaks_form  [POST]
/api/v1/chemspectra/predict/nmr_peaks_form [POST]
```
![Predict NMR with peaks as FORM request](./diagrams/img/PredictNMRByPeaksFormRequest.svg)

#### *b. Predict NMR signals with peaks data as JSON request*
```
Endpoint:
/predict/by_peaks_json  [POST]
/api/v1/chemspectra/predict/nmr_peaks_json [POST]
```
![Predict NMR with peaks as JSON request](./diagrams/img/PredictNMRByPeaksJSON.svg)

#### *c. Predict Mass spectrum*
```
Endpoint:
/predict/ms  [POST]
/api/v1/chemspectra/predict/ms [POST]
```
![Predict Mass spectrum](./diagrams/img/PredictMassSpectrum.svg)

#### *d. Predict Mass spectrum*
```
Endpoint:
/predict/infrared  [POST]
/api/v1/chemspectra/predict/infrared [POST]
```
![Predict Infrared](./diagrams/img/PredictInfrared.svg)
