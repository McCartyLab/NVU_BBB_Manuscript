# NVU_BBB_Manuscript
NVU Remodeling and BBB dysfunction in Primary Gliomas and Brain Metastasis

This repository contains the computational workflows and analysis scripts supporting our manuscript investigating neurovascular unit (NVU) remodeling and blood–brain barrier (BBB) dysfunction across primary gliomas and brain metastases.

## Repository Contents

- scRNA-seq Analysis: Single-cell RNA sequencing processing, quality control, cell-type annotation, and downstream analyses.

- Xenium Large Tissue Analysis: Processing, quality control, cell-type annotation, and spatial transcriptomic analysis of large tissue sections.

- Xenium TMA Analysis: Processing and spatial analysis of Xenium tissue microarray (TMA) datasets.

- NVU State Assignment: Vessel-centered identification and classification of NVU states based on endothelial, pericyte, astrocyte-like, and tumor cell associations, including BBB-related analyses.

- Fiji_ImageJ_Multiplex (Multiplex Image Analysis): Image processing and quantitative analysis of multiplex immunofluorescence data. Custom Fiji/ImageJ macros were used for standardized analysis of COMET multiplex immunofluorescence images, including tumor-region definition, CD31-positive vessel segmentation, vessel and perivascular measurements, and marker-intensity quantification. Each .ijm file was loaded in Fiji using Plugins > Macros > Run, after which the required co-registered channel images, masks or ROI files, and output directory were selected through the on-screen prompts. Analysis-specific macros were applied using consistent parameters across samples within each comparison and generated tabulated measurements, ROI files, and quality-control images, where applicable.

Raw data and Seurat objects will be available at GEO accession GSE349917 once the manuscript is published.




