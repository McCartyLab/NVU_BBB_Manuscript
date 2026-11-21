## Spatially Resolved Blood-Brain Barrier Stability States Across Gliomas and Brain Metastases

NVU Remodeling and BBB dysfunction in Primary Gliomas and Brain Metastasis

This repository contains the computational workflows and analysis scripts supporting our manuscript investigating neurovascular unit (NVU) remodeling and blood–brain barrier (BBB) dysfunction across primary gliomas and brain metastases.

## Repository Contents

Each folder contains scripts used to process specific datasets mentioned under the folder details below:

- scRNAseq_analysis: Single-cell RNA sequencing processing, quality control, cell-type annotation, and downstream analyses of Breast and Melanoma datasets.

- Xenium_TMA_WholeTissue_analysis: Xenium Large Tissue Analysis: Processing, quality control, cell-type annotation, and spatial transcriptomic analysis of large tissue sections. Xenium TMA Analysis: Processing and spatial analysis of Xenium tissue microarray (TMA) datasets. The seurat objects (Whole_tissue_annotated.rds, Tumor_core_annotated.rds) required for running the scripts can be found at GEO accession mentioned below. 

- NVU State Assignment: Vessel-centered identification and classification of NVU states in whole tissue (Tumor_core_annotated.rds) based on endothelial, pericyte, astrocyte-like, and tumor cell associations, including BBB-related analyses.

- Fiji_ImageJ_Multiplex (Multiplex Image Analysis): Image processing and quantitative analysis of multiplex immunofluorescence data. Custom Fiji/ImageJ macros were used for standardized analysis of COMET multiplex immunofluorescence images, including tumor-region definition, CD31-positive vessel segmentation, vessel and perivascular measurements, and marker-intensity quantification. 
	
	- Each ".ijm" file was loaded in Fiji using Plugins > Macros > Run, after which the required co-registered channel images, masks or ROI files, and output directory were selected through the on-screen prompts. 
        
	- Analysis-specific macros were applied using consistent parameters across samples within each comparison and generated tabulated measurements, ROI files, and quality-control images, where applicable.

- Renalyzed_Public_dataset: contains script to reanalyse and perform cross dataset comparison of Bejarano et. al, Xie et. al and dataset generated in this manuscript. 

- Avastin_analysis: contains scripts to analyze perivascular endogenous IgG signal in vehicle- and Avastin-treated tumors.


Raw data and Seurat objects will be available at GEO accession GSE349917 once the GEO processed our request to make it public.



