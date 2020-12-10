# Immune assay performance

## Description
This repository in python3 contains scripts to investigate the role of two factors, ambient shipment temperature and turn-around time (TAT), on immune assay performance, as measured in mean spot count for anti-CD3 positive control.

## Usage

Pre-requisites: two csv files with relevant columns describing sample ambient shipment temperature, timestamps for sample collection and processing, etc.

1. Data preprocessing and joining with create_df_01.py

   User provides the file names in < >


```
python3 create_df_01.py <samples_1> <samples_2> <assay_results> <output_file>
```

2. Analysis

Refer to Jupyter notebook provided.
