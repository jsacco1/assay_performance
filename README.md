# Immune assay performance

## Description
This repository in python3 contains scripts to investigate the role of two factors, ambient shipment temperature and turn-around time (TAT), on immune assay performance, as measured in mean spot count for anti-CD3 positive control. This repository serves only to convey conceptual analytical workflows by example on a real-world application. 

## Usage

Pre-requisites: 
   1. Two csv files with relevant columns describing: sample ambient shipment temperature, timestamps for sample collection and processing, cell viability, etc.
   2. A third csv file with immune assay performance of positive controls taken from samples.
   
   

### 1. Data preprocessing and joining with create_df_01.py

   User provides the file names in < >


```
python3 create_df_01.py <samples_1> <samples_2> <assay_results> <output_file>
```

### 2. Analysis

Refer to Jupyter notebook; with it, one can follow step-by-step: 
   1. exploratory data analysis (EDA),
   2. data preprocessing and joining, 
   3. missingness imputation, 
   4. class imbalance mitigation, 
   5. machine learning, 
   6. model evaluation 

Please reach out to me for any further questions. 
