# treatment-graphical-modeling
Analysis code for "Separating the genetics of disease, treatment, and treatment response using graphical modeling and large-scale electronic health records"

## Dependencies

- https://github.com/medical-genomics-group/ci-gwas
- https://github.com/nickmachnik/bicorr
- https://github.com/ippas/wes-qc-loftee-annot
  
## Structure
The analysis is divided into the following sections:

### 1. bp-data-extraction
This creates a database in the UBK DNANexus RAP platform with all of the prescriptions, diagnoses and GP records, extracts and preliminarly cleans BP (blood pressure) and LDL measurments together with important medications. There is also a notebook that prepares dicionaries of diagnostic codes and medication based on tables provided by the UKB.

### 2. bp-data-aggregation 
This part includes cleaning and aggregation of the extracted data into the pre- post- design into all of the sensitivity analyses.

### 3. wes-analysis
This part uses the loftee_annot sofware to evaluate loss of function variants from UKB WES data.

### 4. simulation

### 5. ci-gwas 
This part contains scripts for input prep, submission scripts for graphical modeling as well as initial result aggregation.

### 6. standard-gwas
Scripts for standard GWAS (Regenie) for comparison with the graphical modeling.

### 6. arb-response
Pharmacogenomic followup for ARB-response.

### 7. ukb-analysis-figure
Code for main and supplementary figures as well as the downstream data analisis.

