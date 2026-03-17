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

### 7. arb-response
Pharmacogenomic followup for ARB-response.

### 8. ukb-analysis-figure
Code for main and supplementary figures as well as the downstream data analysis.

# How to apply our framework
Folders above provide code used for the publication which has multiple modelling setups.

Below is a walkthrough of the steps one has to go through to analyse their own pre-post data with our framework. The method itself is documented in a [dedicated repository](https://github.com/medical-genomics-group/ci-gwas). Here we are using the cuskss (summary statistic) version.

#### STEP 1 trait data aggregation

For a pre-post design raw data after QC should be aggregated per individual per row. If traits from multiple age groups for an individual are included these should also be in the same row, see [this file](https://github.com/medical-genomics-group/treatment-graphical-modeling/blob/main/bp-data-aggregation/bp_aggregation.py) for helper functions on the aggregation. Binary (e.g. treatment indicator) phenotypes can be mixed in the same table with continous phenotypes (e.g. measurements). Continous traits should be standardised and adjusted for selected covariates that one does not want to explicitly model ([see functions here](https://github.com/medical-genomics-group/treatment-graphical-modeling/blob/main/bp-data-aggregation/adjust_and_pxp.py)).

For each of the traits you want to model decide on the time index. If all traits should be conditioned on each other in the analysis then give them all a time index of 1. If there is some time-dependency involved (e.g. pre- treatment mesurements, treatment indicators, post-treatment measurements) these sets of traits can be given subsequent time indices.

Note: Phenotypes can have missing values (coded as `nan`), but for the method to work correctly each pair of phenotypes should have some overlapping samples. For binary phenotypes this overlap should incldue both cases as controls, as otherwise it is not possible to compute trait x trait correlations.

#### STEP 2 genetic data preparation

Genetic data should be in plink bed/bim/fam format and have accompanying LD matrices computed per block. First, determine block edges using [`ci-gwas.py block`](https://github.com/medical-genomics-group/ci-gwas/blob/main/ci-gwas.py) (run this per chromosome). Then compute your LD matrices in blocks [with this script](https://github.com/medical-genomics-group/treatment-graphical-modeling/blob/main/ci-gwas/mxm/calc_mxm_per_block.py).

Note: If the LD matrices come from a different source than the genetic data the variant id and positions (order) have to be identical in both the original genetic data and in the LD matrices.

#### STEP 3 CI-GWAS input preparation

To run cuskss one has to input the following files:

1. trait x trait correlation and standard error matrices (`pxp` and `pxp-se` in the cuskss input)
2. marker x trait correlation and standard error matrices (`mxp` and `mxp-se` in the cuskss input)

Both can be computed using the [bicorr python package](https://github.com/nickmachnik/bicorr), [see examples for pxp here](https://github.com/medical-genomics-group/treatment-graphical-modeling/blob/main/bp-data-aggregation/adjust_and_pxp.py) [and for mxp here](https://github.com/medical-genomics-group/treatment-graphical-modeling/blob/main/ci-gwas/mxp/calc_mxp.py).

3. per-block LD matrices (`mxm` in the cuskss input)
4. a file with the list of blocks and their edges (`blockfile` in the input)

LD matrices calculation is described in STEP 2 above. The blockfile is an output of `ci-gwas.py block`.

5. a time index file (`time-index` in the input)

A text file, where for each of the variables in the pxp there should be one row in this file with the index.

Note: The order of variants in the LD matrices and in the mxp files has to be identical. The order of traits in the pxp rows and columns and in the mxp columns needs to be identical.

#### STEP 4 CI-GWAS

Running CI-GWAS cuskss consists of a series of steps:
1. first cuskss run: [see here for example submission script and parameters](https://github.com/medical-genomics-group/treatment-graphical-modeling/blob/main/ci-gwas/submit-main-analysis/1_submit.sh)

Then one has to prepare inputs for the second run.

2. run [merge-block-outputs](https://github.com/medical-genomics-group/treatment-graphical-modeling/blob/main/ci-gwas/submit-main-analysis/2_submit_merge_blocks.sh)
3. [recompute mxm](https://github.com/medical-genomics-group/treatment-graphical-modeling/blob/main/ci-gwas/submit-main-analysis/calc_mxm_merged.py)

4. second cuskss run [see here for an example submission script](https://github.com/medical-genomics-group/treatment-graphical-modeling/blob/main/ci-gwas/submit-main-analysis/4_submit_after_merge.sh)
5. Create final output files and compute FDR [see this script](https://github.com/medical-genomics-group/treatment-graphical-modeling/blob/main/ci-gwas/submit-main-analysis/create_table_cuskss.py). To compute FDR you will need to compute effective samples size from `mxp` correlation and se files [see this script for implementation](https://github.com/medical-genomics-group/treatment-graphical-modeling/blob/main/ci-gwas/mxp/compute_effective_n.py)

The final table will have all of the variants that remained in the graph together with their traits and parameters.



