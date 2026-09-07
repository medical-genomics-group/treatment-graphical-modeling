# Technical reference

This file lists all of the scripts in this repository with a short description of what it does, what the user needs to set (e.g. replacement of placeholders with paths to files) and what are the inputs and outputs are.

**Recurring inputs.** Several tables below refer to files produced earlier in the pipeline:

| Name used below | Produced by |
| --- | --- |
| *record database* | `bp-data-extraction/2-drug-database-construction.ipynb` (on RAP) |
| *raw BP / LDL table* | `bp-data-extraction/3-bp-db.ipynb`, `5-ldl-fetch.ipynb` |
| *pre/post table* | `bp-data-aggregation/record_aggregation.ipynb` |
| *pxp* | trait x trait correlation and SE matrices, `cor_pxp.tsv` / `ses_pxp.tsv` |
| *mxp* | marker x trait correlation and SE matrices, `c{chr}_cors.tsv` / `c{chr}_sds.tsv` |
| *mxm* | per-block marker x marker LD matrices, `c{chr}_b{block}.ld.bin` |
| *CI-GWAS results table* | `ci-gwas/submit-main-analysis/create_table_cuskss.py`, i.e. `cuskss_selected_markers.tsv` |

**Genotype data.** PLINK bed/bim/fam is assumed throughout, both per-chromosome
(`c{1..22}.bed`) and merged (`autosomes.bed`). **Variant identifiers and order must be identical
between the bedfiles, the LD matrices and the mxp files.**

---

## Order of execution

```
1. bp-data-extraction      (RAP)        raw records  ->  BP / LDL / CVD tables
2. bp-data-aggregation     (local)      pre/post design, covariate adjustment, pxp
3. wes-analysis            (RAP)        LoF annotation
4. ci-gwas                 (HPC + GPU)  mxm, mxp, time index, cuskss runs, results table
5. standard-gwas           (HPC)        Regenie comparison
6. arb-response            (RAP + HPC)  pharmacogenomic follow-up
7. ukb-analysis-figures    (local)      figures and tables
```

---

## 1. bp-data-extraction

| File | What it does | Inputs | Outputs | You must set |
| --- | --- | --- | --- | --- |
| `lookups-download.sh` | Downloads the UKB primary-care code lookup tables | public UKB URL | `primarycare_codings.zip` and its contents | nothing *(run locally)* |
| `pheno_prep.py` | Helper module: builds recoding dictionaries for BP measurements, heart-disease diagnoses and drug codes from the UKB lookup tables | lookup tables, `bp_drugs.tsv` | *(imported)* | lookup files must sit in `data/lkps/` |
| `1-pheno_prep.ipynb` | Runs `pheno_prep`, writes the code dictionaries | lookup tables | `recoding_dict.json` | nothing *(local, relative paths)* |
| `recoding_dict.json` | Committed output of the above; dictionaries used downstream | - | - | - |
| `recode_anno_filter.py` | Hail/Spark module: BP parsing and QC, prescription annotation, CVD flags, LDL extraction | RAP database tables | see notebooks 2 to 6 | RAP database URI is passed in by the notebooks |
| `2-drug-database-construction.ipynb` | Builds the main record database (hospital, death, GP) as a Spark database on RAP | UKB dispensed datasets | *record database* (`dnax://`) | RAP project and database name |
| `3-bp-db.ipynb` | Extracts and cleans BP measurements, joins BP-lowering medication | *record database*, `recoding_dict.json` | `bp-no-zero.tsv`, `bp-with-zero.tsv` | RAP database URI |
| `4-heart-disease.ipynb` | Extracts cardiovascular diagnoses | *record database* | `cvd.tsv` | RAP database URI |
| `5-ldl-fetch.ipynb` | Extracts and cleans LDL measurements, joins statins | *record database* | `ldl.tsv` | RAP database URI |
| `6-get-unrelated-ind.ipynb` | Derives the unrelated-individual sets | UKB relatedness file | `up-to-second-degree.tsv`, `up-to-third-degree.tsv` | `<PATH_TO_UKB_REL_DAT>` |

## 2. bp-data-aggregation

| File | What it does | Inputs | Outputs | You must set |
| --- | --- | --- | --- | --- |
| `bp_aggregation.py` | windowing of measurements around the first prescription, mapping substances to drug classes, year-of-birth and age-at-first/last-measurement metadata, QC, and the pre/post aggregation | raw BP / LDL / CVD tables | see notebooks below | nothing *(paths passed as arguments)* |
| `record_aggregation.ipynb` | Runs the aggregation for every modelling setup in the paper | `bp-no-zero.tsv`, `ldl.tsv`, `cvd.tsv` | `bp_pre_post_{lo}to{hi}*.tsv`, `ldl_pre_post_*.tsv` one *pre/post table* per setup | `BP_RAW_TSV`, `LDL_RAW_TSV` (empty strings in the first cell) |
| `adjust_and_pxp.py` | covariate adjustment by OLS residualisation, z-scoring, and trait x trait correlations (via `bicorr`) | *pre/post table*, covariate file, `.fam` | see notebooks below | nothing *(paths passed as arguments)* |
| `adjust-and-make-pxp.ipynb` | Runs the adjustment and builds the *pxp* matrices for every setup | *pre/post tables*, covariates, `.fam` | per setup: `cor_pxp.tsv`, `ses_pxp.tsv`, `cor_pxp_heatmap.png`, `manifest.tsv`; plus the adjusted per-trait files `{setup}__{trait}.tsv` consumed by the mxp step | `BASE_DIR`, `COVAR_FILE`, `OUT_PXP_ROOT`, `FAM_PATH` |
| `analysis_and_figs/plots.py` | the BP-characterisation panels (time since therapy start, number of measurements, GP vs assessment-centre source, measurement patterns, SBP/DBP densities) | unaggregated BP table | see notebook below | - |
| `analysis_and_figs/bp_characterisation_figure.ipynb` | Renders those panels | unaggregated BP table | `plot_a` ... `plot_f` SVGs | `bp_tsv` |
| `analysis_and_figs/drug_action_models.ipynb` | Mixed-model estimates of per-class drug effects on BP | *pre/post table* | `bp_effects.csv` | paths in the first cell |


## 3. wes-analysis

Loss-of-function annotation. Runs on **RAP**. Installable as the `loftee_annot` package - see
[`wes-lof-annot/README.md`](wes-analysis/wes-lof-annot/README.md).

| File | What it does | Inputs | Outputs | You must set |
| --- | --- | --- | --- | --- |
| `wes-lof-annot/preprocessing/install_vep.py`, `install_vep.sh` | Installs Ensembl VEP and LOFTEE on the RAP cluster. Console script: `install_vep` | - | VEP installation | - |
| `wes-lof-annot/analysis/utils/load_spark.py` | Spark / Hail session setup for the RAP project | - | *(imported)* | `<PROJECT_DIR>`, `<RESOURCES_DIR>` |
| `wes-lof-annot/analysis/utils/dxpathlib.py` | `PathDx`, a `pathlib` wrapper for `dnax://` and `/mnt/project` paths | - | *(imported)* | - |
| `wes-lof-annot/analysis/utils/variant_filtering.py` | `VCFFilter`, the variant QC filters | - | *(imported)* | - |
| `wes-lof-annot/analysis/utils/vep-config.json` | VEP / LOFTEE configuration | - | - | - |
| `wes-lof-annot/analysis/cmd/split_vep.py` | Two console scripts. `annotate_vcf <chrom>` runs VEP and LOFTEE over the population-level pVCFs and writes a Hail MatrixTable; `rare_variants_table <chrom> <eids>` aggregates high-confidence LoF variants into a table | UKB population-level exome pVCFs, EID list | annotated MatrixTable, per-chromosome LoF table | chromosome and EID path as CLI arguments |
| `joint-data-prep/gene-pos.py` | Gene coordinates in three passes: GRCh37 GTF, GRCh38 GTF plus UCSC liftover, manual annotation | GRCh37 and GRCh38 GTFs, `hgnc.txt`, `after-liftover.bed`, `manually-filled-genes.csv` | `genes_pos.csv`, `to_fill_manually.csv` | GTF and HGNC paths (edited in place); the liftover done externally at [UCSC](https://genome.ucsc.edu/cgi-bin/hgLiftOver) |


## 4. ci-gwas

**Run `ci-gwas.py block` from the
[CI-GWAS repository](https://github.com/medical-genomics-group/ci-gwas) first** - to obtain the `blockfile`

### 4a. Marker x marker LD (`mxm/`)

| File | What it does | Inputs | Outputs | You must set |
| --- | --- | --- | --- | --- |
| `mxm/calc_mxm_per_block.py` | For each LD block on one chromosome, calls PLINK `--r triangle bin4` to write that block's LD matrix | blockfile, PLINK bedfiles | `c{chr}_b{block}.ld.bin` | CLI arguments: `chrom blockfile bfile_stem outdir plink_binary` |

### 4b. Marker x trait correlations (`mxp/`)

| File | What it does | Inputs | Outputs | You must set |
| --- | --- | --- | --- | --- |
| `mxp/make_args.py` | writes one job-array line per setup x chromosome | adjusted per-trait files | `args.txt` | `path/to/phenotype/files`, `path/to/mxp/output` |
| `mxp/calc_mxp.py` | Computes marker x trait correlations and SEs for one setup x chromosome, in 10,000-variant chunks, using `bicorr` (biserial) or Pearson for `*ADJ` traits | bedfiles, `cor_pxp.tsv` (for trait order), adjusted per-trait files | `c{chr}_cors.tsv`, `c{chr}_sds.tsv` | `path/to/bedfiles`, `path/to/pxp` |
| `mxp/submit_mxp.sh` | SLURM array wrapper around `calc_mxp.py` | `args.txt` | as above | array size, and the args filename (currently `args-for-rerun.txt`) |
| `mxp/split_datasets_mpx.py` | Derives the sensitivity-analysis mxp folders from the pooled ones | pooled mxp folders | one mxp folder per setup | `path/to/ci-gwas` |
| `mxp/split_mxp.ipynb` | Runs the three functions in `split_datasets_mpx.py` | as above | as above | - |
| `mxp/join_chr_mxp.py` | Concatenates the 22 per-chromosome mxp files into one, checking the row count against the `.bim` | `c{1..22}_cors.tsv` / `_sds.tsv` | `merged_mxp_cors.tsv`, `merged_mxp_sds.tsv` | `--mxp-root`, `--bim` |
| `mxp/submit_merge_mxp.sh` | SLURM wrapper for the merge | as above | as above | `path/to/mxp/pruned`, `path/to/bedfiles/autosomes.bim` |
| `mxp/compute_effective_n.py` | Effective sample size per trait from the median correlation and SE; needed for the FDR step | `merged_mxp_cors.tsv`, `merged_mxp_sds.tsv` | `effective_n.json` per setup | `--root` (defaults to `path/to/mxp/main`) |

### 4c. Trait x trait correlations (`pxp/`) and the time index

| File | What it does | Inputs | Outputs | You must set |
| --- | --- | --- | --- | --- |
| `pxp/split_pxp.py` | derives sensitivity-analysis pxp by row and column selection. *(The pxp matrices are built in `bp-data-aggregation/adjust-and-make-pxp.ipynb`.)* | pooled pxp folders | `cor_pxp.tsv`, `ses_pxp.tsv` per setup | `path/to/ci-gwas` |
| `time-index/time_index.py` | assigns each trait a time index from its name | `cor_pxp.tsv` per setup | `{setup}_time_index.txt` | one or two CLI arguments: `PXP_ROOT [OUTDIR]` |

### 4d. Running the model (`submit-main-analysis/`)

| File | What it does | Inputs | Outputs | You must set |
| --- | --- | --- | --- | --- |
| `make_args.py` | Writes one job-array line per chromosome x block x setup | blockfiles | `args-2.txt` | `path/to/blockfiles/`, the `GEN_SETUPS` list |
| `preflight_check.py` | validates a run before it starts | pxp, mxp directory, time index | `preflight_table.tsv` | - *(called by `1_submit.sh` with explicit flags)* |
| `1_submit.sh` **(GPU)** | first `cuskss` run, one job per LD block | mxm, mxp, pxp, blockfile, time index | per-block cusk output under `out/{setup}/{runid}/` | `CI_GWAS_ROOT`, `MXM_ROOT`, `DATA_ROOT`, `CIGWAS_BIN`; array size; `ALPHA_EXP`, `MAX_LEVEL`, `MAX_LEVEL_TWO`, `NUM_SAMPLES`, `DEPTH` |
| `2_submit_merge_blocks.sh` | `ci-gwas.py merge-block-outputs` - merges the per-block graphs | per-block output, autosome blockfile | `merged_blocks.*` | `CI_GWAS_ROOT`, `CIGWAS_BIN`, `GEN_SETUPS` |
| `calc_mxm_merged.py` | Recomputes LD across the surviving markers only (PLINK `--r triangle bin4` on the merged marker set) | `merged_blocks_{scm,sam}.mtx`, `merged_blocks.ixs`, pxp, bedfiles | `merged_blocks.mxm.ld.bin`, `merged_blocks_selected_markers.tsv`, `merged_blocks.rsids` | `path/to/bedfiles/autosomes[.bim]`, `path/to/plink`; two CLI arguments: `CUSK_DIR PXP_PATH` |
| `3_submit_mxm_merge.sh` | SLURM wrapper for the above | as above | as above | `CI_GWAS_ROOT`, `GEN_SETUPS` |
| `4_submit_after_merge.sh` **(GPU)** | second `cuskss` run | `merged_blocks.mxm.ld.bin`, `merged_blocks.ixs`, merged mxp, pxp, time index | `cuskss_merged.*` | `CI_GWAS_ROOT`, `CIGWAS_BIN`, `GEN_SETUPS`, `MAX_LEVEL_TWO`, `NUM_SAMPLES` |
| `create_table_cuskss.py` | builds the final results table | `cuskss_merged_{scm,sam,szm}.mtx`, `cuskss_merged.ixs`, pxp, `effective_n.json`, `.bim` | **`cuskss_selected_markers.tsv`** (the *CI-GWAS results table*), `cuskss.rsids` | `path/to/bedfiles/autosomes.bim`, `path/to/mxp/main`; two CLI arguments: `CUSK_DIR PXP_PATH` | 
| `5_make_final_tables.sh` | SLURM wrapper for the above| as above | as above | `CI_GWAS_ROOT`, `GEN_SETUPS` |

## 5. standard-gwas


| File | What it does | Inputs | Outputs | You must set |
| --- | --- | --- | --- | --- |
| `prep_regenie_inputs.py` | builds the run table (one row per trait x scope x covariate set) | *pre/post tables*, covariates, `.fam` | `pairs_all.tsv`, per-run phenotype and covariate files | - *(paths passed in from the notebook)* |
| `regenie-pheno-prep.ipynb` | runs the above | *pre/post tables*, UKB covariates, `.fam` | as above | `<BP_POOLED_TSV>`, `<BP_AGE5_TSV>`, `<LDL_POOLED_TSV>`, `<LDL_AGE5_TSV>`, `<UKB_COVARS_COV>`, `<UKB_COVARS_NAMECOV>`, `<FAM_FILE>`, `<OUT_ROOT>` |
| `regenie_step_1_run.sh` | Regenie step 1 | `pairs_all.tsv`, step-1 bedfiles | `{run}_step1*` | `<BED_STEP1>`, `<OUTROOT>`, conda environment name, array size |
| `regenie_step_2_normal.sh` | Regenie step 2 | step-1 output, step-2 bedfiles | per-chromosome association results | `<BED_STEP2>`, `<OUTROOT>`, array size |
| `regenie_step_2_int.sh` | As above, with the genotype x drug interaction | as above | as above | `<BED_STEP2>`, `<OUTROOT>`, array size |
| `gwas_postproc.py` | checks every run for missing chromosomes | step-2 output tree, `pairs_all.tsv` | check and summary frames | - |
| `result-aggregation.ipynb` | collects the step-2 results and compares them with the CI-GWAS hits | step-2 output, *CI-GWAS results table* | comparison tables | `<OUTROOT>`, `<CI_GWAS_FULL_RESULTS_TSV>` |
| `cigwas_regenie_panels_AB.ipynb` | makes CI-GWAS vs Regenie figure | *CI-GWAS results table*, mxp files, step-2 output | figure panels | `<STEP2_ROOT>`, `<CI_GWAS_FULL_RESULTS_TSV>`, `<SBP_MXP_PATH>`, `<SBP1_MXP_PATH>` |

## 6. arb-response


| File | What it does | Inputs | Outputs | You must set |
| --- | --- | --- | --- | --- |
| `data-extraction-rap/preprocessing/extract_from_xlsx.py` | Splits the UKB `all_lkps_maps` workbook into per-coding CSVs | `all_lkps_maps_v4.xlsx` | per-coding CSVs in `../data/lkps/` | *(relative paths; put the workbook in `data/lkps/`)* |
| `data-extraction-rap/preprocessing/get_lists_of_drugs.py` | Builds code lists for each drug class | lookup CSVs | `{class}_{coding}_drug.csv` | - |
| `data-extraction-rap/preprocessing/merge_lists_of_drugs.py` | concatenates the per-class lists | the above | `all_{coding}_drug.csv` | - |
| `data-extraction-rap/data/` | committed outputs of those three scripts, plus the brand-name and drug-type dictionaries | - | - | - |
| `data-extraction-rap/preprocessing/1_create_database.ipynb` | builds the ARB record database on RAP | UKB dispensed datasets, the code lists above | RAP database | `<NOTEBOOKS_DIR>`, RAP database name |
| `data-extraction-rap/preprocessing/2_add_dos_and_tabs.ipynb` | Adds dose and tablet-count information to the prescription records | the above | annotated prescriptions | `<NOTEBOOKS_DIR>` |
| `drug-response-pheno-prep.ipynb` | Selects the SNPs to carry forward and builds the drug-response phenotype | `5_age_groups.csv`, `combined_SBP.csv` | `snps_to_keep.txt`, response phenotype table | input filenames in the first cells *(relative)* |
| `plink_filter.slurm` | Extracts those SNPs from the merged bedfiles | `snps_to_keep.txt`, `autosomes.*` | `filtered_for_arb_response.*` | bedfile stem and output name |
| `arb-analysis-bp-snps.ipynb` | Association tests of each BP SNP against each drug-response phenotype | filtered bedfiles, response phenotypes | `agg_arb_bwith_bp_snps_and_bp.csv` | `<OUTPUT_DIR>` |
| `kcnip4-wgs.ipynb` | *KCNIP4* follow-up in the DRAGEN population-level WGS data | UKB WGS pVCFs (chr4) | region MatrixTable and variant summaries | `<PROJECT_FOLDER>` |

## 7. ukb-analysis-figures


| Notebook | Produces | Main inputs | You must set |
| --- | --- | --- | --- |
| `table-1-snp-counts-figure-3-a-b.ipynb` | Table 1 (SNP counts) and Figure 3a-b (stacked bars, 50 kb Venns) | *CI-GWAS results table*, pruned and unpruned CI-GWAS output | `<PATH_TO_CI_GWAS_OUT>`, `<PATH_TO_CI_GWAS_OUT_PRUNED>`, `<PATH_TO_CUSKSS_SELECTED_MARKERS_TSV>`, `<PATH_TO_MAIN_SETUP_RESULTS_TSV>`, plus four output paths |
| `figure-3-c.ipynb` | Figure 3c (enrichment balloon plot) and the two enrichment supplementary tables | Enrichr gene-set libraries (GO, Reactome, ENCODE/ChEA, GWAS Catalog), enrichment results | `<PATH_TO_ENRICHMENT_DIR>`, four library paths, three output paths |
| `figure-4-a-c-and-supp.ipynb` | Figures 4a and 4c (drug-target and Open Targets overlap) | *CI-GWAS results table* per setup, ARB SNP and LoF tables, Open Targets exports | eleven `<PATH_TO_*>` inputs plus `<PATH_TO_STACKED_SNP_LOF_PNG>` |
| `figure-4-b.ipynb` | Figure 4b | BP table with zero-day records, drug class table, PLINK genotypes | `<PATH_TO_BP_WITH_ZERO_TSV>`, `<PATH_TO_DRUG_CLASSES_TSV>`, `<PATH_TO_PLINK_ROOT>`, plus two output paths |
| `figure-4-d.ipynb` | Figure 4d (drug-response GLM bar plots) | drug-response GLM input table | `<PATH_TO_FOR_GLM_DRUG_RESPONSE_TSV>` plus two output paths |
| `figure-supp-gwas.ipynb` | Supplementary figure: CI-GWAS vs Regenie | *CI-GWAS results table*, Regenie step-2 output, mxp files | `<PATH_TO_CI_GWAS_FULL_RESULTS_TSV>`, `<PATH_TO_REGENIE_STEP2_ROOT>`, `<PATH_TO_SBP_MXP_TSV>`, `<PATH_TO_SBP_AGE5_WITH_STATINS_MXP_TSV>`, `<PATH_TO_CACHE_STEP2_SIG_DIR>`, `<PATH_TO_FIGURE_PANELS_SVG>` |
| `figure-sup-external-gwas.ipynb` | Supplementary figure: overlap with published BP GWAS | *CI-GWAS results table*, external SBP/DBP summary statistics, LD-proxy SNPs | four `<PATH_TO_*>` inputs plus `<PATH_TO_EXTERNAL_GWAS_PNG>` |
| `figure-supp-gene-based-venn.ipynb` | Supplementary gene-level Venn diagrams | gene lists per setup | `<PATH_TO_GENE_LISTS_TSV>` plus four output paths |
| `table-sup-samplesize.ipynb` | Supplementary sample-size table | *pre/post tables* (pooled and age-split), merged bedfiles | `<PATH_TO_BP_PRE_POST_TSV>`, `<PATH_TO_BP_PRE_POST_AGE5_WITH_STATINS_TSV>`, `<PATH_TO_MERGED_BEDFILES_DIR>`, `<PATH_TO_SAMPLE_SIZE_CSV>` |
| `table-sup-ldl.ipynb` | Supplementary LDL results table | CI-GWAS LDL output | `<PATH_TO_CI_GWAS_OUT>`, `<PATH_TO_CUSKSS_SELECTED_MARKERS_TSV>` |
| `table-sup-interaction-for-post.ipynb` | Supplementary genotype x treatment interaction table, post-treatment traits | *CI-GWAS results table*, Regenie interaction output, PLINK genotypes | `<PATH_TO_CI_GWAS_FULL_RESULTS_TSV>`, `<PATH_TO_REGENIE_STEP2_ROOT>`, `<PATH_TO_REGENIE_PHENO_ROOT>`, `<PATH_TO_PLINK_ROOT>`, plus two output paths |
| `table-supp-interaction-drugs.ipynb` | As above, for the drug-class traits | as above | the same six `<PATH_TO_*>` |
| `table-sup-interaction-LDL.ipynb` | As above, for LDL | as above | eight `<PATH_TO_*>` |

---

