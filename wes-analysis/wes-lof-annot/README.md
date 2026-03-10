# Installation

```
# pip install wheel
pip install git+https://github.com/gosborcz/ukb-bp-rep#subdirectory=wes_lof_annot
```

# Usage

## Annotate VCF files with loftee

Chromosome 6
```
install_vep
annotate_vcf 6
```

## Aggregate loss-of-function variants into csv file

```
rare_variants_table 6 [path/to/eids]
```
