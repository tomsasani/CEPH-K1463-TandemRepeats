# Analysis of tandem repeats in the 4-gen CEPH K1463 family

Guide to files and scripts for reproducing numbers and figures:


### Unfiltered TRGT genotype CSVs:

```
/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/{ASSEMBLY}_v1.0_50bp_merge/493ef25/trgt-denovo/
```

...where `{ASSEMBLY}` is either `GRCh38` or `CHM13v2`.

### Raw loci annotations:

```
/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/{ASSEMBLY}_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed.gz
```

...where {ASSEMBLY} is either `human_GRCh38_no_alt_analysis_set` or `chm13v2.0_maskedY_rCRS` (note the lowercase).

### Minimally-filtered *de novo* calls in G3:

```
/scratch/ucgd/lustre-work/quinlan/u1006375/CEPH-K1463-TandemRepeats/tr_validation/data/denovos/{SAMPLE}.{ASSEMBLY}.prefiltered.tsv
```

The following filters were applied to the [unfiltered calls](#unfiltered-trgt-genotype-csvs):

* `denovo_coverage >= 1`


### Minimally-filtered *de novo* calls in G3, merged with transmission, grandparental evidence, and phasing:

```
/scratch/ucgd/lustre-work/quinlan/u1006375/CEPH-K1463-TandemRepeats/tr_validation/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.tsv
```

#### Additional columns and expected values

These files have a lot of extra columns, but the ones to pay attention to are:

* `children_with_denovo_allele`
    * only present in files for `2189` and `2216`
    * comma-separated list of sample IDs that likely inherited the *de novo* allele length (taken directly from Tom Mokveld)
* `gp_ev`
    * pipe-separated list of evidence for the *de novo* AL in grandparents (`pgf` is paternal grandfather, `pgm` is paternal grandmother, etc.). `N` means no evidence for the *de novo* AL, `Y`  means there was evidence for the *de novo* AL in that grandparent, and `?` means we didn't have a genotype in that grandparent
* `phase_summary`
    * for *de novos* we were able to phase using HiPhase output, this column will be something like `dad:47` or `mom:19`, indicating the likely parent-of-origin as well as the number of informative sites surrounding the *de novo* that were consistent with the specified parent-of-origin

### Orthogonal validation of minimally-filtered calls in G3:

```
/scratch/ucgd/lustre-work/quinlan/u1006375/CEPH-K1463-TandemRepeats/tr_validation/data/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.read_support.csv
```

This file will contain orthogonal validation data **only for the minimally-filtered *de novo* calls with at least 10 orthogonal sequencing reads in mom, dad, and the kid**. 

#### Additional columns and expected values

* `exp_allele_diff_denovo`
    * expected difference (in base pairs) between the *de novo* allele length (defined by TRGT-dn) and the length of the reference allele
    * for example, if the TR locus in the reference is `GAGA` and the *de novo* allele is `GAGAGA`, this column will be `2`
* `exp_allele_diff_non_denovo`
    * expected difference (in base pairs) between the non-*de novo* allele length (defined by TRGT-dn) and the length of the reference allele
* `mom_evidence`
    * pipe-separated list of read evidence for the observed "diffs" in the orthogonal sequencing reads from the mom in the trio
    * for example, this column might contain something like this: `-16:10|0:6`
        * this indicates that there were 10 reads supporting a "diff" of `-16` (that is, 10 reads contained a total of 16 deleted bp w/r/t the reference sequence) and there were 6 reads supporting a "diff" of `0` (that is, 6 reads perfectly matched the reference sequence)
        * you can compare this column to the `exp_allele_diff` columns above to figure out how many reads in this individual supported the expected *de novo* and non-*de novo* "diffs"
* `dad_evidence`
    * same as above, but for dad's reads
* `kid_evidence`
    * same as above, but for the kid's reads




