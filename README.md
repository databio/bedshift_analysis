# bedshift analysis

## Contents of the repository

This repository contains results from [Bedshift: perturbation of genomic interval sets](https://doi.org/10.1101/2020.11.11.378554), a paper describing application of the [bedshift](http://bedshift.databio.org) tool to explore the behavior of region set similarity metrics on simulated data.

There are 4 experiment folders. Each of these is organized as a [PEP](http://pep.databio.org) that includes a CSV file containing bedshift parameter sets:

- [pep_demo](pep_demo) provides a quick demo that runs quickly for testing purposes.
- [pep_basic](pep_basic) provides a set of 36 parameter sets to test combinatorial perturbations on a single input file
- [pep_main](pep_main) runs the same parameter sets, but across 3 different input files
- [pep_universe](pep_universe) runs the same parameter sets with the original file from pep_basic, but using 3 different universes

For each experiment, you will find a table of perturbations, with one row per perturbation, and a configuration file. For example, in the `pep_main` project, look at [sample_table.csv](/pep_main/sample_table.csv). The [config file](/pep_main/config.yaml) points to this.

This repository also contains 2 pipeline interface files the [piface_bedshift.yaml](piface_bedshift.yaml), describes how to run `bedshift` on the specified bed file with the listed perturbation parameters; and the [piface_similarity_scores.yaml](piface_similarity_scores.yaml) file runs our included similarity score calculation script (in [src](/src)).


## Setup and configuration

### Refgenie setup

To use refgenie to grab the chrom sizes file, just do:

```
pip install refgenie
export REFGENIE="refgenie_config.yaml"
refgenie init -c $REFGENIE
refgenie pull hg38/fasta
```

### Environment setup

```
export CODE=/path/to/directory
```

Where `directory` is the directory containing this repo, `bedshift_analysis`.

## Demo

### Run bedshift

Start with the PEP in the [pep_demo](/pep_demo) folder. Everything needed to run this is stored in this repository, and it's a short example to show that you have everything set up correctly. Run it like this.

```
looper run pep_demo/project_config.yaml
```

Or, to use bulker, just do `looper run pep_demo/project_config.yaml -p bulker_local`. This will produce an output folder in `results/pep_demo/bedshifted_regions`. 

### Aggregate bedshift scores

After this is complete, you can aggregate the results and generate the scores with `runp` like this:

```
looper runp project_config.yaml
```

Or `looper runp project_config.yaml -p bulker_local`. This will create score files in `pep_demo/results/scores`.

### Make plots

Once the scores are generated, you can reproduce the plots we made in the paper by following the R script in [src/plot_summary_results.R](src/plot_summary_results.R).

## Real project

There are 3 real projects in here:

- `pep_basic` is a basic run, which tests a single BED file against a single universe
- `pep_main` is the main project, which tests 3 BED files, using just 1 universe.
- `pep_universe` tests how different universes behave, using 1 BED file and 3 different universes.

All projects can be run the same way.

### Download data

Run `./src/download_data.sh` to download all the data required to run all projects.

The bed files will be downloaded from bedbase:

- CTCF TF ChIP-seq on human HCT116: http://bedbase.org/bedsplash/713f58a6497a9168a326123919672ebe
- H3K4me3 Histone ChIP-seq on human GM12864 http://bedbase.org/bedsplash/0f84fea95b736ec99914bc66e74ab6e0
- DNase-seq on human stromal cell of bone marrow: http://bedbase.org/bedsplash/c75ea5133f825d779a02be41a529342e

The universes are downloaded from SCREEN:

- https://screen.encodeproject.org/

You can also see more information about the universes at [bedbase.org](http://bedbase.org): 

- 1) GRCh38-ccREs (Primary universe, http://bedbase.org/bedsplash/f31d4aa5e499637f28a338d5768e4ad5)
- 2) DNase-H3K4me3 (http://bedbase.org/bedsplash/87d8e916fc910254aa61d5a1611622b7)
- 3) CTCF-only (http://bedbase.org/bedsplash/fd007af75d7d0ef5c2e76d8e94916dcf)
- 4) PLS (http://bedbase.org/bedsplash/d0b72e9adfeaf3ec37f303f16ad36bc4)


### Run bedshift

 Here are some example looper run commands:

```
looper run pep_main/project_config.yaml -p local
looper run pep_universe/project_config.yaml -p local
looper run pep_universe/project_config.yaml --lumpn 10
looper run pep_universe/project_config.yaml --lumpn 10 --sel-attr trial --sel-incl 0 1 2
```

### Aggregate results

```
looper runp pep_universe/project_config.yaml -l 1 -p local
looper runp pep_main/project_config.yaml
```
