&nbsp;
&nbsp;
&nbsp;
<p align="center">
🚧🚨 <em>Under development</em> 🚨🚧
</p>

# HMFtools

> Currently only the GRIDSS/PURPLE/LINX TN workflow is supported.

## Table of contents

* [Installation](#installation)
* [Usage](#usage)
* [Outputs](#outputs)
* [Run tests](#run-tests)
* [Requirements](#requirements)
* [Reference data](#reference-data)
* [License](#license)

## Installation

```bash
git clone https://github.com/scwatts/nextflow_pipelines && cd ./nextflow_pipelines/hmftools/
```

## Usage

Configure inputs and pipeline (see [`nextflow.config`](nextflow.config) for example) then execute:

```bash
# Pull in reference data
git clone https://github.com/umccr/reference_data -b dev reference_data_gitrepo/ && cd reference_data_gitrepo/
dvc pull reference_data/{genomes,hmftools}/ -r storage-s3 && cd ../
ln -s reference_data_gitrepo/reference_data/ reference_data

# Install modules
mkdir -p ./modules/
nf-core modules install
nf-core modules -g scwatts/nextflow_modules -b main install

# Launch pipeline
nextflow run ./main.nf -profile docker --input /path/to/samplesheet.tsv --outdir ./output/
```

## Outputs

### Directories

| Name                  | Contents                              |
| ---                   | ---                                   |
| `gridss/`             | GRIDSS final output                   |
| `gripss/`             | Filtered SVs                          |
| `linx_annotation/`    | Data for clustered and annotated SVs  |
| `linx_visualiser/`    | Plots for clustered and annotated SVs |
| `nextflow/`           | Pipeline config, logs, and reports    |
| `nextflow/work/`      | Intermediary work files               |
| `purple/`             | CNA calls, purity, ploidy             |

### Useful files

| Name                                      | Description                               |
| ---                                       | ---                                       |
| `<tumor_name>_linx.html`                  | gpgr LINX Rmd report                      |
| `linx_annotation/*tsv`                    | Grouped and annotated SV events           |
| `linx_visualiser/plot/*png`               | SV event plots                            |
| `purple/plot/*png`                        | Purity, ploidy, circos, etc plots         |
| `purple/<tumor_name>.<vcf_type>.vcf.gz`   | VCF provided to and annotated by PURPLE   |
| `gridss/sv_annotated.vcf.gz`              | Final GRIDSS output SV VCF                |
| `gripps/<prefix>.gripps.filtered.vcf.gz`  | Hard filtered SV VCF                      |
| `gripps/<prefix>.gripss.vcf.gz`           | Soft filtered SV VCF                      |
| `nextflow/nextflow_log.txt`               | Pipeline log file                         |
| `nextflow/nextflow.config`                | Pipeline configuration used in run        |
| `nextflow/reports/timeline.html`          | Stage execution durations as a timeline   |

## Run tests

```bash
# Pull in reference data
git clone https://github.com/umccr/reference_data -b dev reference_data_gitrepo/ && cd reference_data_gitrepo/
dvc pull reference_data/{genomes,hmftools}/ -r storage-s3 && cd ../
ln -s reference_data_gitrepo/reference_data/ reference_data

# Install modules
mkdir -p ./modules/
nf-core modules install
nf-core modules -g scwatts/nextflow_modules -b main install

# Run test
nextflow run ./main.nf -profile docker,test --outdir output/
```

## Requirements

> Software versions only indicate what is currently in use rather than strict requirements

### Pipeline

* [Nextflow](https://www.nextflow.io/) (v22.04.3)
* [nf-core](https://nf-co.re) (v2.4.1)
* [Docker](https://www.docker.com/get-started) (v20.10.11)
* [BCFtools](https://www.htslib.org/) (v1.15.1)

### GPL toolkit

* [AMBER](https://github.com/hartwigmedical/hmftools/blob/master/amber/) (v3.9)
* [COBALT](https://github.com/hartwigmedical/hmftools/blob/master/cobalt/) (v1.13)
* [GRIDSS](https://github.com/PapenfussLab/gridss) (v2.13.2)
* [GRIPSS](https://github.com/hartwigmedical/hmftools/blob/master/gripss/) (v2.1)
* [PURPLE](https://github.com/hartwigmedical/hmftools/blob/master/purple/) (v3.4)
* [LINX](https://github.com/hartwigmedical/hmftools/blob/master/linx/) (v1.19)

## License

Software and code in this repository are under [GNU General Public License
v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) unless otherwise indicated.
