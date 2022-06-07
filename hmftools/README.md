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
* [Requirements](#requirements)
* [Reference data](#reference-data)
* [License](#license)

## Installation

```bash
git clone https://github.com/scwatts/nextflow_pipelines && cd ./nextflow_pipelines/hmftools/
```

## Usage

First you'll need to obtain reference data as described [here](#reference-data). Then create a configuration file (for
an example see: [`nextflow.config`](nextflow.config)). To execute the pipeline:

```bash
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

## Requirements

> Software versions only indicate what is currently in use rather than strict requirements

### Pipeline

* [Docker](https://www.docker.com/get-started) (v20.10.11)
* [Nextflow](https://www.nextflow.io/) (v22.04.3)
* [BCFtools](https://www.htslib.org/) (v1.15.1)

### GPL toolkit

* [AMBER](https://github.com/hartwigmedical/hmftools/blob/master/amber/) (v3.9)
* [COBALT](https://github.com/hartwigmedical/hmftools/blob/master/cobalt/) (v1.13)
* [GRIDSS](https://github.com/PapenfussLab/gridss) (v2.13.2)
* [GRIPSS](https://github.com/hartwigmedical/hmftools/blob/master/gripss/) (v2.1)
* [PURPLE](https://github.com/hartwigmedical/hmftools/blob/master/purple/) (v3.4)
* [LINX](https://github.com/hartwigmedical/hmftools/blob/master/linx/) (v1.19)

## Reference data

The GPL toolkit requires a number of reference files. These can be obtained from the HMF Nextcloud instance
[here](https://nextcloud.hartwigmedicalfoundation.nl/s/LTiKTd8XxBqwaiC?path=%2FHMFTools-Resources). Alternatively, I've
precompiled the required files on S3, located at `s3://umccr-refdata-dev/gpl-nf/`.

## License

Software and code in this repository are under [GNU General Public License
v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) unless otherwise indicated.
