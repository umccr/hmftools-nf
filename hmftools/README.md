&nbsp;
&nbsp;
&nbsp;
<p align="center">
🚧🚨 <em>Under development</em> 🚨🚧
</p>

# HMFtools

> TDB

## Table of contents

* [Installation](#installation)
* [Usage](#usage)
* [Outputs](#outputs)
* [Run tests](#run-tests)
* [Requirements](#requirements)
* [Reference data](#reference-data)
* [License](#license)

## Overview

> TDB

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

## Pipeline modes

> TBD

Several modes of execution are offered and can be accessed using the `--mode` argument.

| Name                  | Description                           |
| ---                   | ---                                   |
| `full`                |                                       |
| `gridss`              |                                       |
| `purple`              |                                       |
| `linx`                |                                       |
| `gridss-purple-linx`  |                                       |
| `lilac`               |                                       |
| `teal`                |                                       |

> Currently only T/N inputs are offered


## Outputs

### Directories

> TBD

### Useful files

> TDB

## Run tests
The internal pipeline logic can be tested using a stub run. No actual outputs are generated in a stub run and
only the `stub` section of a process is executed. Hence, it completes in a short amount of time but does not test
validity of actual processes beyond declared inputs and outputs.

```bash
nextflow run main.nf --input ./assets/samplesheet.tsv --outdir output/ --max_memory '1.GB' -stub
```

A more comprehensive test that involves both the process `script` section and internal pipeline logic can be run with
the supplied downsampled test dataset. This of course takes longer (on the order of minutes) but replicates a full run
and generates 'real' process outputs.

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
nextflow run ./main.nf -profile docker,test --outdir output/ --max_memory '6.GB'
```

## Requirements

> TDB

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
