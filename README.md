# BioSimulators-PySCeS
BioSimulators-compliant command-line interface to the [PySCeS](http://pysces.sourceforge.net/) simulation program.

This command-line interface and Docker image enable users to use PySCeS to execute [COMBINE/OMEX archives](https://combinearchive.org/) that describe one or more simulation experiments (in [SED-ML format](https://sed-ml.org)) of one or more models (in [SBML format](http://sbml.org])).

A list of the algorithms and algorithm parameters supported by PySCeS is available at [BioSimulators](https://biosimulators.org/simulators/pysces).

A simple web application and web service for using PySCeS to execute COMBINE/OMEX archives is also available at [runBioSimulations](https://run.biosimulations.org).

## Contents
* [Installation](#installation)
* [Usage](#usage)
* [License](#license)
* [Development team](#development-team)
* [Questions and comments](#questions-and-comments)

## Installation

### Install Python package
```
pip install git+https://github.com/biosimulators/Biosimulators_PySCeS
```

### Install Docker image
```
docker pull ghcr.io/biosimulators/pysces
```

## Local usage
```
usage: pysces [-h] [-d] [-q] -i ARCHIVE [-o OUT_DIR] [-v]

BioSimulators-compliant command-line interface to the PySCeS simulation program <http://pysces.sourceforge.net/>.

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           full application debug mode
  -q, --quiet           suppress all console output
  -i ARCHIVE, --archive ARCHIVE
                        Path to OMEX file which contains one or more SED-ML-
                        encoded simulation experiments
  -o OUT_DIR, --out-dir OUT_DIR
                        Directory to save outputs
  -v, --version         show program's version number and exit
```

## Usage through Docker container
```
docker run \
  --tty \
  --rm \
  --mount type=bind,source="$(pwd)"/tests/fixtures,target=/root/in,readonly \
  --mount type=bind,source="$(pwd)"/tests/results,target=/root/out \
  ghcr.io/biosimulators/pysces:latest \
    -i /root/in/BIOMD0000000297.omex \
    -o /root/out
```

## License
This package is released under the [MIT license](LICENSE).

## Development team
This package was developed by the [Center for Reproducible Biomedical Modeling](http://reproduciblebiomodels.org) and the [Karr Lab](https://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York.

## Questions and comments
Please contact the [BioSimulators Team](mailto:info@biosimulators.org) with any questions or comments.