[![Latest release](https://img.shields.io/github/v/tag/biosimulators/Biosimulators_PySCeS)](https://github.com/biosimulations/Biosimulators_PySCeS/releases)
[![PyPI](https://img.shields.io/pypi/v/biosimulators_pysces)](https://pypi.org/project/biosimulators_pysces/)
[![CI status](https://github.com/biosimulators/Biosimulators_PySCeS/workflows/Continuous%20integration/badge.svg)](https://github.com/biosimulators/Biosimulators_PySCeS/actions?query=workflow%3A%22Continuous+integration%22)
[![Test coverage](https://codecov.io/gh/biosimulators/Biosimulators_PySCeS/branch/dev/graph/badge.svg)](https://codecov.io/gh/biosimulators/Biosimulators_PySCeS)
[![All Contributors](https://img.shields.io/github/all-contributors/biosimulators/Biosimulators_PySCeS/HEAD)](#contributors-)

# BioSimulators-PySCeS
BioSimulators-compliant command-line interface to the [PySCeS](http://pysces.sourceforge.net/) simulation program.

This command-line interface and Docker image enable users to use PySCeS to execute [COMBINE/OMEX archives](https://combinearchive.org/) that describe one or more simulation experiments (in [SED-ML format](https://sed-ml.org)) of one or more models (in [SBML format](http://sbml.org])).

A list of the algorithms and algorithm parameters supported by PySCeS is available at [BioSimulators](https://biosimulators.org/simulators/pysces).

A simple web application and web service for using PySCeS to execute COMBINE/OMEX archives is also available at [runBioSimulations](https://run.biosimulations.org).

## Installation

### Install Python package
```
pip install biosimulators-pysces
```

### Install Docker image
```
docker pull ghcr.io/biosimulators/pysces
```

## Usage

### Local usage
```
usage: biosimulators-pysces [-h] [-d] [-q] -i ARCHIVE [-o OUT_DIR] [-v]

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

### Usage through Docker container
The entrypoint to the Docker image supports the same command-line interface described above.

For example, the following command could be used to use the Docker image to execute the COMBINE/OMEX archive `./modeling-study.omex` and save its outputs to `./`.

```
docker run \
  --tty \
  --rm \
  --mount type=bind,source="$(pwd)",target=/root/in,readonly \
  --mount type=bind,source="$(pwd)",target=/root/out \
  ghcr.io/biosimulators/pysces:latest \
    -i /root/in/modeling-study.omex \
    -o /root/out
```

## Documentation
Documentation is available at https://docs.biosimulators.org/Biosimulators_PySCeS/.

## License
This package is released under the [MIT license](LICENSE).

## Development team
This package was developed by the [Center for Reproducible Biomedical Modeling](http://reproduciblebiomodels.org) and the [Karr Lab](https://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York with assistance from the contributors listed [here](CONTRIBUTORS.md).

## Questions and comments
Please contact the [BioSimulators Team](mailto:info@biosimulators.org) with any questions or comments.
