# Base OS
FROM ubuntu:20.04

# metadata
LABEL base_image="ubuntu:20.04"
LABEL version="1.0.0"
LABEL software="PySCeS"
LABEL software.version="0.9.8"
LABEL about.summary="Simulation and analysis tools for modelling biological systems"
LABEL about.home="http://pysces.sourceforge.net/"
LABEL about.documentation="https://pythonhosted.org/PySCeS/"
LABEL about.license_file="https://github.com/PySCeS/pysces/blob/master/LICENSE.txt"
LABEL about.license="BSD-3-Clause"
LABEL about.tags="BioSimulators,mathematical model,kinetic model,simulation,systems biology,computational biology,SBML,SED-ML,COMBINE,OMEX"
LABEL extra.identifiers.biotools="pysces"
LABEL maintainer="BioSimulators Team <info@biosimulators.org>"

# Install requirements
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
    && pip3 install -U pip \
    && pip3 install -U setuptools \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Copy code for command-line interface into image and install it
COPY . /root/biosimulators_pysces
RUN pip3 install /root/biosimulators_pysces

# Entrypoint
ENTRYPOINT ["pysces"]
CMD []
