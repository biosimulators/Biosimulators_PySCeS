# Base OS
FROM ghcr.io/biosimulators/biosimulators_pysces/pysces_base:latest

ARG VERSION=0.1.18
ARG SIMULATOR_VERSION="0.9.9"

# metadata
LABEL \
    org.opencontainers.image.title="PySCeS" \
    org.opencontainers.image.version="${SIMULATOR_VERSION}" \
    org.opencontainers.image.description="Simulation and analysis tools for modelling biological systems" \
    org.opencontainers.image.url="http://pysces.sourceforge.net/" \
    org.opencontainers.image.documentation="https://pythonhosted.org/PySCeS/" \
    org.opencontainers.image.source="https://github.com/biosimulators/Biosimulators_PySCeS" \
    org.opencontainers.image.authors="BioSimulators Team <info@biosimulators.org>" \
    org.opencontainers.image.vendor="BioSimulators Team" \
    org.opencontainers.image.licenses="BSD-3-Clause" \
    \
    base_image="python:3.9-slim-buster" \
    version="${VERSION}" \
    software="PySCeS" \
    software.version="${SIMULATOR_VERSION}" \
    about.summary="Simulation and analysis tools for modelling biological systems" \
    about.home="http://pysces.sourceforge.net/" \
    about.documentation="https://pythonhosted.org/PySCeS/" \
    about.license_file="https://github.com/PySCeS/pysces/blob/master/LICENSE.txt" \
    about.license="SPDX:BSD-3-Clause" \
    about.tags="BioSimulators,mathematical model,kinetic model,simulation,systems biology,computational biology,SBML,SED-ML,COMBINE,OMEX" \
    extra.identifiers.biotools="pysces" \
    maintainer="BioSimulators Team <info@biosimulators.org>"

# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_PySCeS
RUN pip install /root/Biosimulators_PySCeS \
    && mkdir /Pysces \
    && mkdir /Pysces/psc \
    && mv /root/Biosimulators_PySCeS/.pys_usercfg.Dockerfile.ini /Pysces/.pys_usercfg.ini \
    && rm -rf /root/Biosimulators_PySCeS

# supported environment variables
ENV ALGORITHM_SUBSTITUTION_POLICY=SIMILAR_VARIABLES \
    VERBOSE=0 \
    MPLBACKEND=PDF

# Entrypoint
ENTRYPOINT ["biosimulators-pysces"]
CMD []
