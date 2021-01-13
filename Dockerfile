# Base OS
FROM continuumio/miniconda3:4.9.2

ARG VERSION=0.1.6
ARG SIMULATOR_VERSION="0.9.9"
ARG PYTHON_VERSION=3.7

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
    base_image="continuumio/miniconda3:4.9.2" \
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

# Install requirements
ENV CONDA_ENV=py \
    PATH=/opt/conda/envs/py/bin:$PATH \
    MPLBACKEND=PDF
RUN conda update -y -n base -c defaults conda \
    && conda create -y -n ${CONDA_ENV} python=${PYTHON_VERSION} \
    && conda install --name ${CONDA_ENV} -y -c conda-forge \
        assimulo \
    && /bin/bash -c "source activate ${CONDA_ENV}" \
    \
    && apt-get update -y \
    && apt-get install -y --no-install-recommends \
        gfortran \
    && pip install python-libsbml \
    && pip install matplotlib \
    \
    && cd /tmp \
    && git clone --branch assimulo https://github.com/PySCeS/pysces.git \
    && cd pysces \
    && pip install . \
    \
    && cd /tmp \
    && rm -rf pysces \
    \
    && apt-get remove -y \
        gfortran \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists

# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_PySCeS
RUN pip install /root/Biosimulators_PySCeS \
    && mkdir /Pysces \
    && mkdir /Pysces/psc \
    && mv /root/Biosimulators_PySCeS/.pys_usercfg.Dockerfile.ini /Pysces/.pys_usercfg.ini \
    && rm -rf /root/Biosimulators_PySCeS

# supported environment variables
ENV ALGORITHM_SUBSTITUTION_POLICY=SIMILAR_VARIABLES \
    VERBOSE=0

# Entrypoint
ENTRYPOINT ["pysces"]
CMD []
