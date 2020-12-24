# Base OS
FROM continuumio/miniconda3:4.9.2

ARG VERSION=0.0.1
ARG SIMULATOR_VERSION="0.9.8"

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
ENV PYTHON_VERSION=3.8 \
    CONDA_ENV=py3 \
    PATH=/opt/conda/envs/${CONDA_ENV}/bin:${PATH}
RUN conda update -y -n base -c defaults conda \
    && conda create -y -n ${CONDA_ENV} python=${PYTHON_VERSION} \
    && conda install --name ${CONDA_ENV} -y -c bgoli -c conda-forge \
        python-libsbml \
        pysces==${SIMULATOR_VERSION} \
    && conda activate ${CONDA_ENV} \
    && apt-get update -y \
    && apt-get install -y --no-install-recommends \
        gcc \
        build-essential \
    && pip install pysundials \
    && apt-get remove -y \
        gcc \
        build-essential \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists


ARG sundials_version=2.3.0
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        gfortran \
        libopenblas-base \
        libopenblas-dev \
        wget \
    \
    && cd /tmp \
    && wget https://computation.llnl.gov/projects/sundials/download/sundials-${sundials_version}.tar.gz \
    && tar xzf sundials-${sundials_version}.tar.gz \
    && cd sundials-${sundials_version} \
    && mkdir build \
    && cd build \
    && cmake \
        -DEXAMPLES_ENABLE=OFF \
        -DLAPACK_ENABLE=ON \
        -DSUNDIALS_INDEX_TYPE=int32_t \
        .. \
    && make \
    && make install \
    \
    && cd /tmp \
    && rm sundials-${sundials_version}.tar.gz \
    && rm -r sundials-${sundials_version} \
    \
    && apt-get remove -y \
        build-essential \
        cmake \
        libopenblas-dev \
        wget \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*


# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_pysces
RUN pip3 install /root/Biosimulators_pysces

# Entrypoint
ENTRYPOINT ["pysces"]
CMD []
