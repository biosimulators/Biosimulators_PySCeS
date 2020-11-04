# Base OS
FROM continuumio/miniconda3:4.8.2

# metadata
LABEL base_image="continuumio/miniconda3:4.8.2"
LABEL version="0.0.1"
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
ENV CONDA_ENV=py37 \
    PATH=/opt/conda/envs/${CONDA_ENV}/bin:${PATH}
RUN conda update -y -n base -c defaults conda \
    && conda create -y -n ${CONDA_ENV} python=3.7 \
    && conda install --name ${CONDA_ENV} -y -c bgoli -c sbmlteam -c conda-forge \
        pysces \
        python-libsbml
    && /bin/bash -c "source activate ${CONDA_ENV}" \
    && pip install pysundials

# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_pysces
RUN pip3 install /root/Biosimulators_pysces

# Entrypoint
ENTRYPOINT ["pysces"]
CMD []
