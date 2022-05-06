FROM julia:latest
RUN apt-get -y update
RUN apt-get -y install --no-install-recommends git make
RUN apt-get autoclean

WORKDIR /root/.julia/dev
COPY ./CoordinateConventions CoordinateConventions/
COPY ./EPEDNN EPEDNN/
COPY ./Equilibrium Equilibrium/
COPY ./FUSE FUSE/
COPY ./FiniteElementHermite FiniteElementHermite/
COPY ./FusionMaterials FusionMaterials/
COPY ./IMAS IMAS/
COPY ./IMASDD IMASDD/
COPY ./QED QED/
COPY ./TAUENN TAUENN/
COPY ./TGLFNN TGLFNN/
COPY ./VacuumFields VacuumFields/

WORKDIR /root/.julia/dev/FUSE
RUN cat Makefile
RUN make develop_no_registry
RUN make precompile