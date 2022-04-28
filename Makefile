JULIA_PKG_REGDIR ?= $(HOME)/.julia/registries
JULIA_PKG_DEVDIR ?= $(HOME)/.julia/dev
CURRENTDIR = $(shell pwd)

all:
	@echo 'FUSE makefile help'
	@echo ''
	@echo ' - make develop  : install FUSE and its PTP dependencies to $(JULIA_PKG_DEVDIR)'
	@echo ' - make update   : git pull FUSE and its PTP dependencies'
	@echo ''

registry:
	julia -e 'using Pkg;Pkg.add("Revise")' # call this first to make sure General registry gets installed
	if [ ! -d "$(JULIA_PKG_REGDIR)" ]; then mkdir -p $(JULIA_PKG_REGDIR); fi
	cd $(JULIA_PKG_REGDIR);\
	if [ ! -d "$(JULIA_PKG_REGDIR)/GAregistry" ]; then git clone git@github.com:ProjectTorreyPines/GAregistry.git GAregistry ; fi

develop_no_registry:
	julia -e '\
using Pkg;\
Pkg.activate(".");\
Pkg.develop(["IMAS", "IMASDD", "CoordinateConventions", "FusionMaterials", "VacuumFields", "Equilibrium", "TAUENN", "EPEDNN", "TGLFNN", "QED", "FiniteElementHermite"]);\
Pkg.activate();\
Pkg.develop(["FUSE", "IMAS", "IMASDD", "CoordinateConventions", "FusionMaterials", "VacuumFields", "Equilibrium", "TAUENN", "EPEDNN", "TGLFNN", "QED", "FiniteElementHermite"]);\
'

develop: registry develop_no_registry

IJulia:
	julia -e '\
using Pkg;\
Pkg.add("IJulia");\
Pkg.build("IJulia");\
'

precompile:
	julia -e 'using Pkg; Pkg.activate("."); Pkg.precompile()'

update_all:
	make -j 100 update_FUSE update_IMAS update_IMASDD update_CoordinateConventions update_FusionMaterials update_VacuumFields update_Equilibrium update_TAUENN update_EPEDNN update_TGLFNN update_QED update_FiniteElementHermite

update: develop update_all precompile

update_FUSE:
	cd $(JULIA_PKG_DEVDIR)/FUSE; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

update_IMAS:
	cd $(JULIA_PKG_DEVDIR)/IMAS; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

update_IMASDD:
	cd $(JULIA_PKG_DEVDIR)/IMASDD; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

update_CoordinateConventions:
	cd $(JULIA_PKG_DEVDIR)/CoordinateConventions; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

update_FusionMaterials:
	cd $(JULIA_PKG_DEVDIR)/FusionMaterials; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

update_VacuumFields:
	cd $(JULIA_PKG_DEVDIR)/VacuumFields; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

update_Equilibrium:
	cd $(JULIA_PKG_DEVDIR)/Equilibrium; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

update_TAUENN:
	cd $(JULIA_PKG_DEVDIR)/TAUENN; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

update_TGLFNN:
	cd $(JULIA_PKG_DEVDIR)/TGLFNN; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

update_EPEDNN:
	cd $(JULIA_PKG_DEVDIR)/EPEDNN; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

update_QED:
	cd $(JULIA_PKG_DEVDIR)/QED; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

update_FiniteElementHermite:
	cd $(JULIA_PKG_DEVDIR)/FiniteElementHermite; git fetch; git pull; julia -e 'using Pkg; Pkg.activate("."); Pkg.resolve()'

docker_image:
	cd docker; sudo docker build -t julia_fuse .

docker_volume:
	docker volume create FUSE

docker_data:
	docker run --rm -v FUSE:/root julia_fuse mkdir -p /root/.julia
	docker run --rm -v FUSE:/root julia_fuse rm -rf /root/.julia/dev
	docker container create --name temp -v FUSE:/root alpine
	docker cp $(HOME)/.julia/dev/. temp:/root/.julia/dev
	#docker cp $(HOME)/.julia/dev/FUSE/Makefile temp:/root/.julia/dev/FUSE/Makefile
	docker rm temp
	docker run --rm -v FUSE:/root julia_fuse bash -c "cd /root/.julia/dev/FUSE && make develop_no_registry && make precompile"

docker_run:
	docker run -it --rm -v FUSE:/root julia_fuse

.PHONY:
