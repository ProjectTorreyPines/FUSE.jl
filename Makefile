JULIA_PKG_REGDIR ?= $(HOME)/.julia/registries
JULIA_PKG_DEVDIR ?= $(HOME)/.julia/dev
CURRENTDIR = $(shell pwd)

FUSE_BROKER := 35.247.84.151

DOCKER_PLATFORM := $(shell uname -m)

define clone_update_repo
    if [ ! -d "$(JULIA_PKG_DEVDIR)" ]; then mkdir -p $(JULIA_PKG_DEVDIR); fi
	cd $(JULIA_PKG_DEVDIR);\
	if [ ! -d "$(JULIA_PKG_DEVDIR)/$(1)" ]; then git clone --single-branch git@github.com:ProjectTorreyPines/$(1).jl.git $(1) ; else cd $(1) && git pull origin `git rev-parse --abbrev-ref HEAD` ; fi
endef

all:
	@echo 'FUSE makefile help'
	@echo ''
	@echo ' - make install      : install FUSE and its PTP dependencies to $(JULIA_PKG_DEVDIR)'
	@echo ' - make update       : git pull FUSE and its PTP dependencies'
	@echo ' - make IJulia       : Install IJulia'
	@echo ' - make dd           : regenerate IMADDD.dd.jl file'
	@echo ' - make docker_image : generate a new FUSE docker image'
	@echo ' - make docker_run   : run the FUSE docker image'
	@echo ''

registry:
	julia -e 'using Pkg;Pkg.add("Revise")' # call this first to make sure General registry gets installed
	if [ ! -d "$(JULIA_PKG_REGDIR)" ]; then mkdir -p $(JULIA_PKG_REGDIR); fi
	cd $(JULIA_PKG_REGDIR);\
	if [ ! -d "$(JULIA_PKG_REGDIR)/GAregistry" ]; then git clone git@github.com:ProjectTorreyPines/GAregistry.git GAregistry ; fi

install_no_registry:
	julia -e '\
using Pkg;\
Pkg.activate(".");\
Pkg.develop(["IMAS", "IMASDD", "CoordinateConventions", "MillerExtendedHarmonic", "FusionMaterials", "VacuumFields", "Equilibrium", "TAUENN", "EPEDNN", "TGLFNN", "QED", "FiniteElementHermite", "Fortran90Namelists", "CHEASE", "EFIT", "NNeutronics", "Broker", "ZMQ"]);\
Pkg.activate();\
Pkg.develop(["FUSE", "IMAS", "IMASDD", "CoordinateConventions", "MillerExtendedHarmonic", "FusionMaterials", "VacuumFields", "Equilibrium", "TAUENN", "EPEDNN", "TGLFNN", "QED", "FiniteElementHermite", "Fortran90Namelists", "CHEASE", "EFIT", "NNeutronics", "Broker", "ZMQ"]);\
'

rm_manifests:
	find .. -name "Manifest.toml" -exec rm -rf \{\} \;

install: clone_update_all install_no_registry precompile

sysimage:
	julia -e '\
using Pkg;\
Pkg.add("PackageCompiler");\
Pkg.add("IJulia");\
import PackageCompiler;\
Pkg.activate(".");\
using FUSE;\
PackageCompiler.create_sysimage(["FUSE"], sysimage_path="FUSEsysimage.so");\
'

sysimage_ijulia:
	julia -e '\
import IJulia;\
IJulia.installkernel("Julia FUSEsysimage", "--sysimage=$(shell pwd)/FUSEsysimage.so", "--trace-compile=stderr");\
'

IJulia:
	julia -e '\
using Pkg;\
Pkg.add(["Revise", "JuliaFormatter", "Test", "Plots", "IJulia", "Interact"]);\
Pkg.build("IJulia");\
'
	python3 -m pip install --upgrade webio_jupyter_extension

precompile:
	julia -e 'using Pkg; Pkg.precompile()'

clone_update_all:
	make -j 100 FUSE IMAS IMASDD CoordinateConventions MillerExtendedHarmonic FusionMaterials VacuumFields Equilibrium TAUENN EPEDNN TGLFNN QED FiniteElementHermite Fortran90Namelists CHEASE EFIT NNeutronics Broker ZMQ

update: install clone_update_all precompile

FUSE:
	$(call clone_update_repo,$@)

IMAS:
	$(call clone_update_repo,$@)

IMASDD:
	$(call clone_update_repo,$@)

CoordinateConventions:
	$(call clone_update_repo,$@)

MillerExtendedHarmonic:
	$(call clone_update_repo,$@)

FusionMaterials:
	$(call clone_update_repo,$@)

VacuumFields:
	$(call clone_update_repo,$@)

Equilibrium:
	$(call clone_update_repo,$@)

TAUENN:
	$(call clone_update_repo,$@)

TGLFNN:
	$(call clone_update_repo,$@)

EPEDNN:
	$(call clone_update_repo,$@)

QED:
	$(call clone_update_repo,$@)

FiniteElementHermite:
	$(call clone_update_repo,$@)

Fortran90Namelists:
	$(call clone_update_repo,$@)

CHEASE:
	$(call clone_update_repo,$@)

EFIT:
	$(call clone_update_repo,$@)

NNeutronics:
	$(call clone_update_repo,$@)

Broker:
	$(call clone_update_repo,$@)

ZMQ:
	$(call clone_update_repo,$@)

docker_network:
	docker network create fuse-net

# build a new FUSE docker image
docker_image:
	rm -rf ../Dockerfile
	cp docker/Dockerfile_base ../Dockerfile
	sed 's/_PLATFORM_/$(DOCKER_PLATFORM)/g' ../Dockerfile > tmp
	mv tmp ../Dockerfile
	cat ../Dockerfile
	cp .gitignore ../.dockerignore
	cd .. ; sudo docker build --platform=linux/$(DOCKER_PLATFORM) -t julia_fuse_$(DOCKER_PLATFORM) .

# update the base FUSE docker image
docker_update:
	rm -rf ../Dockerfile
	cp docker/Dockerfile_update ../Dockerfile
	sed 's/_PLATFORM_/$(DOCKER_PLATFORM)/g' ../Dockerfile > tmp
	mv tmp ../Dockerfile
	cat ../Dockerfile
	cp .gitignore ../.dockerignore
	cd .. ; sudo docker build --platform=linux/$(DOCKER_PLATFORM) -t julia_fuse_$(DOCKER_PLATFORM)_updated .

# run FUSE worker container
broker_run:
	docker run -it --rm --platform=linux/$(DOCKER_PLATFORM) --network=fuse-net -p 6666:6666 -p 7777:7777 -p 8888:8888 -p 9999:9999 --name=broker julia_fuse_$(DOCKER_PLATFORM)_updated python3 Broker/src/broker.py

# run FUSE worker container
worker_run:
	docker run -it --rm --platform=linux/$(DOCKER_PLATFORM) --network=fuse-net julia_fuse_$(DOCKER_PLATFORM)_updated julia -e 'import FUSE; FUSE.worker_start("$(FUSE_BROKER)")'

# test FUSE client container 
client_test:
	docker run -it --rm --platform=linux/$(DOCKER_PLATFORM) --network=fuse-net julia_fuse_$(DOCKER_PLATFORM)_updated julia -e 'import FUSE; FUSE.client_test("$(FUSE_BROKER)")'

cleanup:
	julia -e 'using Pkg; using Dates; Pkg.gc(; collect_delay=Dates.Day(0))'

html:
	cd docs; julia make.jl

web:
	if [ ! -d "$(PWD)/docs/pages" ]; then cd docs; git clone --single-branch -b gh-pages git@github.com:ProjectTorreyPines/FUSE.jl.git pages; fi
	#cd docs/pages; git reset --hard 049da2c703ad7fc552c13bfe0651da677e3c7f58
	cd docs/pages; git pull
	cd docs; cp -rf build/* pages/
	cd docs/pages; touch .nojekyll
	cd docs/pages; git add -A; git commit --allow-empty -m "documentation"; git push --force

examples: .PHONY
	cd docs; julia notebooks_to_html.jl

dd:
	julia ../IMASDD/src/generate_dd.jl

.PHONY:
