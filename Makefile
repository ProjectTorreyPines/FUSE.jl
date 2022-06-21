JULIA_PKG_REGDIR ?= $(HOME)/.julia/registries
JULIA_PKG_DEVDIR ?= $(HOME)/.julia/dev
CURRENTDIR = $(shell pwd)

define clone_update_repo
    if [ ! -d "$(JULIA_PKG_DEVDIR)" ]; then mkdir -p $(JULIA_PKG_DEVDIR); fi
	cd $(JULIA_PKG_DEVDIR);\
	if [ ! -d "$(JULIA_PKG_DEVDIR)/$(1)" ]; then git clone --single-branch git@github.com:ProjectTorreyPines/$(1).jl.git $(1) ; else cd $(1) && git pull origin `git rev-parse --abbrev-ref HEAD` ; fi
endef

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
Pkg.develop(["IMAS", "IMASDD", "CoordinateConventions", "MillerExtendedHarmonic", "FusionMaterials", "VacuumFields", "Equilibrium", "TAUENN", "EPEDNN", "TGLFNN", "QED", "FiniteElementHermite", "CHEASE"]);\
Pkg.activate();\
Pkg.develop(["FUSE", "IMAS", "IMASDD", "CoordinateConventions", "MillerExtendedHarmonic", "FusionMaterials", "VacuumFields", "Equilibrium", "TAUENN", "EPEDNN", "TGLFNN", "QED", "FiniteElementHermite", "CHEASE"]);\
'

develop: clone_update_all develop_no_registry precompile

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
Pkg.add(["Revise", "JuliaFormatter", "Test", "Plots", "IJulia"]);\
Pkg.build("IJulia");\
'
	python3 -m pip install --upgrade webio_jupyter_extension

precompile:
	julia -e 'using Pkg; Pkg.activate("."); Pkg.precompile()'

clone_update_all:
	make -j 100 FUSE IMAS IMASDD CoordinateConventions MillerExtendedHarmonic FusionMaterials VacuumFields Equilibrium TAUENN EPEDNN TGLFNN QED CHEASE FiniteElementHermite

update: develop clone_update_all precompile

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

CHEASE:
	$(call clone_update_repo,$@)

FiniteElementHermite:
	$(call clone_update_repo,$@)

docker_image:
	rm -rf ../Dockerfile
	cp docker/Dockerfile ..
	cp .gitignore ../.dockerignore
	cd .. ; sudo docker build -t julia_fuse .

docker_fresh:
	rm -rf ../Dockerfile
	cp docker/Dockerfile_fresh ../Dockerfile
	cp .gitignore ../.dockerignore
	cd .. ; sudo docker build --no-cache --progress=plain -t julia_fuse .

docker_volume:
	docker volume create FUSE
	docker run --rm -v FUSE:/root julia_fuse mkdir -p /root/.julia
	docker run --rm -v FUSE:/root julia_fuse rm -rf /root/.julia/dev
	docker container create --name temp -v FUSE:/root alpine
	docker cp $(HOME)/.julia/dev/. temp:/root/.julia/dev
	#docker cp $(HOME)/.julia/dev/FUSE/Makefile temp:/root/.julia/dev/FUSE/Makefile
	docker rm temp
	docker run --rm -v FUSE:/root julia_fuse bash -c "cd /root/.julia/dev/FUSE && make develop_no_registry && make precompile"

docker_run:
	docker run -it --rm julia_fuse julia

docker_run_volume:
	docker run -it --rm -v FUSE:/root julia_fuse

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

.PHONY:
