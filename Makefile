JULIA_PKG_REGDIR ?= $(HOME)/.julia/registries
JULIA_PKG_DEVDIR ?= $(HOME)/.julia/dev
CURRENTDIR := $(shell pwd)
TODAY := $(shell date +'%Y-%m-%d')

PTP_PACKAGES := $(shell find ../*/.git/config -exec grep ProjectTorreyPines \{\} \; | cut -d'/' -f 2 | cut -d'.' -f 1 | tr '\n' ' ')

# use command line interface for git to work nicely with private repos
export JULIA_PKG_USE_CLI_GIT := true

# define SERIAL environmental variable to run update serially
ifdef SERIAL
	PARALLELISM := -j 1
else
	PARALLELISM := -j 100
endif

DOCKER_PLATFORM := $(shell uname -m)
DOCKER_PLATFORM := amd64
#DOCKER_PLATFORM := arm64

define clone_update_repo
    if [ ! -d "$(JULIA_PKG_DEVDIR)" ]; then mkdir -p $(JULIA_PKG_DEVDIR); fi
	cd $(JULIA_PKG_DEVDIR);\
	if [ ! -d "$(JULIA_PKG_DEVDIR)/$(1)" ]; then git clone --single-branch git@github.com:ProjectTorreyPines/$(1).jl.git $(1) ; else cd $(1) && git pull origin `git rev-parse --abbrev-ref HEAD` ; fi
endef

all:
	@echo ''
	@echo '  ███████╗██╗   ██╗███████╗███████╗'
	@echo '  ██╔════╝██║   ██║██╔════╝██╔════╝'
	@echo '  █████╗  ██║   ██║███████╗█████╗  '
	@echo '  ██╔══╝  ██║   ██║╚════██║██╔══╝  '
	@echo '  ██║     ╚██████╔╝███████║███████╗'
	@echo '  ╚═╝      ╚═════╝ ╚══════╝╚══════╝'
	@echo ''
	@echo ' - make install      : install FUSE and its dependencies to $(JULIA_PKG_DEVDIR)'
	@echo ' - make update       : git pull FUSE and its dependencies'
	@echo ' - make IJulia       : Install IJulia'
	@echo ' - make dd           : regenerate IMADDD.dd.jl file'
	@echo ' - make html         : generate documentation (FUSE/docs/build/index.html)'
	@echo ''

registry:
	julia -e 'using Pkg;Pkg.add("Revise")' # call this first to make sure General registry gets installed
	if [ ! -d "$(JULIA_PKG_REGDIR)" ]; then mkdir -p $(JULIA_PKG_REGDIR); fi
	cd $(JULIA_PKG_REGDIR);\
	if [ ! -d "$(JULIA_PKG_REGDIR)/GAregistry" ]; then git clone git@github.com:ProjectTorreyPines/GAregistry.git GAregistry ; fi

register:
	$(foreach package,$(PTP_PACKAGES),julia -e 'println("$(package)");using LocalRegistry; register("$(package)", registry="GAregistry")';)

develop:
	julia -e '\
using Pkg;\
Pkg.activate(".");\
Pkg.develop(["IMAS", "IMASDD", "CoordinateConventions", "MillerExtendedHarmonic", "FusionMaterials", "VacuumFields", "Equilibrium", "TAUENN", "EPEDNN", "TGLFNN", "QED", "FiniteElementHermite", "Fortran90Namelists", "CHEASE", "NNeutronics", "SimulationParameters"]);\
Pkg.activate();\
Pkg.develop(["FUSE", "IMAS", "IMASDD", "CoordinateConventions", "MillerExtendedHarmonic", "FusionMaterials", "VacuumFields", "Equilibrium", "TAUENN", "EPEDNN", "TGLFNN", "QED", "FiniteElementHermite", "Fortran90Namelists", "CHEASE", "NNeutronics", "SimulationParameters"]);\
'

rm_manifests:
	find .. -name "Manifest.toml" -exec rm -rf \{\} \;

install_no_registry: clone_update_all develop precompile

install_via_registry: registry develop precompile

install: install_no_registry

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
Pkg.add(["Revise", "JuliaFormatter", "Test", "Plots", "IJulia", "WebIO", "Interact"]);\
Pkg.build("IJulia");\
'
	python3 -m pip install --upgrade webio_jupyter_extension

precompile:
	julia -e 'using Pkg; Pkg.precompile()'

clone_update_all:
	make $(PARALLELISM) FUSE IMAS IMASDD CoordinateConventions MillerExtendedHarmonic FusionMaterials VacuumFields Equilibrium TAUENN EPEDNN TGLFNN QED FiniteElementHermite Fortran90Namelists CHEASE NNeutronics SimulationParameters

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

NNeutronics:
	$(call clone_update_repo,$@)

SimulationParameters:
	$(call clone_update_repo,$@)

docker_clean:
	rm -rf ../Dockerfile
	cp docker/Dockerfile_clean ../Dockerfile
	sed 's/_PLATFORM_/$(DOCKER_PLATFORM)/g' ../Dockerfile > tmp
	mv tmp ../Dockerfile
	cat ../Dockerfile
	cp .gitignore ../.dockerignore
	cd .. ; sudo docker build --platform=linux/$(DOCKER_PLATFORM) -t julia_clean_$(DOCKER_PLATFORM) .

# build a new FUSE docker base image
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

# upload docker image to gke
docker_upload:
	docker tag julia_fuse_amd64_updated gcr.io/sdsc-20220719-60951/fuse
	docker push gcr.io/sdsc-20220719-60951/fuse

cleanup:
	julia -e 'using Pkg; using Dates; Pkg.gc(; collect_delay=Dates.Day(0))'

html:
	cd docs; julia make.jl

web_push:
	cd docs/pages; git reset --hard 049da2c703ad7fc552c13bfe0651da677e3c7f58
	cd docs; cp -rf build/* pages/
	cd docs/pages; echo "fuse.help" > CNAME ### this is to set the custom domain name for gh-pages
	cd docs/pages; touch .nojekyll
	cd docs/pages; git add -A; git commit --allow-empty -m "documentation"; git push --force

web:
	if [ ! -d "$(PWD)/docs/pages" ]; then cd docs; git clone --single-branch -b gh-pages git@github.com:ProjectTorreyPines/FUSE.jl.git pages; fi
	cd docs/pages; git pull
	make web_push

web_ci:
	git clone $(PWD) docs/pages
	cp .git/config docs/pages/.git/config
	cd docs/pages; git fetch; git checkout gh-pages
	cd docs/pages; git config user.email "fuse@fusion.gat.com"
	cd docs/pages; git config user.name "FUSE-BOT"
	make web_push

clean_examples:
	cd docs/src; rm -rf example_*.md

all_examples: clean_examples examples

examples: .PHONY
	cd docs; julia notebooks_to_html.jl --execute

all_blank_examples: clean_examples blank_examples

blank_examples:
	cd docs; julia notebooks_to_html.jl

daily_example:
	cd docs; julia notebooks_to_html.jl --daily --execute --canfail

daily_example_commit:
	git checkout -b examples_$(TODAY)
	git add -A
	git config user.email "fuse-bot@fusion.gat.com"
	git config user.name "fuse bot"
	git config push.autoSetupRemote true
	git commit --allow-empty -m "example of the day"
	git push --set-upstream origin examples_$(TODAY)

dd:
	julia ../IMASDD/src/generate_dd.jl

.PHONY:
