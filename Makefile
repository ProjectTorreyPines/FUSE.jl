all: header branch

help: header
	@echo ' - make install      : install FUSE and its dependencies to $(JULIA_PKG_DEVDIR)'
	@echo ' - make update       : git pull FUSE and its TorreyPines dependencies'
	@echo ' - make update_all   : git pull FUSE and all of its dependencies'
	@echo ' - make IJulia       : Install IJulia'
	@echo ' - make dd           : regenerate IMADDD.dd.jl file'
	@echo ' - make html         : generate documentation (FUSE/docs/build/index.html)'
	@echo ''

header:
	@echo ''
	@echo '  ███████╗██╗   ██╗███████╗███████╗'
	@echo '  ██╔════╝██║   ██║██╔════╝██╔════╝'
	@echo '  █████╗  ██║   ██║███████╗█████╗  '
	@echo '  ██╔══╝  ██║   ██║╚════██║██╔══╝  '
	@echo '  ██║     ╚██████╔╝███████║███████╗'
	@echo '  ╚═╝      ╚═════╝ ╚══════╝╚══════╝'
	@echo ''

# =========================

realpath = $(shell cd $(dir $(1)); pwd)/$(notdir $(1))
JULIA_DIR ?= $(call realpath,$(HOME)/.julia)
JULIA_CONF := $(JULIA_DIR)/config/startup.jl
JULIA_PKG_REGDIR ?= $(JULIA_DIR)/registries
JULIA_PKG_DEVDIR ?= $(JULIA_DIR)/dev
CURRENTDIR := $(shell (pwd -P))
TODAY := $(shell date +'%Y-%m-%d')
export JULIA_NUM_THREADS ?= $(shell julia -e "println(length(Sys.cpu_info()))")

ifdef GITHUB_HEAD_REF
  # For pull_request events
  FUSE_LOCAL_BRANCH=$(GITHUB_HEAD_REF)
else
  # For push events and others
  FUSE_LOCAL_BRANCH=$(shell echo $(GITHUB_REF) | sed 's/refs\/heads\///')
endif

FUSE_PACKAGES_MAKEFILE := ADAS BoundaryPlasmaModels CHEASE CoordinateConventions EPEDNN FiniteElementHermite Fortran90Namelists FuseUtils FusionMaterials FXP IMAS IMASDD MXHEquilibrium MeshTools MillerExtendedHarmonic NEO NNeutronics QED SimulationParameters TEQUILA TGLFNN TJLF VacuumFields XSteamTP ThermalSystemModels
FUSE_PACKAGES_MAKEFILE := $(sort $(FUSE_PACKAGES_MAKEFILE))
FUSE_PACKAGES := $(shell echo '$(FUSE_PACKAGES_MAKEFILE)' | awk '{printf("[\"%s\"", $$1); for (i=2; i<=NF; i++) printf(", \"%s\"", $$i); print "]"}')
DEV_PACKAGES := $(shell find ../*/.git/config -exec grep ProjectTorreyPines \{\} \; | cut -d'/' -f 2 | cut -d'.' -f 1 | tr '\n' ' ')

# use command line interface for git to work nicely with private repos
export JULIA_PKG_USE_CLI_GIT := true

# define SERIAL environmental variable to run update serially
ifdef SERIAL
	PARALLELISM := -j 1
else
	PARALLELISM := -j 100
endif

define clone_pull_repo
	@ if [ ! -d "$(JULIA_PKG_DEVDIR)" ]; then mkdir -p $(JULIA_PKG_DEVDIR); fi
	@ cd $(JULIA_PKG_DEVDIR); if [ ! -d "$(JULIA_PKG_DEVDIR)/$(1)" ]; then git clone git@github.com:ProjectTorreyPines/$(1).jl.git $(1) ; else cd $(1) && git pull origin `git rev-parse --abbrev-ref HEAD` ; fi
endef

define feature_or_master_julia
function feature_or_master(package, feature_branch) ;\
	token = "$(PTP_READ_TOKEN)" ;\
	url = "https://api.github.com/repos/ProjectTorreyPines/$$(package).jl/branches/$$(feature_branch)" ;\
	;\
	curl_cmd = `curl -s -o /dev/null -w "%{http_code}" -L -H "Authorization: Bearer $$token" -H "Accept: application/vnd.github+json" -H "X-GitHub-Api-Version: 2022-11-28" $$url` ;\
	;\
	http_status = chomp(read(`$$curl_cmd`, String)) ;\
	;\
	if http_status == "200" ;\
		return feature_branch ;\
	elseif http_status == "404" ;\
		return "master" ;\
	else ;\
		error("GitHub API returned status code: $$http_status") ;\
	end ;\
end
endef

# =========================

# simple test to see how many threads julia will run on (set by JULIA_NUM_THREADS)
threads:
	@echo "set the JULIA_NUM_THREADS environment variable"
	julia -e "println(Threads.nthreads())"

# remove everything under $HOME/.julia besides $HOME/.julia/dev
nuke_julia:
	mv $(JULIA_PKG_DEVDIR) $(call realpath,$(JULIA_DIR))/../asddsaasddsa
	rm -rf $(JULIA_DIR)
	mkdir -p $(JULIA_DIR)
	mv $(call realpath,$(JULIA_DIR))/../asddsaasddsa $(JULIA_PKG_DEVDIR)

# install the FuseRegistry to the list of Julia registries
registry:
	julia -e 'using Pkg; Pkg.add("Revise"); pkg"registry add git@github.com:ProjectTorreyPines/FuseRegistry.jl.git"; Pkg.add("LocalRegistry")'

# register a package to FuseRegistry
# >> make register repo=IMASDD
register:
	julia -e '\
using Pkg;\
Pkg.Registry.update("FuseRegistry");\
Pkg.activate();\
using LocalRegistry;\
LocalRegistry.is_dirty(path, gitconfig)= false; register("$(repo)", registry="FuseRegistry")'

# register all packages with FuseRegistry
all_register:
	$(foreach package,FUSE $(FUSE_PACKAGES_MAKEFILE), $(MAKE) register repo=$(package);)

# install FUSE packages in global environment to easily develop and test changes made across multiple packages at once 
develop:
	julia -e '\
fuse_packages = $(FUSE_PACKAGES);\
println(fuse_packages);\
using Pkg;\
Pkg.activate();\
Pkg.develop([["FUSE"] ; fuse_packages]);\
Pkg.add(["JuliaFormatter", "Test", "Plots"]);\
Pkg.activate(".");\
Pkg.develop(fuse_packages);\
'
	make revise

# load Revise when Julia starts up
revise:
	julia -e 'using Pkg; Pkg.add("Revise")'
	mkdir -p $(JULIA_DIR)/config
	touch $(JULIA_CONF)
	grep -v -F -x "using Revise" "$(JULIA_CONF)" > "$(JULIA_CONF).tmp" || true
	echo "using Revise" | cat - "$(JULIA_CONF).tmp" > "$(JULIA_CONF)"
	rm -f "$(JULIA_CONF).tmp"

# list branches of all the ProjectTorreyPines packages used by FUSE
branch: .PHONY
	@cd $(CURRENTDIR); $(foreach package,FUSE ServeFUSE $(FUSE_PACKAGES_MAKEFILE),printf "%25s" "$(package)"; echo ":  `cd ../$(package); git rev-parse --abbrev-ref HEAD | sed 's/$$/ \*/' | sed 's/^master \*$$/master/'`";)

# install (add) FUSE via HTTPS and $PTP_READ_TOKEN
# looks for same branch name for all repositories otherwise falls back to master
https_add:
	julia -e ';\
$(feature_or_master_julia);\
fuse_packages = $(FUSE_PACKAGES);\
println(fuse_packages);\
using Pkg;\
Pkg.activate(".");\
dependencies = Pkg.PackageSpec[];\
for package in fuse_packages;\
	branch = feature_or_master(package, "$(FUSE_LOCAL_BRANCH)");\
	println("$$(package) $$(branch)");\
	push!(dependencies, Pkg.PackageSpec(url="https://project-torrey-pines:$(PTP_READ_TOKEN)@github.com/ProjectTorreyPines/"*package*".jl.git", rev=branch));\
end;\
Pkg.add(dependencies)'

# install (dev) FUSE via HTTPS and $PTP_READ_TOKEN (needed for documentation)
# all repos are on the master branch (this should be the case for generating documentation)
https_dev:
	@mkdir -p ~/.julia/dev
	@ln -sf $(PWD) ~/.julia/dev/FUSE
	@julia -e ';\
fuse_packages = $(FUSE_PACKAGES);\
using Pkg;\
Pkg.activate(".");\
dependencies = Pkg.PackageSpec[];\
for package in fuse_packages;\
	push!(dependencies, Pkg.PackageSpec(url="https://project-torrey-pines:$(PTP_READ_TOKEN)@github.com/ProjectTorreyPines/"*package*".jl.git"));\
end;\
Pkg.develop(dependencies);\
Pkg.develop(fuse_packages);\
Pkg.activate("./docs");\
Pkg.develop(["FUSE"; fuse_packages])'

# install FUSE without using the registry
install_no_registry: forward_compatibility clone_pull_all develop special_dependencies

# install FUSE using the registry (requires registry to be up-to-date, which most likely are not! Don't use!)
install_via_registry: forward_compatibility registry develop special_dependencies

# install used by CI (add packages, do not dev them)
install_ci_add: https_add special_dependencies
install_ci_dev: https_dev special_dependencies

# set default install method
install: install_no_registry

# dependencies that are not in the global registry
special_dependencies: #DISABLED
#	@julia -e ';\
using Pkg;\
dependencies = Pkg.PackageSpec[PackageSpec(url="https://github.com/IanButterworth/Flux.jl", rev="ib/cudaext")];\
Pkg.add(dependencies);\
'

# update_all, a shorthand for install and precompile
update_all: install
	julia -e 'using Pkg; Pkg.resolve(); Pkg.activate("."); Pkg.resolve(); Pkg.update(); Pkg.precompile()'

# update, a synonim of clone_pull and develop
update: forward_compatibility clone_pull_all develop resolve

# resolve the current environment (eg. after manually adding a new package)
resolve:
	julia -e 'using Pkg; Pkg.resolve(); Pkg.activate("."); Pkg.resolve(); Pkg.precompile()'

# delete local packages that have become obsolete
forward_compatibility:
	julia -e '\
using Pkg;\
for package in ("Equilibrium", "Broker", "ZMQ", "FUSE_GA");\
	try; Pkg.activate();    Pkg.rm(package); catch; end;\
	try; Pkg.activate("."); Pkg.rm(package); catch; end;\
	try; Pkg.activate("./docs"); Pkg.rm(package); catch; end;\
end;\
'
# undo --single-branch clones of git repos
undo_single_branch:
	$(foreach package,$(FUSE_PACKAGES_MAKEFILE),cd ../$(package)/; echo `pwd`; git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"; git fetch origin;)

# clone and update all FUSE packages
clone_pull_all: branch
	@ if [ ! -d "$(JULIA_PKG_DEVDIR)" ]; then mkdir -p $(JULIA_PKG_DEVDIR); fi
	make -i $(PARALLELISM) FUSE ServeFUSE $(FUSE_PACKAGES_MAKEFILE)

playground: .PHONY
	if [ -d playground ] && [ ! -f playground/.gitattributes ]; then mv playground playground_private ; fi
	if [ ! -d "playground" ]; then git clone git@github.com:ProjectTorreyPines/FUSE_playground.git playground ; else cd playground && git pull origin `git rev-parse --abbrev-ref HEAD` ; fi

ADAS:
	$(call clone_pull_repo,$@)

FUSE:
	$(call clone_pull_repo,$@)

IMAS:
	$(call clone_pull_repo,$@)

IMASDD:
	$(call clone_pull_repo,$@)

CoordinateConventions:
	$(call clone_pull_repo,$@)

MillerExtendedHarmonic:
	$(call clone_pull_repo,$@)

FuseUtils:
	$(call clone_pull_repo,$@)

FusionMaterials:
	$(call clone_pull_repo,$@)

FXP:
	$(call clone_pull_repo,$@)

VacuumFields:
	$(call clone_pull_repo,$@)

MXHEquilibrium:
	$(call clone_pull_repo,$@)

MeshTools:
	$(call clone_pull_repo,$@)

TGLFNN:
	$(call clone_pull_repo,$@)

TJLF:
	$(call clone_pull_repo,$@)

EPEDNN:
	$(call clone_pull_repo,$@)

QED:
	$(call clone_pull_repo,$@)

FiniteElementHermite:
	$(call clone_pull_repo,$@)

Fortran90Namelists:
	$(call clone_pull_repo,$@)

CHEASE:
	$(call clone_pull_repo,$@)

TEQUILA:
	$(call clone_pull_repo,$@)

NNeutronics:
	$(call clone_pull_repo,$@)

SimulationParameters:
	$(call clone_pull_repo,$@)

BoundaryPlasmaModels:
	$(call clone_pull_repo,$@)

NEO:
	$(call clone_pull_repo,$@)

XSteamTP:
	$(call clone_pull_repo,$@)

ThermalSystemModels:
	@find ~/.julia/environments/ -type f -name "*.toml" -exec sed -i.bak 's/ThermalSystem_Models/ThermalSystemModels/g' {} \; -exec rm -f {}.bak \;
	@find ../ -type f -name "*.toml" -exec sed -i.bak 's/ThermalSystem_Models/ThermalSystemModels/g' {} \; -exec rm -f {}.bak \;
	@rm -rf ../ThermalSystem_Models
	$(call clone_pull_repo,$@)

ServeFUSE:
	$(call clone_pull_repo,$@)
	julia -e '\
fuse_packages = $(FUSE_PACKAGES);\
println(fuse_packages);\
using Pkg;\
Pkg.activate("../ServeFUSE/task_tiller");\
Pkg.develop([["FUSE"] ; fuse_packages]);\
'

# Install IJulia
IJulia:
	julia -e '\
using Pkg;\
Pkg.add(["Plots", "IJulia", "WebIO", "Interact"]);\
Pkg.build("IJulia");\
import IJulia;\
n=get(ENV, "JULIA_NUM_THREADS", string(length(Sys.cpu_info())));\
IJulia.installkernel("Julia ("*n*" threads)"; env=Dict("JULIA_NUM_THREADS"=>n));\
'
	jupyter kernelspec list
	python3 -m pip install --upgrade webio_jupyter_extension
	jupyter nbextension list
	jupyter labextension list

# Install PyCall
PyCall:
	julia -e '\
ENV["PYTHON"]="";\
using Pkg;\
Pkg.add("PyCall");\
Pkg.build("PyCall");\
'

# setup ./docs environment to build documentation
develop_docs:
	julia -e '\
fuse_packages = $(FUSE_PACKAGES);\
println(fuse_packages);\
using Pkg;\
Pkg.activate("./docs");\
Pkg.develop([["FUSE"] ; fuse_packages]);\
Pkg.add(["Documenter", "Pkg", "ProgressMeter", "InteractiveUtils"]);\
'

# generate documentation
html: develop_docs
	cd docs; julia make.jl

# push documentation to the web
web_push:
	cd docs/pages; git reset --hard 049da2c703ad7fc552c13bfe0651da677e3c7f58
	cd docs; cp -rf build/* pages/
	cd docs/pages; echo "fuse.help" > CNAME ### this is to set the custom domain name for gh-pages
	cd docs/pages; touch .nojekyll
	cd docs/pages; git add -A; git commit --allow-empty -m "documentation"; git push --force

# push documentation to the web (user entry point)
web:
	if [ ! -d "$(PWD)/docs/pages" ]; then cd docs; git clone --single-branch -b gh-pages git@github.com:ProjectTorreyPines/FUSE.jl.git pages; fi
	make web_push

# push documentation to the web (CI entry point)
web_ci:
	git clone $(PWD) docs/pages
	cp .git/config docs/pages/.git/config
	cd docs/pages; git fetch; git checkout gh-pages
	cd docs/pages; git config user.email "fuse@fusion.gat.com"
	cd docs/pages; git config user.name "FUSE-BOT"
	make web_push

# clean all examples
clean_examples:
	cd docs/src; rm -rf example_*.md

# clean and run all examples
all_examples: clean_examples examples

# run all examples
examples: .PHONY
	cd docs; julia notebooks_to_md.jl --execute

# clean and convert examples to md without executing
all_blank_examples: clean_examples blank_examples

# convert examples to md without executing
blank_examples:
	cd docs; julia notebooks_to_md.jl

# run daily example to md
daily_example:
	cd docs; julia notebooks_to_md.jl --daily --execute --canfail

# commit daily example md (this must only be run by the CI)
ifdef GITHUB_ACTION
daily_example_ci_commit:
	git config user.email "fuse-bot@fusion.gat.com"
	git config user.name "fuse bot"
	git config push.autoSetupRemote true
	git checkout -b examples_$(TODAY)
	git add -A
	git commit --allow-empty -m "example of the day"
	git push --set-upstream origin examples_$(TODAY)
endif

# commit manifest (this must only be run by the CI)
ifdef GITHUB_ACTION
manifest_ci_commit:
	git config user.email "fuse-bot@fusion.gat.com"
	git config user.name "fuse bot"
	git config push.autoSetupRemote true
	git fetch
	git checkout manifest
	git merge master
	@sed 's/https:\/\/project-torrey-pines:$(PTP_READ_TOKEN)@/https:\\/\\//g' Manifest.toml > Manifest_CI.toml
	git add Manifest_CI.toml
	git commit --allow-empty -m "Manifest $(TODAY)"
	git push --set-upstream origin manifest
endif

# run julia using the Manifest_CI.toml
manifest_ci:
	@TEMP_DIR=$$(mktemp -d /var/tmp/manifest_ci.XXXXXX) && echo $$TEMP_DIR &&\
	sed "s/git@/https:\\/\\//g" Manifest_CI.toml > $$TEMP_DIR/Manifest.toml && \
	cd $$TEMP_DIR && julia -i -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'

# remove all Manifest.toml files
rm_manifests:
	find ..  -maxdepth 3 -type f -name "Manifest.toml" -exec rm -f \{\} \;

# update dd from the json files
dd:
	julia generate_dd/generate_dd.jl

# generates init_expressions.json file, which lists entries that are
# always expected to be expressions when coming out of init()
init_expressions:
	julia -e 'import FUSE; FUSE.init_expressions(;save=true)'

# create an empty commit
empty_commit:
	git reset HEAD
	git commit --allow-empty -m 'empty commit'

# GitHub merge of `branch` into `master` for a series of repos
# >> make branch_master branch=my_branch repos='FUSE IMAS IMASDD'
branch_master:
ifeq ($(repos),)
	$(error repos variable is not set)
endif
	$(foreach repo,$(repos), \
curl -X POST \
-H "Authorization: token $$(security find-generic-password -a orso82 -s GITHUB_TOKEN -w)" \
-H "Accept: application/vnd.github.v3+json" \
https://api.github.com/repos/ProjectTorreyPines/$(repo).jl/merges \
-d '{"base": "master", "head": "$(branch)", "commit_message": "merging $(branch) into master"}';)

# update LICENSE, NOTICE.md, github workflows, docs, juliaformatter and gitignore in preparation of public release
# The starting information is taken from IMASDD.jl and moved to the target repo
# >> make apache repo=CHEASE
# in addition, one must add the DOCUMENTER_KEY to the repo
# https://m3g.github.io/JuliaNotes.jl/stable/publish_docs/#How-to-deploy-the-documentation-of-a-project
apache:
ifeq ($(repo),)
	$(error repo variable is not set)
endif
	echo $(repo)
	cp ../IMASDD/LICENSE ../$(repo)/ ;\
\
cp ../IMASDD/NOTICE.md ../$(repo)/ ;\
sed -i.bak "s/IMASDD/$(repo)/g" ../$(repo)/NOTICE.md && rm ../$(repo)/NOTICE.md.bak ;\
\
mkdir -p ../$(repo)/.github/workflows ;\
cp ../IMASDD/.github/workflows/make_docs.yml ../$(repo)/.github/workflows/ ;\
cp ../IMASDD/.github/workflows/CompatHelper.yml ../$(repo)/.github/workflows/ ;\
cp ../IMASDD/.github/workflows/TagBot.yml ../$(repo)/.github/workflows/ ;\
\
cp ../IMASDD/README.md ../$(repo)/ ;\
sed -i.bak "s/IMASDD/$(repo)/g" ../$(repo)/README.md && rm ../$(repo)/README.md.bak ;\
\
cp -R ../IMASDD/docs ../$(repo)/ ;\
sed -i.bak "s/IMASDD/$(repo)/g" ../$(repo)/docs/make.jl && rm ../$(repo)/docs/make.jl.bak ;\
echo "# $(repo).jl" > ../$(repo)/docs/src/index.md ;\
rm ../$(repo)/docs/Manifest.toml ;\
rm -rf ../$(repo)/docs/build ;\
\
cp -R ../IMASDD/.JuliaFormatter.toml ../$(repo)/ ;\
\
cp -R ../IMASDD/.gitignore ../$(repo)/ ;\
\
julia -e 'import Pkg; Pkg.add("DocumenterTools"); import DocumenterTools; DocumenterTools.genkeys()'

# loop over all FUSE packages
all_apache:
	$(foreach package,FUSE $(FUSE_PACKAGES_MAKEFILE), $(MAKE) apache repo=$(package);)

# bump patch version of a repo
# >> make bump_version repo=IMAS
bump_version:
ifeq ($(repo),)
	$(error repo variable is not set)
endif
	echo $(repo) ;\
new_version=$$(awk '/^version =/ {split($$3, a, "\""); split(a[2], v, "."); v[3]++; printf "%d.%d.%d", v[1], v[2], v[3]}' ../$(repo)/Project.toml) ;\
awk '/^version =/ {split($$3, a, "\""); split(a[2], v, "."); v[3]++; printf "version = \"%d.%d.%d\"\n", v[1], v[2], v[3]; next} {print}' ../$(repo)/Project.toml > ../$(repo)/Project.tmp && mv ../$(repo)/Project.tmp ../$(repo)/Project.toml ;\
git -C ../$(repo) add Project.toml ;\
git -C ../$(repo) commit -m "v$${new_version}" ;\
git -C ../$(repo) tag -d v$${new_version} ;\
git -C ../$(repo) tag -a v$${new_version} -m "v$${new_version}" ;\
git -C ../$(repo) push

# loop over FUSE packages and bump up their versions
# >> make bump_versions repos="IMAS IMASDD"
bump_versions:
ifeq ($(repos),)
	$(error repos variable is not set)
endif
	$(foreach package,$(repos), $(MAKE) bump_version repo=$(package);)

# loop over all FUSE packages and bump up their versions
all_bump_version:
	$(foreach package,FUSE $(FUSE_PACKAGES_MAKEFILE), $(MAKE) bump_version repo=$(package);)


.PHONY:
