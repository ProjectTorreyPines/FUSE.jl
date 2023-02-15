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
	julia -e 'using Pkg; Pkg.Registry.update("GAregistry")'

register:
	$(foreach package,$(PTP_PACKAGES),julia -e 'println("$(package)"); using Pkg; Pkg.Registry.update("GAregistry"); Pkg.activate("."); using LocalRegistry; register("$(package)", registry="GAregistry")';)

forward_compatibility:
	julia -e '\
using Pkg;\
for package in ["Equilibrium", "Broker", "ZMQ"];\
	try;\
		Pkg.activate();\
		Pkg.rm(package);\
	catch;\
	end;\
	try;\
		Pkg.activate(".");\
		Pkg.rm(package);\
	catch;\
	end;\
end;\
'

develop:
	# install in global environment to easily develop and test changes made across multiple packages at once
	julia -e '\
fuse_packages = ["IMAS", "IMASDD", "CoordinateConventions", "MillerExtendedHarmonic", "FusionMaterials", "VacuumFields", "MXHEquilibrium", "MeshTools", "TAUENN", "EPEDNN", "TGLFNN", "QED", "FiniteElementHermite", "Fortran90Namelists", "CHEASE", "NNeutronics", "SimulationParameters"];\
using Pkg;\
Pkg.activate(".");\
Pkg.develop(fuse_packages);\
Pkg.activate();\
Pkg.develop(["FUSE"; fuse_packages]);\
Pkg.add(["Revise", "JuliaFormatter", "Test"]);\
'

rm_manifests:
	find .. -name "Manifest.toml" -exec rm -rf \{\} \;

install_no_registry: forward_compatibility clone_update_all develop

install_via_registry: forward_compatibility registry develop

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
Pkg.add(["Plots", "IJulia", "WebIO", "Interact"]);\
Pkg.build("IJulia");\
'
	python3 -m pip install --upgrade webio_jupyter_extension

precompile:
	julia -e 'using Pkg; Pkg.resolve(); Pkg.activate("."); Pkg.instantiate(); Pkg.update(); Pkg.precompile()'

clone_update_all:
	make -i $(PARALLELISM) FUSE IMAS IMASDD CoordinateConventions MillerExtendedHarmonic FusionMaterials VacuumFields MXHEquilibrium MeshTools TAUENN EPEDNN TGLFNN QED FiniteElementHermite Fortran90Namelists CHEASE NNeutronics SimulationParameters

update: install_no_registry precompile

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

MXHEquilibrium:
	$(call clone_update_repo,$@)

MeshTools:
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
	cd docs; julia notebooks_to_md.jl --execute

all_blank_examples: clean_examples blank_examples

blank_examples:
	cd docs; julia notebooks_to_md.jl

daily_example:
	cd docs; julia notebooks_to_md.jl --daily --execute --canfail

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

milestones:
    # NOTE, to get this running: `gh auth login` then `gh extension install valeriobelli/gh-milestone`
	$(foreach package,$(PTP_PACKAGES),cd  ../$(package);\
	gh milestone create --due-date 2023-05-01 --title "Stability: Scaling Laws" --description "BOE 2.02.01.01.04.01: STATIONARY MODEL DEVELOPMENT -> Plasma Core -> Stability -> Scaling laws";\
	gh milestone create --due-date 2025-06-01 --title "Stability: Low-n Ideal MHD" --description "BOE 2.02.01.01.04.02: STATIONARY MODEL DEVELOPMENT -> Plasma Core -> Stability -> Low-n ideal MHD";\
	gh milestone create --due-date 2023-03-01 --title "Fast-Ions: Slowing-Down" --description "BOE 2.02.01.01.05.01: STATIONARY MODEL DEVELOPMENT -> Plasma Core -> Fast-Ions -> Slowing-down";\
	gh milestone create --due-date 2024-04-01 --title "Fast-Ions: Diffusion" --description "BOE 2.02.01.01.05.02: STATIONARY MODEL DEVELOPMENT -> Plasma Core -> Fast-Ions -> Diffusion";\
	gh milestone create --due-date 2023-09-01 --title "Stationary Fueling: Pellet Ablation" --description "BOE 2.02.01.01.07.01: STATIONARY MODEL DEVELOPMENT -> Plasma Core -> Fueling -> Pellet ablation w/ drift";\
	gh milestone create --due-date 2024-01-01 --title "Stationary Fueling: Gas Fueling Rate" --description "BOE 2.02.01.01.07.02: STATIONARY MODEL DEVELOPMENT -> Plasma Core -> Fueling -> Time-averaged gas fueling rate";\
	gh milestone create --due-date 2023-04-01 --title "Scrape-off Layer: Two-Point Model" --description "BOE 2.02.01.01.08.01: STATIONARY MODEL DEVELOPMENT -> Plasma Core -> Scrape-off Layer -> 2-point model";\
	gh milestone create --due-date 2025-04-01 --title "Scrape-off Layer: Neural-Net 1D" --description "BOE 2.02.01.01.08.02: STATIONARY MODEL DEVELOPMENT -> Plasma Core -> Scrape-off Layer -> NN 1D model";\
	gh milestone create --due-date 2023-03-01 --title "Tokamak Assembly: Ports" --description "BOE 2.02.01.02.01.01: STATIONARY MODEL DEVELOPMENT -> Build -> Tokamak Assembly -> Ports";\
	gh milestone create --due-date 2023-05-01 --title "Tokamak Assembly: Maintenance" --description "BOE 2.02.01.02.01.02: STATIONARY MODEL DEVELOPMENT -> Build -> Tokamak Assembly -> Maintenance";\
	gh milestone create --due-date 2024-05-01 --title "Coils: Winding Layouts" --description "BOE 2.02.01.02.03.01: STATIONARY MODEL DEVELOPMENT -> Build -> Coils -> Winding layouts";\
	gh milestone create --due-date 2024-06-01 --title "Coils: Quench Protection" --description "BOE 2.02.01.02.03.02: STATIONARY MODEL DEVELOPMENT -> Build -> Coils -> Quench protection";\
	gh milestone create --due-date 2024-07-01 --title "Coils: Support Strategies" --description "BOE 2.02.01.02.03.03: STATIONARY MODEL DEVELOPMENT -> Build -> Coils -> Support strategies";\
	gh milestone create --due-date 2023-07-01 --title "Blanket: Neutron Damage" --description "BOE 2.02.01.02.04.02: STATIONARY MODEL DEVELOPMENT -> Build -> Blanket -> Neutron damage";\
	gh milestone create --due-date 2025-07-01 --title "Blanket: Materials Comparison" --description "BOE 2.02.01.02.04.03: STATIONARY MODEL DEVELOPMENT -> Build -> Blanket -> Materials comparison";\
	gh milestone create --due-date 2024-04-01 --title "Blanket: Lifetime Estimation" --description "BOE 2.02.01.02.04.04: STATIONARY MODEL DEVELOPMENT -> Build -> Blanket -> Lifetime estimation";\
	gh milestone create --due-date 2023-07-01 --title "Divertor: Geometry" --description "BOE 2.02.01.02.05.02: STATIONARY MODEL DEVELOPMENT -> Build -> Divertor -> Geometry";\
	gh milestone create --due-date 2023-07-01 --title "Divertor: Neutron Damage" --description "BOE 2.02.01.02.05.03: STATIONARY MODEL DEVELOPMENT -> Build -> Divertor -> Neutron damage";\
	gh milestone create --due-date 2023-10-01 --title "Divertor: Heat-Flux Damage" --description "BOE 2.02.01.02.05.04: STATIONARY MODEL DEVELOPMENT -> Build -> Divertor -> Heat-flux damage";\
	gh milestone create --due-date 2024-01-01 --title "Divertor: Lifetime Estimation" --description "BOE 2.02.01.02.05.05: STATIONARY MODEL DEVELOPMENT -> Build -> Divertor -> Lifetime estimation";\
	gh milestone create --due-date 2024-08-01 --title "RF Power: Reduced ECH" --description "BOE 2.02.01.03.01.02: STATIONARY MODEL DEVELOPMENT -> H & CD -> RF Power -> Reduced ECH";\
	gh milestone create --due-date 2025-05-01 --title "RF Power: Beam Tracing" --description "BOE 2.02.01.03.01.03: STATIONARY MODEL DEVELOPMENT -> H & CD -> RF Power -> Beam tracing";\
	gh milestone create --due-date 2024-03-01 --title "Diagnostics: Magnetics" --description "BOE 2.02.01.04.01: STATIONARY MODEL DEVELOPMENT -> Diagnostics -> Magnetics";\
	gh milestone create --due-date 2024-05-01 --title "Diagnostics: Reflectometer" --description "BOE 2.02.01.04.02: STATIONARY MODEL DEVELOPMENT -> Diagnostics -> Reflectometer";\
	gh milestone create --due-date 2024-07-01 --title "Diagnostics: ECE" --description "BOE 2.02.01.04.03: STATIONARY MODEL DEVELOPMENT -> Diagnostics -> ECE";\
	gh milestone create --due-date 2024-04-01 --title "Balance of Plant: Tritium Plant" --description "BOE 2.02.01.05.02: STATIONARY MODEL DEVELOPMENT -> BoP -> Tritium plant";\
	gh milestone create --due-date 2024-04-01 --title "Balance of Plant: Electricity Generation" --description "BOE 2.02.01.05.03: STATIONARY MODEL DEVELOPMENT -> BoP -> Electricity generation";\
	gh milestone create --due-date 2024-01-01 --title "Balance of Plant: H&CD Efficiencies" --description "BOE 2.02.01.05.04: STATIONARY MODEL DEVELOPMENT -> BoP -> H & CD efficiencies";\
	gh milestone create --due-date 2023-05-01 --title "Balance of Plant: Coupling/Transmission Efficiencies" --description "BOE 2.02.01.05.05: STATIONARY MODEL DEVELOPMENT -> BoP -> Coupling and transmission efficiencies";\
	gh milestone create --due-date 2024-03-01 --title "Balance of Plant: Vacuum Pumps" --description "BOE 2.02.01.05.06: STATIONARY MODEL DEVELOPMENT -> BoP -> Vacuum pumps";\
	gh milestone create --due-date 2023-03-01 --title "Balance of Plant: Availability" --description "BOE 2.02.01.05.07: STATIONARY MODEL DEVELOPMENT -> BoP -> Availability";\
	gh milestone create --due-date 2023-07-01 --title "Site: Costing" --description "BOE 2.02.01.06.01: STATIONARY MODEL DEVELOPMENT -> Site -> Costing";\
	gh milestone create --due-date 2024-01-01 --title "Site: Risk" --description "BOE 2.02.01.06.02: STATIONARY MODEL DEVELOPMENT -> Site -> Risk";\
	gh milestone create --due-date 2023-07-01 --title "Site: Infrastructure" --description "BOE 2.02.01.06.03: STATIONARY MODEL DEVELOPMENT -> Site -> Infrastructure";\
	gh milestone create --due-date 2025-04-01 --title "Time-Dependent Breakdown" --description "BOE 2.02.02.01.01: TIME-DEPENDENT MODEL DEVELOPMENT -> Plasma Core -> Breakdown";\
	gh milestone create --due-date 2023-03-01 --title "Time-Dependent Transport: QED" --description "BOE 2.02.02.01.02.01: TIME-DEPENDENT MODEL DEVELOPMENT -> Plasma Core -> Transport -> QED";\
	gh milestone create --due-date 2023-10-01 --title "Time-Dependent Transport: TORTAS" --description "BOE 2.02.02.01.02.02: TIME-DEPENDENT MODEL DEVELOPMENT -> Plasma Core -> Transport -> TORTAS";\
	gh milestone create --due-date 2024-05-01 --title "Time-Dependent Gas Fueling Rate" --description "BOE 2.02.02.01.03.01: TIME-DEPENDENT MODEL DEVELOPMENT -> Plasma Core -> Fueling -> Gas fueling rate";\
	gh milestone create --due-date 2024-03-01 --title "Time-Dependent Pellet Ablation" --description "BOE 2.02.02.01.03.02: TIME-DEPENDENT MODEL DEVELOPMENT -> Plasma Core -> Fueling -> Pellet ablation w/ drift";\
	gh milestone create --due-date 2025-04-01 --title "Time-Dependent ELM" --description "BOE 2.02.02.01.04.01: TIME-DEPENDENT MODEL DEVELOPMENT -> Plasma Core -> Pedestal -> ELM";\
	gh milestone create --due-date 2025-07-01 --title "Time-Dependent L-H Transition" --description "BOE 2.02.02.01.04.02: TIME-DEPENDENT MODEL DEVELOPMENT -> Plasma Core -> Pedestal -> L-H transition";\
	gh milestone create --due-date 2024-03-01 --title "Pulse Design: Workflows" --description "BOE 2.02.02.02.01: TIME-DEPENDENT MODEL DEVELOPMENT -> Pulse Design -> Workflows";\
	gh milestone create --due-date 2024-04-01 --title "Pulse Design: Waveform Definition" --description "BOE 2.02.02.02.02: TIME-DEPENDENT MODEL DEVELOPMENT -> Pulse Design -> Waveform definition";\
	gh milestone create --due-date 2023-12-01 --title "Pulse Design: Shape Evolution" --description "BOE 2.02.02.02.03: TIME-DEPENDENT MODEL DEVELOPMENT -> Pulse Design -> Shape evolution";\
	gh milestone create --due-date 2025-07-01 --title "Time-Dependent Neutronics: Transmutation" --description "BOE 2.02.02.03.01: TIME-DEPENDENT MODEL DEVELOPMENT -> Neutronics -> Transmutation";\
	gh milestone create --due-date 2025-07-01 --title "Time-Dependent Neutronics: Activation/Shutdown Dose Rate" --description "BOE 2.02.02.03.02: TIME-DEPENDENT MODEL DEVELOPMENT -> Neutronics -> Activation/shutdown dose rate";\
	gh milestone create --due-date 2024-07-01 --title "Control: Core Plasma" --description "BOE 2.02.03.01.01: CONTROL INTEGRATION -> Internal Feedback -> Core plasma";\
	gh milestone create --due-date 2024-09-01 --title "Control: Equilibrium" --description "BOE 2.02.03.01.02: CONTROL INTEGRATION -> Internal Feedback -> Equilibrium";\
	gh milestone create --due-date 2024-11-01 --title "Control: Burn" --description "BOE 2.02.03.01.03: CONTROL INTEGRATION -> Internal Feedback -> Burn";\
	gh milestone create --due-date 2025-01-01 --title "Control: Scenario Access" --description "BOE 2.02.03.01.04: CONTROL INTEGRATION -> Internal Feedback -> Scenario access";\
	gh milestone create --due-date 2025-07-01 --title "Control: TokSys Two-Way" --description "BOE 2.02.03.02.01.02: CONTROL INTEGRATION -> External Coupling -> TokSys -> Two-way";\
	gh milestone create --due-date 2024-05-01 --title "Control: GATM One-Way" --description "BOE 2.02.03.02.02.01: CONTROL INTEGRATION -> External Coupling -> GATM -> One-way";\
	gh milestone create --due-date 2026-01-01 --title "Control: GATM Two-Way" --description "BOE 2.02.03.02.02.02: CONTROL INTEGRATION -> External Coupling -> GATM -> Two-way";\
	gh milestone create --due-date 2023-04-01 --title "Support: Cloud Execution" --description "BOE 2.02.04.01.01: SUPPORT, TESTING, AND UPGRADES -> Framework -> Cloud execution";\
	gh milestone create --due-date 2023-04-01 --title "Support: UQ Workflows" --description "BOE 2.02.04.01.02: SUPPORT, TESTING, AND UPGRADES -> Framework -> UQ workflows";)

.PHONY:
