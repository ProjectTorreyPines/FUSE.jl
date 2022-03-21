JULIA_PKG_REGDIR ?= $(HOME)/.julia/registries
JULIA_PKG_DEVDIR ?= $(HOME)/.julia/dev
CURRENTDIR = $(shell pwd)

all:
	@echo 'FUSE makefile help'
	@echo ''
	@echo ' - make develop  : install FUSE and its PTP dependencies to $(JULIA_PKG_DEVDIR)'
	@echo ' - make update   : git pull FUSE and its PTP dependencies'
	@echo ''

sysimage:
	julia -e '\
using Pkg;\
Pkg.add("PackageCompiler");\
Pkg.add("IJulia");\
import PackageCompiler;\
Pkg.activate(".");\
PackageCompiler.create_sysimage(["Contour", "DataStructures", "EFIT", "ForwardDiff", "Interpolations", "JSON", "LibGEOS", "LinearAlgebra", "ModelingToolkit", "NumericalIntegration", "Optim", "OrdinaryDiffEq", "Plots", "PolygonOps", "Printf", "Random", "Revise", "StaticArrays", "Statistics", "Test"], sysimage_path="FUSEsysimage.so");\
import IJulia;\
IJulia.installkernel("Julia FUSEsysimage", "--sysimage=$(shell pwd)/FUSEsysimage.so", "--trace-compile=stderr");\
'

registry:
	cd $(JULIA_PKG_REGDIR);\
	if [ ! -d "$(JULIA_PKG_REGDIR)/GAregistry" ]; then git clone git@github.com:ProjectTorreyPines/GAregistry.git GAregistry ; fi

develop: registry
	julia -e '\
using Pkg;\
Pkg.activate(".");\
Pkg.develop(["IMAS", "IMASDD", "CoordinateConventions", "FusionMaterials", "VacuumFields", "Equilibrium", "TAUENN", "EPEDNN", "TGLFNN", "QED"]);\
'

IJulia:
	julia -e '\
using Pkg;\
Pkg.add("IJulia");\
Pkg.build("IJulia");\
'

update: develop
	make -j 100 update_FUSE update_IMAS update_IMASDD update_CoordinateConventions update_FusionMaterials update_VacuumFields update_Equilibrium update_TAUENN update_EPEDNN update_TGLFNN update_QED

update_FUSE:
	cd $(JULIA_PKG_DEVDIR)/FUSE; git fetch; git pull

update_IMAS:
	cd $(JULIA_PKG_DEVDIR)/IMAS; git fetch; git pull

update_IMASDD:
	cd $(JULIA_PKG_DEVDIR)/IMASDD; git fetch; git pull

update_CoordinateConventions:
	cd $(JULIA_PKG_DEVDIR)/CoordinateConventions; git fetch; git pull

update_FusionMaterials:
	cd $(JULIA_PKG_DEVDIR)/FusionMaterials; git fetch; git pull

update_VacuumFields:
	cd $(JULIA_PKG_DEVDIR)/VacuumFields; git fetch; git pull

update_Equilibrium:
	cd $(JULIA_PKG_DEVDIR)/Equilibrium; git fetch; git pull

update_TAUENN:
	cd $(JULIA_PKG_DEVDIR)/TAUENN; git fetch; git pull

update_TGLFNN:
	cd $(JULIA_PKG_DEVDIR)/TGLFNN; git fetch; git pull

update_EPEDNN:
	cd $(JULIA_PKG_DEVDIR)/EPEDNN; git fetch; git pull

update_QED:
	cd $(JULIA_PKG_DEVDIR)/QED; git fetch; git pull

.PHONY:
