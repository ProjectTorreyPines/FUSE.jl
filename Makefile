JULIA_PKG_DEVDIR ?= $(HOME)/.julia/dev
CURRENTDIR = $(shell pwd)

all:
	@echo 'FUSE makefile help'
	@echo ''
	@echo ' - make install  : install FUSE and all of its dependencies'
	@echo ' - make update   : update FUSE and all of its dependencies'
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

install: install_FUSE install_IJulia
	julia -e '\
using Pkg;\
Pkg.activate();\
Pkg.develop(["FUSE", "IMAS", "IMASDD", "CoordinateConventions", "FusionMaterials", "VacuumFields", "Equilibrium", "TAUENN", "EPEDNN", "TGLFNN", "QED"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_IJulia:
	julia -e '\
using Pkg;\
Pkg.add("IJulia");\
Pkg.build("IJulia");\
'

install_FUSE: install_IMAS install_IMASDD install_FusionMaterials install_VacuumFields install_Equilibrium install_TAUENN install_QED
	if [ ! -d "$(JULIA_PKG_DEVDIR)/FUSE" ]; then ln -s $(CURRENTDIR) $(JULIA_PKG_DEVDIR)/FUSE; fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/FUSE");\
Pkg.develop(["IMAS", "IMASDD", "CoordinateConventions", "FusionMaterials", "VacuumFields", "Equilibrium", "TAUENN", "EPEDNN", "TGLFNN", "QED"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_IMAS: install_IMASDD install_CoordinateConventions
	if [ ! -d "$(JULIA_PKG_DEVDIR)/IMAS" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/IMAS.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/IMAS");\
Pkg.develop(["IMASDD", "CoordinateConventions"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_IMASDD:
	if [ ! -d "$(JULIA_PKG_DEVDIR)/IMASDD" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/IMASDD.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/IMASDD");\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_CoordinateConventions:
	if [ ! -d "$(JULIA_PKG_DEVDIR)/CoordinateConventions" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/CoordinateConventions.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/CoordinateConventions");\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_FusionMaterials:
	if [ ! -d "$(JULIA_PKG_DEVDIR)/FusionMaterials" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/FusionMaterials.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/FusionMaterials");\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_VacuumFields: install_Equilibrium
	if [ ! -d "$(JULIA_PKG_DEVDIR)/VacuumFields" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/VacuumFields.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/VacuumFields");\
Pkg.develop(["Equilibrium", "CoordinateConventions"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_Equilibrium: install_CoordinateConventions
	if [ ! -d "$(JULIA_PKG_DEVDIR)/Equilibrium" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/Equilibrium.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/Equilibrium");\
Pkg.develop(["CoordinateConventions"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_TAUENN: install_IMAS install_IMASDD install_CoordinateConventions install_TGLFNN install_EPEDNN
	if [ ! -d "$(JULIA_PKG_DEVDIR)/TAUENN" ]; then ln -s $(CURRENTDIR) $(JULIA_PKG_DEVDIR)/TAUENN; fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/TAUENN");\
Pkg.develop(["IMAS", "IMASDD", "CoordinateConventions", "TGLFNN", "EPEDNN"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_TGLFNN:
	if [ ! -d "$(JULIA_PKG_DEVDIR)/TGLFNN" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/TGLFNN.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_EPEDNN:
	if [ ! -d "$(JULIA_PKG_DEVDIR)/EPEDNN" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/EPEDNN.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_QED:
	if [ ! -d "$(JULIA_PKG_DEVDIR)/QED" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:ProjectTorreyPines/QED.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

update: update_FUSE update_IMAS update_IMASDD update_CoordinateConventions update_FusionMaterials update_VacuumFields update_Equilibrium update_TAUENN update_EPEDNN update_TGLFNN update_QED
	make install

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
