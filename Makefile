JULIA_PKG_DEVDIR ?= $(HOME)/.julia/dev
CURRENTDIR = $(shell pwd)

all:
	@echo 'FUSE makefile help'
	@echo ''
	@echo ' - make install  : install FUSE and all of its dependencies
	@echo ' - make update   : update FUSE and all of its dependencies
	@echo ''

install: install_FUSE
	julia -e '\
using Pkg;\
Pkg.activate();\
Pkg.develop(["FUSE", "Equilibrium", "IMAS", "AD_GS"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_FUSE: install_Equilibrium install_IMAS install_AD_GS
	if [ ! -d "$(JULIA_PKG_DEVDIR)/FUSE" ]; then ln -s $(CURRENTDIR) $(JULIA_PKG_DEVDIR)/FUSE; fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/FUSE");\
Pkg.develop(["Equilibrium", "IMAS", "AD_GS"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_IMAS: install_Equilibrium
	if [ ! -d "$(JULIA_PKG_DEVDIR)/IMAS" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:JuliaFusion/IMAS.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/IMAS");\
Pkg.develop(["Equilibrium"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_AD_GS: install_Equilibrium
	if [ ! -d "$(JULIA_PKG_DEVDIR)/AD_GS" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:JuliaFusion/AD_GS.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/AD_GS");\
Pkg.develop(["Equilibrium"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

install_Equilibrium:
	if [ ! -d "$(JULIA_PKG_DEVDIR)/EFIT" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:JuliaFusion/EFIT.jl.git");';\
	fi
	if [ ! -d "$(JULIA_PKG_DEVDIR)/Equilibrium" ]; then\
		julia -e 'using Pkg; Pkg.develop(url="git@github.com:JuliaFusion/Equilibrium.jl.git");';\
	fi
	julia -e '\
using Pkg;\
Pkg.activate("$(JULIA_PKG_DEVDIR)/Equilibrium");\
Pkg.develop(["EFIT"]);\
Pkg.resolve();\
try Pkg.upgrade_manifest() catch end;\
'

update: update_Equilibrium update_IMAS update_AD_GS install

update_FUSE:
	cd $(JULIA_PKG_DEVDIR)/FUSE; git fetch; git pull

update_IMAS:
	cd $(JULIA_PKG_DEVDIR)/IMAS; git fetch; git pull

update_AD_GS:
	cd $(JULIA_PKG_DEVDIR)/AD_GS; git fetch; git pull

update_Equilibrium:
	cd $(JULIA_PKG_DEVDIR)/Equilibrium; git fetch; git pull

.PHONY:
