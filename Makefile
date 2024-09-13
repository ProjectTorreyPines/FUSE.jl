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

FUSE_PACKAGES_MAKEFILE := ADAS BoundaryPlasmaModels CHEASE CoordinateConventions EPEDNN FiniteElementHermite Fortran90Namelists FRESCO FuseUtils FusionMaterials FuseExchangeProtocol IMAS IMASdd MXHEquilibrium MeshTools MillerExtendedHarmonic NEO NNeutronics QED RABBIT SimulationParameters TEQUILA TGLFNN TJLF VacuumFields XSteam ThermalSystemModels
FUSE_PACKAGES_MAKEFILE := $(sort $(FUSE_PACKAGES_MAKEFILE))
FUSE_PACKAGES := $(shell echo '$(FUSE_PACKAGES_MAKEFILE)' | awk '{printf("[\"%s\"", $$1); for (i=2; i<=NF; i++) printf(", \"%s\"", $$i); print "]"}')
DEV_PACKAGES := $(shell find ../*/.git/config -exec grep ProjectTorreyPines \{\} /dev/null \; | cut -d'/' -f 2)

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
	@ cd $(JULIA_PKG_DEVDIR); if [ ! -d "$(JULIA_PKG_DEVDIR)/$(1)" ]; then git clone git@github.com:ProjectTorreyPines/$(1).jl.git $(1) ; else cd $(1) && git pull 2>&1 | sed 's/^/$(1): /'; fi
endef

define feature_or_master_julia
try ;\
	using HTTP ;\
catch ;\
	using Pkg ;\
	Pkg.add("HTTP") ;\
	using HTTP ;\
end ;\
;\
function feature_or_master(package, feature_branch) ;\
    token = "$(PTP_READ_TOKEN)" ;\
    url = "https://api.github.com/repos/ProjectTorreyPines/$$(package).jl/branches/$$(feature_branch)" ;\
    headers = ["Authorization" => "Bearer $$(token)", "Accept" => "application/vnd.github+json", "X-GitHub-Api-Version" => "2022-11-28"] ;\
    response = HTTP.get(url, headers; status_exception=false) ;\
    if response.status == 200 ;\
        return feature_branch ;\
    elseif response.status == 404 ;\
        return "master" ;\
    else ;\
        error("GitHub API returned status code: $$(response.status)") ;\
    end ;\
end
endef

ADAS:
	$(call clone_pull_repo,$@)

FUSE:
	$(call clone_pull_repo,$@)

IMAS:
	$(call clone_pull_repo,$@)

IMASdd:
	$(call clone_pull_repo,$@)

CoordinateConventions:
	$(call clone_pull_repo,$@)

MillerExtendedHarmonic:
	$(call clone_pull_repo,$@)

FuseUtils:
	$(call clone_pull_repo,$@)

FusionMaterials:
	$(call clone_pull_repo,$@)

FuseExchangeProtocol:
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

FRESCO:
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

XSteam:
	$(call clone_pull_repo,$@)

RABBIT:
	$(call clone_pull_repo,$@)

ThermalSystemModels:
	$(call clone_pull_repo,$@)

GenerateDD:
	$(call clone_pull_repo,$@)

ServeFUSE:
	$(call clone_pull_repo,$@)
	@julia -e '\
	fuse_packages = $(FUSE_PACKAGES);\
	using Pkg;\
	Pkg.activate("../ServeFUSE/task_tiller");\
	Pkg.develop([["FUSE"] ; fuse_packages]);\
	'

help: header help_info

.PHONY:

# =%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%= #

blank_examples:clean_examples
# Clean and convert examples to md without executing
	cd docs; julia notebooks_to_md.jl

# @devs
branch_master: error_missing_repos_var
# GitHub merge of `branch` into `master` for a series of repos
# >> make branch_master branch=my_branch repos='FUSE IMAS IMASdd'
	@$(foreach repo,$(repos), \
	curl -X POST \
	-H "Authorization: token $$(security find-generic-password -a orso82 -s GITHUB_TOKEN -w)" \
	-H "Accept: application/vnd.github.v3+json" \
	https://api.github.com/repos/ProjectTorreyPines/$(repo).jl/merges \
	-d '{"base": "master", "head": "$(branch)", "commit_message": "merging $(branch) into master"}';)

# @devs
bump_major: error_missing_repo_var error_on_previous_commit_is_a_version_bump versions_used
# Bump major version of a repo
# >> make bump_major repo=IMAS
	@git -C ../$(repo) checkout -- Project.toml ;\
	new_version=$$(awk '/^version =/ {split($$3, a, "\""); split(a[2], v, "."); v[1]++; v[2]=0; v[3]=0; printf "%d.%d.%d", v[1], v[2], v[3]}' ../$(repo)/Project.toml) ;\
	echo bumping $(repo) to version $${new_version}\\n;\
	awk '/^version =/ {split($$3, a, "\""); split(a[2], v, "."); v[1]++; v[2]=0; v[3]=0; printf "version = \"%d.%d.%d\"\n", v[1], v[2], v[3]; next} {print}' ../$(repo)/Project.toml > ../$(repo)/Project.tmp && mv ../$(repo)/Project.tmp ../$(repo)/Project.toml ;\
	make commit_project_toml repo=$(repo) new_version=$${new_version}

# @devs
bump_minor: error_missing_repo_var error_on_previous_commit_is_a_version_bump versions_used
# Bump minor version of a repo
# >> make bump_minor repo=IMAS
	@git -C ../$(repo) checkout -- Project.toml ;\
	new_version=$$(awk '/^version =/ {split($$3, a, "\""); split(a[2], v, "."); v[2]++; v[3]=0; printf "%d.%d.%d", v[1], v[2], v[3]}' ../$(repo)/Project.toml) ;\
	echo bumping $(repo) to version $${new_version}\\n;\
	awk '/^version =/ {split($$3, a, "\""); split(a[2], v, "."); v[2]++; v[3]=0; printf "version = \"%d.%d.%d\"\n", v[1], v[2], v[3]; next} {print}' ../$(repo)/Project.toml > ../$(repo)/Project.tmp && mv ../$(repo)/Project.tmp ../$(repo)/Project.toml ;\
	make commit_project_toml repo=$(repo) new_version=$${new_version}

# @devs
bump_patch: error_missing_repo_var error_on_previous_commit_is_a_version_bump versions_used
# Bump patch version of a repo
# >> make bump_patch repo=IMAS
	@git -C ../$(repo) checkout -- Project.toml ;\
	new_version=$$(awk '/^version =/ {split($$3, a, "\""); split(a[2], v, "."); v[3]++; printf "%d.%d.%d", v[1], v[2], v[3]}' ../$(repo)/Project.toml) ;\
	echo bumping $(repo) to version $${new_version}\\n;\
	awk '/^version =/ {split($$3, a, "\""); split(a[2], v, "."); v[3]++; printf "version = \"%d.%d.%d\"\n", v[1], v[2], v[3]; next} {print}' ../$(repo)/Project.toml > ../$(repo)/Project.tmp && mv ../$(repo)/Project.tmp ../$(repo)/Project.toml ;\
	make commit_project_toml repo=$(repo) new_version=$${new_version}

# @devs
clean_examples:
# Clean all examples
	cd docs/src; rm -rf example_*.md

# @devs
clone_pull_all:
# Clone and update all FUSE packages
	@ if [ ! -d "$(JULIA_PKG_DEVDIR)" ]; then mkdir -p $(JULIA_PKG_DEVDIR); fi
	@make -i $(PARALLELISM) FUSE GenerateDD $(FUSE_PACKAGES_MAKEFILE)

# @devs
combine_compat_toml: error_missing_repo_var
# Combine Project_PR???.toml files using Julia
	@echo "Combining Project.toml files in repository: $(repo)"
	@cd ../$(repo) && \
	julia -e 'using TOML; \
		project_dir = "./"; \
		project_files = filter(x -> startswith(x, "Project_PR") && endswith(x, ".toml"), readdir(project_dir)); \
		combined_compat = TOML.parsefile("Project.toml")["compat"]; \
		for file in project_files; \
			println(file); \
		    project_path = joinpath(project_dir, file); \
		    project_toml = TOML.parsefile(project_path); \
		    compat_entries = project_toml["compat"]; \
		    for (pkg, version) in compat_entries; \
		        combined_compat[pkg] = version; \
		    end; \
		end; \
		project_toml = TOML.parsefile("Project.toml"); \
		project_toml["compat"] = combined_compat |> sort; \
		project_toml["deps"] = project_toml["deps"] |> sort; \
		open("Project.toml", "w") do io;\
			TOML.print(io, project_toml); \
		end; \
		println("Combined Project.toml saved as Project_combined.toml");'

# @devs
commit_project_toml: resolve
# Commit Project.toml with the new version
# We resolve the environment before committing to ensure that the compat versions are all compatible
	git -C ../$(repo) add Project.toml ;\
	git -C ../$(repo) commit -m "v$(new_version)" ;\
	git -C ../$(repo) push

# @devs
compat: error_missing_repo_var download_compat_tomls combine_compat_toml
# Apply compat patches, combine them in the original Project.toml, and cleanup
	@echo ""
	@echo "To remove temporary Project_PR???.toml files and close CompatHelper PRs, do:"
	@echo ""
	@echo "    make compat_cleanup repo=$(repo)"
	@echo ""

# @devs
compat_cleanup: error_missing_repo_var
# Remove temporary Project_PR???.toml files and close CompatHelper PRs
	@echo "Closing corresponding PRs and deleting temporary Project_PR???.toml files"
	cd ../$(repo) && \
	for file in Project_PR*.toml; do \
		pr_number=$$(echo $$file | sed -n 's/.*Project_PR\([0-9]*\).toml/\1/p'); \
		echo "Closing PR #$$pr_number..."; \
		gh pr close $$pr_number --delete-branch; \
		rm $$file; \
	done && \
	echo "CompatHelper PRs closed and temporary files deleted."

# @devs
daily_example:
# Run daily example to md
	cd docs; julia notebooks_to_md.jl --daily --execute --canfail

ifdef GITHUB_ACTION
# @devs
daily_example_ci_commit:
# Commit daily example md (this must only be run by the CI)
	git config user.email "fuse-bot@fusion.gat.com"
	git config user.name "fuse bot"
	git config push.autoSetupRemote true
	git checkout -b examples_$(TODAY)
	git add -A
	git commit --allow-empty -m "example of the day"
	git push --set-upstream origin examples_$(TODAY)
endif

# @devs
dev_deps_tree:
# Print dependency tree of the packages in dev folder
	@julia -e' ;\
	using Pkg ;\
	Pkg.add("AbstractTrees") ;\
	using AbstractTrees ;\
	function AbstractTrees.printnode(io::IO, uuid::Base.UUID) ;\
		dep = get(Pkg.dependencies(), uuid, nothing) ;\
		print(io, dep.name) ;\
	end ;\
	function AbstractTrees.children(uuid::Base.UUID) ;\
		dep = get(Pkg.dependencies(), uuid, nothing) ;\
		dev_deps = Dict([(key,value) for (key,value) in get(Pkg.dependencies(), uuid, nothing).dependencies if value !== nothing && isdir("../$$(get(Pkg.dependencies(), value, nothing).name)")]) ;\
		tmp= sort!(collect(values(dev_deps)), by=x->get(Pkg.dependencies(), x, (name="",)).name) ;\
	end ;\
	AbstractTrees.print_tree(Pkg.project().dependencies["FUSE"]) ;\
	'

# @devs
develop:
# Develop FUSE and related ProjectTorreyPines packages
	@julia -e '\
fuse_packages = $(FUSE_PACKAGES);\
using Pkg;\
Pkg.activate();\
Pkg.develop([["FUSE"] ; fuse_packages]);\
Pkg.add(["JuliaFormatter", "Test", "Plots"]);\
'
	@make install_revise

# @devs
develop_docs:
# Setup ./docs environment to build documentation
	julia -e '\
	fuse_packages = $(FUSE_PACKAGES);\
	using Pkg;\
	Pkg.activate("./docs");\
	Pkg.develop([["FUSE"] ; fuse_packages]);\
	'

# @devs
devs_update:
# Pull changes for all ProjectTorreyPines packages that are in the .julia/dev folder and resolve environment
	@echo $(DEV_PACKAGES)
	@$(foreach repo,$(DEV_PACKAGES), \
	(sh -c "cd $(JULIA_PKG_DEVDIR)/$(repo) && git pull 2>&1 | sed 's/^/$(repo): /'") & \
	)
	make resolve

# @devs
devs_update_all:devs_update
# Update all dev packages and dependencies
	@julia -e 'using Pkg; Pkg.resolve(); Pkg.update(); Pkg.precompile()'

# @devs
docs: develop_docs
# Generate documentation
	cd docs; julia make.jl

# @devs
download_compat_tomls: error_missing_repo_var
# Apply compat patches and save the resulting Project_PR???.toml files
	@echo "Downloading Project.toml files from CompatHelper PRs in repository: $(repo)"
	@cd ../$(repo) && \
	git fetch && \
	pr_numbers_and_shas=$$(gh pr list --state open --json number,headRefOid,title --jq '.[] | select(.title != null and (.title | contains("CompatHelper:"))) | "\(.number):\(.headRefOid)"') && \
	for pr_info in $$pr_numbers_and_shas; do \
		pr_number=$${pr_info%%:*}; \
		commit_sha=$${pr_info##*:}; \
		echo "Fetching Project.toml from commit $$commit_sha for PR #$$pr_number..."; \
		git show $$commit_sha:Project.toml > Project_PR$${pr_number}.toml; \
	done && \
	echo "Project.toml files downloaded."

# @devs
empty_commit:
# Create an empty commit
	@git reset HEAD
	@git commit --allow-empty -m 'empty commit'

# @devs
error_if_project_toml_is_dirty: error_missing_repo_var
# Utility to error if Project.toml in repo has not been manually modified
	@if ! git -C ../$(repo) diff --quiet -- Project.toml ; then \
		echo "Error: Project.toml has uncommitted changes." ;\
		exit 1 ;\
	fi

# @devs
error_missing_repo_var:
# Throw an error if environment variable `repo` is not set
ifeq ($(repo),)
	$(error `repo` variable is not set)
endif

# @devs
error_missing_repos_var:
# Throw an error if environment variable `repos` is not set
ifeq ($(repo),)
	$(error `repos` variable is not set)
endif

# @devs
error_on_previous_commit_is_a_version_bump: error_missing_repo_var
# Utility to error if previous commit was a version bump
	@previous_commit_message=$$(git -C ../$(repo) log -1 --pretty=%B) ;\
	if echo "$$previous_commit_message" | grep -Eq '^v[0-9]+\.[0-9]+\.[0-9]+' ; then \
		echo "Error: The previous commit was a version bump." ;\
		exit 1 ;\
	fi

# @devs
feature_or_master:
# checks if on the packages remote GitHub repos there is a branch with the same name of the local FUSE branch
	julia -e ';\
	$(feature_or_master_julia);\
	fuse_packages = $(FUSE_PACKAGES);\
	for package in fuse_packages;\
		branch = feature_or_master(package, "$(FUSE_LOCAL_BRANCH)");\
        if branch == "master";\
            println(">>> $$(package)");\
        else;\
            println(">>> $$(package) @ $$(branch)");\
        end;\
	end'

# @devs
generate_dd:
# Update dd from the json files in IMASdd
	@julia -e 'using GenerateDD; update_data_structures_from_OMAS(); generate_dd()'

# @devs
help_devs:
# Print developer makefile commands help
	@FILTER="devs" make help_common

# @devs
help_info:
	@printf "\n"
	@printf ">> Use \`fusebot help_user\` to get the users' list of commands\n"
	@printf ">> Use \`fusebot help_devs\` to get the developers' list of commands\n"
	@printf "\n"

# @devs
help_user:
# Print user makefile commands help
	@FILTER="user" make help_common

header:
	@printf "\n"
	@printf "  \033[1;31m███████\033[1;30m╗\033[1;31m██\033[1;30m╗   \033[1;31m██\033[1;30m╗\033[1;31m███████\033[1;30m╗\033[1;31m███████\033[1;30m╗\033[0m\n"
	@printf "  \033[1;31m██\033[1;30m╔════╝\033[1;31m██\033[1;30m║   \033[1;31m██\033[1;30m║\033[1;31m██\033[1;30m╔════╝\033[1;31m██\033[1;30m╔════╝\033[0m\n"
	@printf "  \033[1;31m█████\033[1;30m╗  \033[1;31m██\033[1;30m║   \033[1;31m██\033[1;30m║\033[1;31m███████\033[1;30m╗\033[1;31m█████\033[1;30m╗  \033[0m\n"
	@printf "  \033[1;31m██\033[1;30m╔══╝  \033[1;31m██\033[1;30m║   \033[1;31m██\033[1;30m║╚════\033[1;31m██\033[1;30m║\033[1;31m██\033[1;30m╔══╝  \033[0m\n"
	@printf "  \033[1;31m██\033[1;30m║     ╚\033[1;31m██████\033[1;30m╔╝\033[1;31m███████\033[1;30m║\033[1;31m███████\033[1;30m╗\033[0m\n"
	@printf "  \033[1;30m╚═╝      ╚═════╝ ╚══════╝╚══════╝\033[0m\n"
	@printf "   Project  Torrey  Pines  (PTP)\n"

init_expressions:
# Generates init_expressions.json file, which lists entries that are
# always expected to be expressions when coming out of init()
	julia -e 'import FUSE; FUSE.init_expressions(;save=true)'

# @devs
install: install_no_registry
# Install with default install method (no registry)

# @devs
install_PyCall:
# Install PyCall
	julia -e '\
	ENV["PYTHON"]="";\
	using Pkg;\
	Pkg.add("PyCall");\
	Pkg.build("PyCall");\
	'

# @devs
install_ci_add:
# Install (add) FUSE via HTTPS and $PTP_READ_TOKEN
# Looks for same branch name for all repositories otherwise falls back to master
	julia -e ';\
	$(feature_or_master_julia);\
	fuse_packages = $(FUSE_PACKAGES);\
	using Pkg;\
	Pkg.activate(".");\
	dependencies = Pkg.PackageSpec[];\
	for package in fuse_packages;\
		branch = feature_or_master(package, "$(FUSE_LOCAL_BRANCH)");\
        if branch == "master";\
            println(">>> $$(package)");\
			push!(dependencies, Pkg.PackageSpec(url="https://project-torrey-pines:$(PTP_READ_TOKEN)@github.com/ProjectTorreyPines/"*package*".jl.git"));\
        else;\
            println(">>> $$(package) @ $$(branch)");\
			push!(dependencies, Pkg.PackageSpec(url="https://project-torrey-pines:$(PTP_READ_TOKEN)@github.com/ProjectTorreyPines/"*package*".jl.git", rev=branch));\
        end;\
	end;\
	Pkg.add(dependencies);\
	Pkg.status()'

# @devs
install_ci_dev: install_registry https_dev
# Install used by CI (dev packages, do not add them)

# @devs
install_examples_dev:
# Install FUSE_examples under FUSE/examples
	@if [ ! -d "examples" ]; then git clone git@github.com:ProjectTorreyPines/FUSE_examples.git examples ; else cd examples && git pull; fi

# @devs
install_no_registry: install_registry clone_pull_all develop
# Install FUSE without using the registry

# @devs
install_playground: .PHONY
# Clone FUSE_playground repository under FUSE/playground folder
	if [ -d playground ] && [ ! -f playground/.gitattributes ]; then mv playground playground_private ; fi
	if [ ! -d "playground" ]; then git clone git@github.com:ProjectTorreyPines/FUSE_playground.git playground ; else cd playground && git pull origin `git rev-parse --abbrev-ref HEAD` ; fi

# @devs
install_via_registry: install_registry develop
# Install FUSE using the registry

# @devs
list_open_compats:
# List compat patches PR on GitHub
	@$(foreach repo,$(DEV_PACKAGES), \
		echo ;\
		echo $(repo) ;\
		cd ../$(repo) && \
		git fetch 2> /dev/null && \
		pr_numbers_shas_and_titles=$$(gh pr list --state open --json number,headRefOid,title --jq '.[] | select(.title != null and (.title | contains("CompatHelper:"))) | "\(.number)@\(.headRefOid)@\(.title)"') && \
		IFS=$$'\n'; \
		for pr_info in $$pr_numbers_shas_and_titles; do \
			pr_number=$${pr_info%%@*}; \
			temp=$${pr_info#*@}; \
			commit_sha=$${temp%%@*}; \
			title=$${temp#*@}; \
			echo "$(repo): PR #$$pr_number - $$title"; \
		done; \
	)

# @devs
apache: error_missing_repo_var
# Update LICENSE, NOTICE.md, github workflows, docs, juliaformatter and gitignore in preparation of public release
# The starting information is taken from IMASdd.jl and moved to the target repo
# >> make apache repo=CHEASE
# in addition, one must add the DOCUMENTER_KEY to the repo
# https://m3g.github.io/JuliaNotes.jl/stable/publish_docs/#How-to-deploy-the-documentation-of-a-project
	@echo $(repo)
	@cp ../IMASdd/LICENSE ../$(repo)/ ;\
	\
	cp ../IMASdd/NOTICE.md ../$(repo)/ ;\
	sed -i.bak "s/IMASdd/$(repo)/g" ../$(repo)/NOTICE.md && rm ../$(repo)/NOTICE.md.bak ;\
	\
	mkdir -p ../$(repo)/.github/workflows ;\
	cp ../IMASdd/.github/workflows/make_docs.yml ../$(repo)/.github/workflows/ ;\
	cp ../IMASdd/.github/workflows/runtests.yml ../$(repo)/.github/workflows/ ;\
	cp ../IMASdd/.github/workflows/CompatHelper.yml ../$(repo)/.github/workflows/ ;\
	cp ../IMASdd/.github/workflows/TagBot.yml ../$(repo)/.github/workflows/ ;\
	\
	cp ../IMASdd/README.md ../$(repo)/ ;\
	sed -i.bak "s/IMASdd/$(repo)/g" ../$(repo)/README.md && rm ../$(repo)/README.md.bak ;\
	\
	cp -R ../IMASdd/docs ../$(repo)/ ;\
	sed -i.bak "s/IMASdd/$(repo)/g" ../$(repo)/docs/make.jl && rm ../$(repo)/docs/make.jl.bak ;\
	echo "# $(repo).jl" > ../$(repo)/docs/src/index.md ;\
	rm ../$(repo)/docs/Manifest.toml ;\
	rm -rf ../$(repo)/docs/build ;\
	\
	cp -R ../IMASdd/.JuliaFormatter.toml ../$(repo)/ ;\
	\
	cp -R ../IMASdd/.gitignore ../$(repo)/ ;\
	\
	julia -e 'import Pkg; Pkg.add("DocumenterTools"); import DocumenterTools; DocumenterTools.genkeys()'

# @devs
nuke_julia:
# Remove everything under $HOME/.julia besides $HOME/.julia/dev
	mv $(JULIA_PKG_DEVDIR) $(call realpath,$(JULIA_DIR))/../asddsaasddsa
	rm -rf $(JULIA_DIR)
	mkdir -p $(JULIA_DIR)
	mv $(call realpath,$(JULIA_DIR))/../asddsaasddsa $(JULIA_PKG_DEVDIR)

# @devs
register: error_missing_repo_var
# Register a package to FuseRegistry
# >> make register repo=IMASdd
	@current_branch=$(shell git -C ../$(repo) rev-parse --abbrev-ref HEAD) ;\
	if [ "$$current_branch" != "master" ]; then \
		echo "Error: $(repo) is not on the master branch" ;\
		exit 1 ;\
	fi
	sed -i.bak "s/https:\/\/github.com\//git@github.com:/g" $(JULIA_DIR)/registries/FuseRegistry/.git/config && rm $(JULIA_DIR)/registries/FuseRegistry/.git/config.bak
	julia -e '\
using Pkg;\
Pkg.Registry.update("FuseRegistry");\
Pkg.activate();\
using LocalRegistry;\
LocalRegistry.is_dirty(path, gitconfig)= false; register("$(repo)", registry="FuseRegistry")'
	version=$$(grep '^version' ../$(repo)/Project.toml | sed -E 's/version = "(.*)"/\1/') ;\
	git -C ../$(repo) tag -a v$${version} -m "v$${version}" ;\
	git -C ../$(repo) push origin v$${version} ;\
	gh release create v$${version} --generate-notes --repo ProjectTorreyPines/$(repo).jl

# @devs
register_major: bump_major register
# Bump major version and register: X.Y.Z --> (X+1).0.0

# @devs
register_minor: bump_minor register
# Bump minor version and register: X.Y.Z --> X.(Y+1).0

# @devs
register_patch: bump_patch register
# Bump patch version and register: X.Y.Z --> X.Y.(Z+1)

# @devs
resolve:
# Resolve the current environment (eg. after manually adding a new package)
	@julia -e 'using Pkg; Pkg.resolve(); Pkg.precompile()'

# @devs
rm_manifests:
# Remove all Manifest.toml files in the .julia/dev folder
	@find ..  -maxdepth 3 -type f -name "Manifest.toml" -exec rm -f \{\} \;

# @devs
run_examples: clean_examples
# Clean and run all examples
	cd docs; julia notebooks_to_md.jl --execute

# @devs
status:
# List branches of all the ProjectTorreyPines packages used by FUSE with version, dirty * flag, and commits since the latest tag
	@cd $(CURRENTDIR); \
	packages="$(DEV_PACKAGES)"; \
	sorted_packages=`echo $$packages | tr ' ' '\n' | sort | tr '\n' ' '`; \
	line_count=0; \
	term_width=`tput cols`; \
	for package in $$sorted_packages; do \
		package_dir="../$$package"; \
		if [ ! -f $$package_dir/Project.toml ]; then \
			continue; \
		fi; \
		line_count=$$((line_count + 1)); \
		if [ $$((line_count % 2)) -eq 0 ]; then \
			color="\033[47m"; \
		else \
			color="\033[46m"; \
		fi; \
		reset="\033[0m"; \
		branch=`cd $$package_dir && git rev-parse --abbrev-ref HEAD`; \
		version=`grep -m1 'version =' $$package_dir/Project.toml | awk -F' = ' '{print $$2}' | tr -d '"'`; \
		dirty=`cd $$package_dir && [ -n "$$(git status --porcelain)" ] && echo "dirty" || echo "    "`; \
		latest_tag=`cd $$package_dir && git describe --tags --abbrev=0 2>/dev/null || echo "(no tag)"`; \
		commits_since_tag=""; \
		if [ "$$latest_tag" != "(no tag)" ]; then \
			commits_since_tag=`cd $$package_dir && git log $$latest_tag..HEAD --oneline | wc -l | tr -d ' '`; \
		fi; \
		commit_info=""; \
		if [ -n "$$commits_since_tag" ] && [ $$commits_since_tag -gt 0 ]; then \
			commit_info="($$commits_since_tag commits since latest version tag)"; \
		elif [ "$$latest_tag" = "(no tag)" ]; then \
			commit_info="(no tag)"; \
		fi; \
		line_text=`printf "%25s %10s @ %-15s %-10s %s" "$$package" "$$version" "$$branch" "$$dirty" "$$commit_info"`; \
		line_length=`echo "$$line_text" | wc -c | tr -d ' '`; \
		padding=$$((term_width - line_length)); \
		printf "$$color%s%*s$$reset\n" "$$line_text" $$padding ""; \
	done

# @devs
versions_used: error_missing_repo_var
# Search for [compat] statements in upstream packages
	@echo
	@echo "Check [compat] statements in the Project.toml of the following repos:"
	@echo
	@grep "$$repo = \"" ../*/Project.toml | grep -v -E '[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}' | sed -E 's|../([^/]+)/Project.toml:.*=.*"([^"]+)"|  \1: \2|'
	@echo

# @devs
web:
# Push documentation to the web (user entry point)
	if [ ! -d "$(PWD)/docs/pages" ]; then cd docs; git clone --single-branch -b gh-pages git@github.com:ProjectTorreyPines/FUSE.jl.git pages; fi
	make web_push

# @devs
web_ci:
# Push documentation to the web (CI entry point)
	git clone $(PWD) docs/pages
	cp .git/config docs/pages/.git/config
	cd docs/pages; git fetch; git checkout gh-pages
	cd docs/pages; git config user.email "fuse@fusion.gat.com"
	cd docs/pages; git config user.name "FUSE-BOT"
	make web_push

# @devs
web_push:
# Push documentation to the web
	cd docs/pages; git reset --hard 049da2c703ad7fc552c13bfe0651da677e3c7f58
	cd docs; cp -rf build/* pages/dev/
	cd docs/pages; echo "fuse.help" > CNAME ### this is to set the custom domain name for gh-pages
	cd docs/pages; touch .nojekyll
	cd docs/pages; git add -A; git commit --allow-empty -m "documentation"; git push --force

# @user
help_common:
# Common function for help command, differentiating between @user and @devs and rejecting anything above the marker
	@awk ' \
		BEGIN { \
			show = ENVIRON["FILTER"]; \
			found_marker = 0; \
		} \
		/^# =%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%= #/ { found_marker = 1; next } \
		!found_marker { next } \
		/^[#] @user/ { show = ("user" == ENVIRON["FILTER"]); next } \
		/^[#] @devs/ { show = ("devs" == ENVIRON["FILTER"]); next } \
		/^[a-zA-Z0-9_.-]+:/ { \
			if (show) { \
				if (target != "" && docstring != "") { \
					printf "\033[1m" target "\033[0m"; \
					if (dependencies != "") { \
						printf " : \033[36m" dependencies "\033[0m"; \
					} \
					print ""; \
					print docstring; \
					print ""; \
				} \
				target = $$1; \
				dependencies = ""; \
				if (NF > 1) { \
					dependencies = $$2; \
					for (i = 3; i <= NF; i++) { \
						dependencies = dependencies " " $$i; \
					} \
				} \
				docstring = ""; \
			} \
			next; \
		} \
		/^[ \t]*#/ { \
			if (show && target != "") { \
				comment = substr($$0, match($$0, /#/) + 1); \
				gsub(/^[ \t]+/, "", comment); \
				if (docstring != "") { \
					docstring = docstring "\n    " comment; \
				} else { \
					docstring = "    " comment; \
				} \
			} \
		} \
		END { \
			if (show && target != "" && docstring != "") { \
				printf "\033[1m" target "\033[0m"; \
				if (dependencies != "") { \
					printf " : \033[36m" dependencies "\033[0m"; \
				} \
				print ""; \
				print docstring; \
				print ""; \
			} \
		}' $(MAKEFILE_LIST)

# @user
install_IJulia:
# Install IJulia
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

# @user
install_examples:
# Install FUSE_examples in current folder
	@cd $(PTP_ORIGINAL_DIR) && if [ ! -d "FUSE_examples" ]; then git clone git@github.com:ProjectTorreyPines/FUSE_examples.git ; else cd FUSE_examples && git pull; fi

# @user
install_registry:
# Add the FuseRegistry to the list of Julia registries
	julia -e 'using Pkg; Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git")); Pkg.Registry.add("General");'

# @user
install_revise:
# Setup Revise.jl to automatically be used when Julia starts
	@echo "Setting Revise.jl to run at startup"
	@julia -e 'using Pkg; Pkg.add("Revise")'
	@mkdir -p $(JULIA_DIR)/config
	@touch $(JULIA_CONF)
	@grep -v -F -x "using Revise" "$(JULIA_CONF)" > "$(JULIA_CONF).tmp" || true
	@echo "using Revise" | cat - "$(JULIA_CONF).tmp" > "$(JULIA_CONF)"
	@rm -f "$(JULIA_CONF).tmp"

# @user
self_update:
# updates the `fusebot` executable to the latest version
	@if cmp -s fusebot `which fusebot`; then \
		echo "fusebot is already up to date"; \
	else \
		cp fusebot `which fusebot`; \
		echo "fusebot has been updated"; \
	fi

# @user
threads:
# Simple test to see how many threads julia will run on (set by JULIA_NUM_THREADS)
	@echo "set the JULIA_NUM_THREADS environment variable"
	@julia -e "println(Threads.nthreads())"

# @user
update:
# Update FUSE and its dependencies
	@julia -e 'using Pkg; Pkg.resolve(); Pkg.update("FUSE"; preserve=Pkg.Types.PreserveLevel(2)); Pkg.precompile()'
