R=R
# -> you can do    R=R-devel  make ....

PACKAGE=surveyIndex
VERSION := $(shell sed -n '/^Version: /s///p' surveyIndex/DESCRIPTION)
DATE := $(shell sed -n '/^Date: /s///p' surveyIndex/DESCRIPTION)
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip

SUBDIRS := $(wildcard testmore/*/.)

.PHONY: all testmoreseq testonemore testmore $(SUBDIRS)

all:
	make doc-update
	make build-package
	make install
	make pdf

doc-update:
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\")" | R --slave

build-package:
	R CMD build --resave-data=no $(PACKAGE)

install:
	make build-package
	R CMD INSTALL --preclean $(TARBALL)

unexport TEXINPUTS
pdf:
	rm -f $(PACKAGE).pdf
	R CMD Rd2pdf --no-preview $(PACKAGE)

check:
	R CMD check $(PACKAGE)

##unlock:
##	rm -rf `Rscript --vanilla -e 'writeLines(.Library)'`/00LOCK-

##.PHONY: dox
##dox:
##	sed -i s/^PROJECT_NUMBER.*/PROJECT_NUMBER=v$(VERSION)/g dox/Doxyfile
##	cd dox; doxygen

##dox-clean:
##	cd dox; rm -rf html latex


## Get a rough changelog since most recent github revision tag
## (Use as starting point when updating NEWS file)
## NOTE: Run *after* updating version and date in DESCRIPTION.
changelog:
	echo; \
	echo "------------------------------------------------------------------------"; \
	echo TMB $(VERSION) \($(DATE)\); \
	echo "------------------------------------------------------------------------"; \
	echo; \
	git --no-pager log --format="o %B" `git describe --abbrev=0 --tags`..HEAD | sed s/^-/\ \ -/g

NPROCS:=1
OS:=$(shell uname -s)

ifeq ($(OS),Linux)
  NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)
endif

testmore:
	$(MAKE) -j $(NPROCS) testmoreseq

testmoreseq: $(SUBDIRS)

testonemore:
	@$(MAKE) testmore/$(ARG)/.

$(SUBDIRS):
	@cp testmore/Makefile $@
	@$(MAKE) -i -s -C $@
	@rm -f $@/Makefile
