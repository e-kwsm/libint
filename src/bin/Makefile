TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

-include $(TOPDIR)/Makedirlist

SUBDIRS = libint
ALLSUBDIRS = $(SUBDIRS) test_eri profile

default::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(JOBS); \
	  done

all:: default

ifndef DODEPEND
DODEPENDOPT = "DODEPEND=no"
endif

export::
	$(MAKE) -C libint $(DODEPENDOPT)

install::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

install_inc::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

install_target::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

uninstall::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

clean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

oclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

distclean::
	for dir in $(ALLSUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done
	-rm -rf MakeVars

targetclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

realclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

