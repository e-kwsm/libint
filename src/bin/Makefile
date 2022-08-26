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
	    $(MAKE) -C $${dir} $(DODEPENDOPT) install; \
	  done

install_inc::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) install_inc; \
	  done

install_target::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) install_target; \
	  done

uninstall::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) uninstall; \
	  done

clean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) clean; \
	  done

oclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) oclean; \
	  done

distclean::
	for dir in $(ALLSUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) distclean; \
	  done
	-rm -rf MakeVars

targetclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) targetclean; \
	  done

realclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) realclean; \
	  done

