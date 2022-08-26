TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

-include $(TOPDIR)/Makedirlist

SUBDIRS = libint

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
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

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
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done
	-rm -rf MakeRules MakeVars

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

