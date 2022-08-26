TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

-include $(TOPDIR)/Makedirlist

SUBDIRS = libint

default::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(JOBS) || exit 1; \
	  done

all:: default

ifndef DODEPEND
DODEPENDOPT = "DODEPEND=no"
endif

export::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

install::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

install_inc::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

install_target::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

uninstall::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

clean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

oclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

distclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done
	-rm -rf MakeRules MakeVars

targetclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

realclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

