TOPDIR=..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

-include $(TOPDIR)/Makedirlist
-include $(TOPDIR)/doc/MakeVars

SUBDIRS = classdoc progman

default::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(JOBS); \
	  done

ifndef DODEPEND
DODEPENDOPT = "DODEPEND=no"
endif

install::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

install_target::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

install-pdf:: pdf
	$(MAKE) -C progman install

install-ps:: ps
	$(MAKE) -C progman install

install-dvi:: dvi
	$(MAKE) -C progman install

install-html:: html
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

pdf::
	$(MAKE) -C progman $@

ps::
	$(MAKE) -C progman $@

dvi::
	$(MAKE) -C progman $@

html::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

clean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done

distclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@; \
	  done
	-rm -f MakeRules MakeVars

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

