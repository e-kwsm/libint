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
	    $(MAKE) -C $${dir} $(JOBS) || exit 1; \
	  done

ifndef DODEPEND
DODEPENDOPT = "DODEPEND=no"
endif

install::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

install_target::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
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
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
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
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

clean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done

distclean::
	for dir in $(SUBDIRS); \
	  do \
	    $(MAKE) -C $${dir} $(DODEPENDOPT) $@ || exit 1; \
	  done
	-rm -f MakeRules MakeVars

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

