# Generic Makefile, courtesy of Brian Farris 2013
#
# THERE SHOULD BE NO NEED TO EDIT THIS MAKEFILE.
#
# ALL MACHINE SPECIFIC FLAGS/PATHS ARE SET IN Makefile.in

MAKEFILE_IN = $(PWD)/Makefile.in
include $(MAKEFILE_IN)

APP      = hydro2d

SRCEXT   = c
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)')
SRCDIRS := $(shell find . -name '*.$(SRCEXT)' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))

GIT_VERSION = $(shell git describe --dirty --always --tags)

#Vital flags 
DEBUG    = -g
CFLAGS   += -O3 -Wall -c $(DEBUG) -DVERSION=\"$(GIT_VERSION)\"
LDFLAGS  = -lm

#HDF5 Installation - specified in Makefile.in
CFLAGS += -I$(H5DIR)/include
LDFLAGS += -L$(H5DIR)/lib -lhdf5

.PHONY: all clean distclean


all: $(BINDIR)/$(APP)

$(BINDIR)/$(APP): buildrepo $(OBJS)
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) $(OBJS) $(LDFLAGS) -o $@

$(OBJDIR)/%.o: %.$(SRCEXT)
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

clean:
	$(RM) -r $(OBJDIR)

distclean: clean
	$(RM) -r $(BINDIR)

buildrepo:
	@$(call make-repo)

define make-repo
   for dir in $(SRCDIRS); \
   do \
	mkdir -p $(OBJDIR)/$$dir; \
   done
endef
