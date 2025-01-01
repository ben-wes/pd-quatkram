lib.name = quatkram

class.sources = \
	atan2~.c \
	faccbounce~.c \
	faccwrap~.c \
	faccleak~.c \
	nchans~.c \
	noisen~.c \
	qacc~.c \
	qmul~.c \
	qdiv~.c \
	qvtrans~.c \
	urn~.c \
	mc_route~.c \
	zc~.c \
	tetra2pos.c \
	tabsmear~.c \
	tabredraw.c \
	tabloop~.c \
	sampdel~.c \
	frft~.c \
	mc_conv~.c \
	mc_conv2d~.c \
	${empty}

datafiles = \
	qconj~.pd \
	qinv~.pd \
	qnorm~.pd \
	qnormalize~.pd \
	qmag~.pd \
	qfromgyroaccel~.pd \
	zc~-help.pd \
	tetra2pos-help.pd \
	sampdel~-help.pd \
	atan2~-help.pd \
	faccbounce~-help.pd \
	faccwrap~-help.pd \
	faccleak~-help.pd \
	mc_conv~-help.pd \
	mc_conv2d~-help.pd \
	mc_route~-help.pd \
	noisen~-help.pd \
	urn~-help.pd \
	qacc~-help.pd \
	qmul~-help.pd \
	qdiv~-help.pd \
	${empty}

define forDarwin
  cflags += -I/opt/homebrew/include
  ldflags += -L/opt/homebrew/lib
endef

define forWindows
  cflags += -IC:/msys64/mingw64/include
  ldflags += -LC:/msys64/mingw64/lib
endef

define forLinux
  cflags += -I/usr/include
  ldflags += -L/usr/lib
endef

# Add FFTW dependency
ldlibs += -lfftw3

objectsdir = ./build
PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
