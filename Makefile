lib.name = quatkram

class.sources = \
	acos~.c \
	asin~.c \
	atan~.c \
	atan2~.c \
	faccbounce~.c \
	faccleak~.c \
	faccwrap~.c \
	frft.c \
	frft~.c \
	matrixfb~.c \
	mc_conv~.c \
	mc_conv2d~.c \
	mc_route~.c \
	nchans~.c \
	noisen~.c \
	qacc~.c \
	qdiv~.c \
	qmul~.c \
	qrsqrt~.c \
	qvtrans~.c \
	randomwalk_space~.c \
	randomwalk_sphere~.c \
	sampdel~.c \
	tabloop~.c \
	tabredraw.c \
	tabsmear~.c \
	tan~.c \
	tanh~.c \
	tetra2pos.c \
	urn~.c \
	zc~.c \
	zcflip~.c \
	${empty}

abstractions = \
	gt~.pd \
	gte~.pd \
	hip~.cl.pd \
	lop~.cl.pd \
	lt~.pd \
	lte~.pd \
	mc_abs~.pd \
	mc_hip~.pd \
	mc_lop~.pd \
	qconj~.pd \
	qdiv~-help.pd \
	qfromgyroaccel~.pd \
	qinv~.pd \
	qmag~.pd \
	qnorm~.pd \
	qnormalize~.pd \
	sign~.pd \
	tanh~.pd \
	tetra2pos_abs.pd \
	urn~-help.pd \
	zc~-help.pd \
	${empty}

helpfiles = \
	acos~-help.pd \
	asin~-help.pd \
	atan~-help.pd \
	atan2~-help.pd \
	faccbounce~-help.pd \
	faccleak~-help.pd \
	faccwrap~-help.pd \
	mc_conv~-help.pd \
	mc_conv2d~-help.pd \
	mc_route~-help.pd \
	matrixfb~-help.pd \
	noisen~-help.pd \
	qacc~-help.pd \
	qrsqrt~-help.pd \
	randomwalk_sphere~-help.pd \
	sampdel~-help.pd \
	qmul~-help.pd \
	tan~-help.pd \
	tetra2pos-help.pd \
	${empty}

datafiles = $(abstractions) $(helpfiles)

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

# Suppress all warnings for a clean build
# cflags += -w

# Add FFTW dependency
ldlibs += -lfftw3

objectsdir = ./build
PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
