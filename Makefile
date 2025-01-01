lib.name = quatkram

# Set FFTW paths explicitly for macOS/Homebrew
FFTW_INCLUDE = /opt/homebrew/include
FFTW_LIB = /opt/homebrew/lib

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
	# add nchans~ help

# List externals that need FFTW3
FFTW_SOURCES = \
	frft~.c \
	mc_conv~.c \
	mc_conv2d~.c \


# Add FFTW flags only to specific externals
$(FFTW_SOURCES:.c=.class.sources.o): cflags += -I$(FFTW_INCLUDE)
$(FFTW_SOURCES:.c=.class.sources.o): ldflags += -L$(FFTW_LIB)
$(FFTW_SOURCES:.c=.class.sources.o): ldlibs += -lfftw3

objectsdir = ./build
PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
