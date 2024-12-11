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
	mc_conv~.c \
	mc_conv2d~.c \
	noisen~.c \
	qacc~.c \
	qmul~.c \
	qdiv~.c \
	qvtrans~.c \
	urn~.c \
	mc_route~.c \
	zc~.c \
	frft~.c \
	ambi2dir.c \
	tetra2pos.c \
	${empty}

datafiles = \
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
	qconj~.pd \
	qinv~.pd \
	qnorm~.pd \
	qnormalize~.pd \
	qmag~.pd \
	qfromgyroaccel~.pd \
	zc~-help.pd \
	tetra2pos-help.pd \
	${empty}
	# add nchans~ help

# FFTW3 specific flags
cflags += -I$(FFTW_INCLUDE)
ldflags += -L$(FFTW_LIB)
ldlibs += -lfftw3

objectsdir = ./build
PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
