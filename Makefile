# library name
lib.name = quatkram

class.sources = \
	atan2~.c \
	faccbounce~.c \
	faccwrap~.c \
	mc_conv~.c \
	mc_conv2d~.c \
	noisen~.c \
	qacc~.c \
	qvtrans~.c \
	urn~.c \
	mc_route~.c \
	${empty}

datafiles = \
	atan2~-help.pd \
	mc_conv~-help.pd \
	mc_conv2d~-help.pd \
	mc_route~-help.pd \
	noisen~-help.pd \
	urn~-help.pd \
	${empty}


objectsdir = ./build
PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
