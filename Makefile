# library name
lib.name = quatkram

class.sources = qacc~.c qvrot~.c faccwrap~.c faccbounce~.c noisen~.c mc_conv~.c

datafiles = \
	noisen~-help.pd \
	mc_wrap-help.pd \
	${empty}

objectsdir = ./build
PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
