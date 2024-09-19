# library name
lib.name = quatkram

class.sources = qacc~.c qvrot~.c faccwrap~.c faccbounce~.c noisen~.c

objectsdir = ./build
PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
