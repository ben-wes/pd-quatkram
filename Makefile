# library name
lib.name = quaternion

class.sources = qacc~.c qvrot~.c facc~.c

objectsdir = ./build
PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
