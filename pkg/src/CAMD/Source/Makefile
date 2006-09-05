#-------------------------------------------------------------------------------
# CAMD Makefile for compiling on Unix systems (for GNU make only)
#-------------------------------------------------------------------------------

include $(MkInclude)

PKG_CFLAGS = -DDINT -I../Include -I../../UFconfig

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

CAMD = camd_aat camd_1 camd_2 camd_dump camd_postorder camd_defaults \
	camd_order camd_control camd_info camd_valid camd_preprocess

CAMDI = $(addsuffix .o, $(subst amd_,amd_i_,$(CAMD)))

UFCONFIG = ../../UFconfig/UFconfig.h

INC = ../Include/camd.h ../Include/camd_internal.h $(UFCONFIG)

LIB = ../../CAMD.a

lib: $(LIB)

$(LIB): camd_global.o ${CAMDI}
	$(AR) -rucs $(LIB) camd_global.o ${CAMDI}

camd_global.o: camd_global.c $(INC)
camd_i_%.o: camd_%.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -c $< -o $@

mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o 
