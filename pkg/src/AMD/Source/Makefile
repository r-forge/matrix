include $(MkInclude)

PKG_CFLAGS = -I../Include -I../../UFconfig

AMD = amd_aat amd_1 amd_2 amd_postorder amd_post_tree amd_defaults \
	amd_order amd_control amd_info amd_valid amd_preprocess amd_dump

AMDI = $(addsuffix .o, $(subst amd_,amd_i_,$(AMD)))
AMDL = $(addsuffix .o, $(subst amd_,amd_l_,$(AMD)))

UFCONFIG = ../../UFconfig/UFconfig.h

INC = ../Include/amd.h ../Include/amd_internal.h $(UFCONFIG)

LIB = ../../AMD.a

lib: $(LIB)

$(LIB): amd_global.o $(AMDI) $(AMDL)
	$(AR) -rucs $(LIB) amd_global.o $(AMDI) $(AMDL)

amd_global.o: amd_global.c $(INC)
amd_i_%.o: amd_%.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDINT -c $< -o $@
amd_l_%.o: amd_%.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -DDLONG -c $< -o $@

mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o 
