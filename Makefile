
MPICC:=mpicc

CFLAGS:=
CPPFLAGS:=

#CFLAGS+=-ggdb

SUBDIRS:=atasp/ z_pack/ zmpi_local/ zmpi_tools/

HFILES:=$(wildcard $(addsuffix *.h,$(SUBDIRS)))
CFILES:=$(wildcard $(addsuffix *.c,$(SUBDIRS)))
OFILES:=$(addsuffix .o,$(basename $(CFILES)))

LIB_A:=libzmpi_atasp.a

DEMO:=specific_demo


$(LIB_A): $(OFILES)
	ar r $@ $^

z_pack/%.o: z_pack/%.c $(HFILES) Makefile
	$(MPICC) $(CPPFLAGS) -Iz_pack/ -c -o $@ $<

zmpi_local/%.o: zmpi_local/%.c $(HFILES) Makefile
	$(MPICC) $(CPPFLAGS) -Iz_pack/ -c -o $@ $<

zmpi_tools/%.o: zmpi_tools/%.c $(HFILES) Makefile
	$(MPICC) $(CPPFLAGS) -Iz_pack/ -c -o $@ $<

atasp/%.o: atasp/%.c $(HFILES) Makefile
	$(MPICC) $(CPPFLAGS) -Iz_pack/ -Izmpi_local/ -DHAVE_ZMPI_LOCAL_H -Izmpi_tools/ -DHAVE_ZMPI_TOOLS_H -Iatasp/ -c -o $@ $<

demo: $(DEMO)

$(DEMO): demo/specific_demo.c $(LIB_A)
	$(MPICC) $(CPPFLAGS) -I./ -Izmpi_local/ -DHAVE_ZMPI_LOCAL_H $(CFLAGS) -o $@ $< $(LIB_A) -lm

clean: 
	rm -f $(OFILES) $(LIB_A) $(DEMO)
