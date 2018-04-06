#  
#  Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 Michael Hofmann, Chemnitz University of Technology
#  
#  This file is part of the ZMPI All-to-all Specific Library.
#  
#  The ZMPI-ATASP is free software: you can redistribute it and/or
#  modify it under the terms of the GNU Lesser Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  The ZMPI-ATASP is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser Public License for more details.
#  
#  You should have received a copy of the GNU Lesser Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  


MPICC:=mpicc

CFLAGS:=
CPPFLAGS:=-DHAVE_Z_CONFIG_H

CFLAGS+=-Wall
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
