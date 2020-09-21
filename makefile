# Copyright (C) 2015--2020 Bhalchandra Digambar Thatte

# This file is a part of the project ASSIGN - a program to allocate
# disciplinas to professors. ASSIGN project is released under the GPL V3 license.

all: assign

OBJS = assign.o

CCFLAGS = -std=gnu++11 -Wall -static -O3

clean:
	-$(RM) $(all) $(OBJS) *~

.c.o:
	gcc -Wall -I. -c $<

.cc.o:
	g++ $(CCFLAGS) -I. -c $<

assign: assign.o
	g++ $(CCFLAGS) -L. assign.o -o assign

# install: assign
# 	cp assign ~/.local/bin/assign
