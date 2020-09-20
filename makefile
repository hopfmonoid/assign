# Copyright (C) 2015--2020 Bhalchandra Digambar Thatte

# This file is a part of the project ASSIGN - a program to allocate
# disciplinas to professors

# Author: Bhalchandra Digambar Thatte

# ASSIGN is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# ASSIGN is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
