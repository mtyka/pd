# PD is a free, modular C++ library for biomolecular simulation with a 
# flexible and scriptable Python interface. 
# Copyright (C) 2003-2013 Mike Tyka and Jon Rea
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#!/usr/bin/python

from pd import *
import sys

if len(sys.argv) < 3:
  print "Syntax: pickle2mime infile outfile"
  sys.exit(0)


infile = sys.argv[1]
outfile = sys.argv[2]

psp = readOldPickle(infile)
psp.writeMIME(outfile)

