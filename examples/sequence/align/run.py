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

from pd import *

# 1) Use a 'ForcefieldParameterSet' for naming definitions

ffps = ForcefieldParameterSet()
ffps.readlib("lib/amber03aa.ff")

s1 = BioSequence(ffps)
s1.ParseSequenceString("*s-ASP-GLU-AlA-Gly-ASP-D-GLU-ALA-GLU-A*")
flushed_print("Sequence 1: '")
s1.PrintToScreen()
flushed_print("'\n")

# 2) Use the 'MoleculeNaming' library set for naming definitions

molName = MoleculeNaming()

s2 = BioSequence(molName)
s2.ParseSequenceString("GLU-AlA-G-ASP-D-AlA-GLU")
flushed_print("Sequence 2: '")
s2.PrintToScreen()
flushed_print("'\n")

# 3) The two can be aligned...

expPair = ExpSeqPair(s1,s2)
ali = SimpleAligner(expPair)
ali.Align()
expPair.PrintToScreen()
