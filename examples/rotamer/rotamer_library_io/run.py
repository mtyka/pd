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

## The code below details possible import options.
## The user only ever wants one of these !!

from pd import *


ffps = FFParamSet()

## Use that Charmm forcefield
ffps.readLib( after("-ff","amber03aa.ff")  )

## **NOTE**
## Readlib imports PD format rotamer libraries
## convertLib imports other format rotamer libraries by taking an additional Conversion class

rot = RotamerLibrary(ffps);

## 1) Shetty format
## Import from a shetty file-pair
rot.convertLib("../../../param/rotlib/shetty/scl-B30-occ1.0-rmsd1.0-prop20.0",RotLibConvert_Shetty());
rot.writeLib("shetty.rotamer"); ## write the lib in PD format to this path

## 2) Legacy PD format
## This will clear the previous import and start again
## Import an old-format PD library
rot.convertLib("../../../param/rotlib/legacy.rotamer",RotLibConvert_OldPDFormat());


## 3) Torsional:Cartesian interconversion! (Cartesian application has benchmarked as faster. Torsional models can be imported,
## and are, by default, applied as cartesian coordinated).
rot.writeLib("output.cart.rotamer"); ## it was imported as a cartesian lib, so this will write cartesian definitions.
rot.writeLib("output.torsion.rotamer", WriteTorsional); ## Force output in torsional format
## automatically clears the current library and imports again.
rot.readLib("output.torsion.rotamer"); ## This time import is of torsional definitions!
## We can now force writing cartesian definitions. They should be the same, bar floating-point rounding
## and anomalies due to bond length and angle deviations between the original definitions and the idealised forms used in torsional conversion.
rot.writeLib("output2.cart.rotamer", WriteCartesian); ## (should match "output.torsion.rotamer")

## 4) Dunbrack format
## Read the Dunbrack backbone independent library
rot.convertLib("../../../param/rotlib/dunbrack/bbind02.May.rot",RotLibConvert_Dunbrack_BBInd());

## Read the Dunbrack backbone independent library - doesnt work 
#dunbrackConv = RotLibConvert_Dunbrack_BBDep();
#rot.convertLib("../../param/rotlib/dunbrack/bbdep02.May.rot",dunbrackConv); ## By-chi mode (best, default)
#dunbrackConv.switchImportMode = true; ## Flag a change to by-rotid mode
#rot.convertLib("../../param/rotlib/dunbrack/bbdep02.May.rot",dunbrackConv); ## By-rotid mode (alternative)
