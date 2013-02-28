// PD is a free, modular C++ library for biomolecular simulation with a 
// flexible and scriptable Python interface. 
// Copyright (C) 2003-2013 Mike Tyka and Jon Rea
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __MOVE_H
#define __MOVE_H

#include "object.h"
#include "manipulators/movebase.h"
#include "manipulators/rotamer_applicatorbase.h"
#include "workspace/workspace.fwd.h"

namespace Library
{
    class PD_API AngleSet;
}

namespace Physics
{
    class PD_API Forcefield;
}

namespace Manipulator
{
    //////////////////////////////////////////////////////////////////////
    //
    //    Cartesian moves
    //

    //-------------------------------------------------
    //
    /// \brief   Moves single atoms by a normal displacement
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API CartesianMove: public MoveBase
    {
    public:
        CartesianMove(
            WorkSpace& _wspace,
            double _fprob,
            int    _ndisplace, 
            double _xsigma)
            : MoveBase(_wspace) 
        {
            fprob = _fprob;
            xsigma = _xsigma;
            ndisplace = _ndisplace;
            RepeatInterval = 1;
            m_MoveCount = 0;
        }

        virtual CartesianMove* clone() const { return new CartesianMove(*this); }

        int apply();

        int RepeatInterval;
    protected:
        int m_MoveCount;
        void displaceAtom(size_t i);
        double fprob; // probability based moves
        int ndisplace; // a fixed number of moves

        double xsigma;
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API CartesianBlockMove: public MoveBase{
    public:
        CartesianBlockMove(
            WorkSpace& _wspace,
            double _fglobal, 
            double _xsigma, 
            double _anglesigma)
            :MoveBase(_wspace) 
        {
            fglobal = _fglobal;
            xsigma = _xsigma;
            anglesigma = _anglesigma;
        }

        virtual CartesianBlockMove* clone() const { return new CartesianBlockMove(*this); }

        int apply();

    protected:
        double fglobal;
        double xsigma;
        double anglesigma;
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API MoleculeDisplacement: public MoveBase
    {
    public:
        MoleculeDisplacement(
            WorkSpace& _wspace,
            double _fprob,
            int    _ndisplace, 
            double _xsigma, 
            double _anglesigma)
            : MoveBase(_wspace) 
        {
            fprob = _fprob;
            ndisplace = _ndisplace;
            xsigma = _xsigma;
            anglesigma = _anglesigma;
        }

        virtual MoleculeDisplacement* clone() const { return new MoleculeDisplacement(*this); }

        int apply();

    protected:
        void   displaceMolecule(size_t im);
        double fprob;                        // probability based moves
        int    ndisplace;                    // a fixed number of moves

        double xsigma;
        double anglesigma;
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API TIP3P_Move: public MoveBase
    {
    public:
        TIP3P_Move(
            WorkSpace& _wspace,
            double _fprob,
            int    _ndisplace, 
            double _trans_sigma, 
            double _rot_sigma,
            double _bond_sigma,
            double _angle_sigma)
            :MoveBase(_wspace) 
        {
            fprob = _fprob;
            ndisplace = _ndisplace;

            trans_sigma = _trans_sigma; 
            rot_sigma   = _rot_sigma;
            bond_sigma  = _bond_sigma;
            angle_sigma = _angle_sigma;
            RepeatMoveInterval = 1;
            movecount = 0;
        }

        virtual TIP3P_Move* clone() const { return new TIP3P_Move(*this); }

        int apply();
        int RepeatMoveInterval;
    protected:
        void   applyMoveToMolecule(size_t im);
        double fprob;                        // probability based moves
        int    ndisplace;                    // a fixed number of moves

        double trans_sigma; 
        double rot_sigma;
        double bond_sigma;
        double angle_sigma;

        std::vector <int> movedmolecules;
        int movecount;
    };


    //////////////////////////////////////////////////////////////////////
    //
    //    Torsional moves
    //


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API NormalTorsionalMove : public MoveBase
    {
    public:
        NormalTorsionalMove(WorkSpace& _wspace):MoveBase(_wspace) 
        {
        }

        virtual NormalTorsionalMove* clone() const { return new NormalTorsionalMove(*this); }

        int apply();
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API SidechainTorsionalMove: public MoveBase
    {
    public:
        SidechainTorsionalMove(
            WorkSpace& _wspace,
            double _p120, 
            double _pnormal, 
            double _sdnormal) 
            : MoveBase(_wspace) 
        {
            p120 = _p120;
            pnormal = _pnormal;
            sdnormal = _sdnormal;
        }

        virtual SidechainTorsionalMove* clone() const { return new SidechainTorsionalMove(*this); }

        int apply();

    protected:
        double p120;
        double pnormal;
        double sdnormal;
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API BackboneTorsionalMove: public MoveBase
    {
    public:
        BackboneTorsionalMove(
            WorkSpace& _wspace,
            double _pnormal, 
            double _sdnormal, 
            int _maxresidues)
            :MoveBase(_wspace) 
        {
            pnormal = _pnormal;
            sdnormal = _sdnormal;
            maxresidues = _maxresidues;
        }

        virtual BackboneTorsionalMove* clone() const { return new BackboneTorsionalMove(*this); }

        int apply();

    protected:
        double pnormal;
        double sdnormal;
        int maxresidues;
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API BackbonePropagationMove: public MoveBase{
    public:
        BackbonePropagationMove(
            WorkSpace& _wspace, 
            double _pprop)
            :MoveBase(_wspace), 
            pprop(_pprop)
        {
            THROW( NotImplementedException, "This code is not tested!");
        }

        virtual BackbonePropagationMove* clone() const { return new BackbonePropagationMove(*this); }

        int apply();

    protected:
        double pprop;
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API PeptideGroupMove: public MoveBase
    {
    public:
        PeptideGroupMove(
            WorkSpace& _wspace,
            double _fglobal, 
            double _prop, 
            double _anglesigma)
            :MoveBase(_wspace),
            fglobal(fglobal),
            prop(_prop),
            anglesigma(_anglesigma)
        {

        }

        virtual PeptideGroupMove* clone() const { return new PeptideGroupMove(*this); }

        int apply();

    protected:
        double fglobal;
        double prop;
        double anglesigma;
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API BackbonePhiPsiSetMove: public MoveBase
    {	

        struct PhiPsiPair
        {
            double phi, psi;
            double prob;
        };

    public:
        BackbonePhiPsiSetMove(
            WorkSpace& _wspace,
            const Library::RotamerLibrary& _rotlib,
            Library::AngleSet * _angset,
            double _fchange,
            int _simultaneous,
            int _local,
            int _sidechainfix)
            : MoveBase(_wspace) 
        {
            rotlib = &_rotlib;
            angset = _angset;
            fchange = _fchange;
            simultaneous = _simultaneous;
            local = _local;
            sidechainfix = _sidechainfix;
            endlength = 4;
        }

        virtual ~BackbonePhiPsiSetMove(){}

        virtual BackbonePhiPsiSetMove* clone() const { return new BackbonePhiPsiSetMove(*this); }

        int apply();
        void changeResidueRandomly(int ir,bool weighted = true);
        void changeConsResiduesRandomly(int sir, int nir);

    private:

        const Library::RotamerLibrary *rotlib;
        Library::AngleSet *angset;
        double fchange;
        int local;
        int sidechainfix;

        /// defines the number of residues from either end of the protein
        /// not to be considered a leaver
        int endlength;
        int simultaneous;
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API BlockRotationMove: public MoveBase
    {
    public:
        BlockRotationMove(
            WorkSpace& _wspace,
            Physics::Forcefield& _stericff,
            double _fchange);


        virtual ~BlockRotationMove();

        virtual BlockRotationMove* clone() const { return new BlockRotationMove(*this); }

        int apply();
        int rotateBlock(int ir1, int ir2, double angle);

    private:
        double fchange;
        Physics::Forcefield *stericff;
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author Mike Tyka & Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API SidechainRotamerLibMove : public RotamerApplicatorBase
    {
    public:
        SidechainRotamerLibMove(
            WorkSpace& _wspace,
            const Library::RotamerLibrary& _rotlib,
            double _pprop,
            double _forcedchangecutoff, 
            Manipulator::RotamerMode _mode = Manipulator::ApplyCartesian );

        virtual SidechainRotamerLibMove* clone() const { return new SidechainRotamerLibMove(*this); }

        int apply();
        double pprop;
        int getRotSterics(int ir, double *rotsterics, bool tryallrotamers, double scorecutoff);
        double forcedchangecutoff;

    private:
        const Library::RotamerLibrary *rotlib;
    };
}

#endif

