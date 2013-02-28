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

#ifndef __PDB_H
#define __PDB_H

#include <string>
#include <vector>
#include <fstream>

#include "library/nameconv.h"    // required as it provides a member variable
#include "library/residues.h"     // required as it provides a member variable
#include "sequence/sequence.h"   // required as it provides a member variable
#include "tools/streamwriter.h"
#include "fileio/outtra.h"       // required as it provides a base class
#include "fileio/infile.h"       // required as it provides a base class
#include "workspace/workspace.h" // required as it provides a base class
#include "system/system.h"       // required as it provides a base class

class PD_API Particle;

namespace Sequence
{
    class PD_API ExpSeqPair;
    class PD_API BioSequence;

#ifdef SWIG
    %template(BioSequenceCollection_char) BioSequenceCollection<char>; 
#endif
}

namespace IO
{
    class File_Molecule;
    class PD_API PDB_ImportBase;

    /// Define a **flagged enumeration** of the different named elements that we can import from a PDB file
    enum HeaderTagFilter
    {
        // Header flags
        HEADER = 1ul,
        TITLE  = 2ul,
        COMPND = 4ul,
        SOURCE = 8ul,
        KEYWDS = 16ul,
        EXPDTA = 32ul,
        AUTHOR = 64ul,
        REVDAT = 128ul,
        JRNL   = 256ul,
        REMARK = 512ul,
        FORMUL = 1024ul,
        HELIX  = 2048ul,
        SHEET  = 4096ul,
        HET    = 8192ul,
        DBREF  = 16384ul,

        // System Information Flags
        SEQRES = 32768ul,
        MODRES = 65536ul,
        SSBOND = 131072ul,
        CONECT = 262144ul,

        // Unknown ???
        UNKNOWN_TAGS = 524288ul,

        // Now define sets of these basic elements that can be imported
        NoHeader = 0, // Just give me the positions bitch ;-)
        SystemInfo = SEQRES | MODRES | SSBOND | CONECT,
        Minimal = SystemInfo | HEADER | TITLE | KEYWDS | EXPDTA,
        Comprehensive = Minimal | AUTHOR | REVDAT | JRNL | HELIX | SHEET | DBREF,
        AllKnown = Comprehensive | COMPND | SOURCE | REMARK | FORMUL | HET,
        Everything = AllKnown | UNKNOWN_TAGS,

        DefaultHeader = Minimal  // What you would normally want to import
    };


    // These are the current experimental methods you would expect to find in a PDB file.
    enum PDB_ExpMethod
    {
        Crystalographic,
        InfraredMicroscopy,
        FiberDiffraction,
        ElectronMicroscopy,
        SingleCrystalElectronMicroscopy,
        NMR,
        UnknownMethod,
        Undefined
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API PDB_SEQRES: public Sequence::BioSequenceCollection<char>
    {
    public:
        PDB_SEQRES();

        void setAlias(const Library::AliasMapper& _ParserNameMapper);
        const Library::AliasMapper& getAlias() const;

        void clear();
        void loadFromFile( const std::string &_fileName );
        virtual void parseLine( const std::string & _line);
        const Library::AliasMapper* getNameMapper() const;

    private:
        const Library::AliasMapper* m_ParserNameMapper;

        // Per instance variable holders - saving the previous state of the parseLine() function call
        void resetStatics(); ///< Reset the member state-variables to their original states for parsing
    };


    //-------------------------------------------------
    //
    /// \brief Hold a PDB ATOM line 
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// Official PDB File Format v2.2:
    /// http://www.pdb.org/pdb/file_formats/pdb/pdbguide2.2/guide2.2_frame.html
    /// PDB File Data Column Descriptors:
    /// lineType (1-6): either 'ATOM ' or 'HETATM'
    /// Atom Number (7-11)
    /// Atom Name (13-16)
    /// AltLoc Indicator (17)
    /// Residue Name (18-20)
    /// Chain ID (22)
    /// Residue Number (23-26)
    /// Insertion Code (27)
    /// X Coordinate (31-38)
    /// Y Coordinate (39-46)
    /// Z Coordinate (47-54)
    /// Occupancy (55-60)
    /// Temperature Factor (61-66)
    /// Segment Identifier (73-76)
    /// Element symbol (77-78)
    /// Atomic Charge (79-80)
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API PDBAtomLine : public FileParticle // struct representing a single PDB line
    {
    public:
        enum AtomLineType
        {
            ATOM,
            HETATM,
            TER // Special terminal record - signifies the end of a polymer
        };

        // additional PDB-Specific Members not present in base class 'FileParticle'
        AtomLineType lineType; 
        char altLoc; 
        double occupancy;    
        double TempFactor;   
        std::string segID;   
        std::string element;
        std::string charge;  

        PDBAtomLine();
        void init();

        bool loadAtom(const Particle &atom, int atomIndex);
        bool loadLine(const std::string& _Line);
        void saveLine(std::ostream &_file) const;
    };


    //-------------------------------------------------
    //
    /// \brief  The "MODRES" pdb file line Type
    ///
    /// \details 
    /// Defines the 'normal residue Type' associated with rather bizare names that are given,
    /// to residues with small alterations.
    /// 
    /// COLUMNS        DATA TYPE       FIELD          DEFINITION
    /// --------------------------------------------------------------------------------
    ///  1 -  6        Record name     "MODRES"
    ///  8 - 11        IDcode          idCode         ID code of this entry.
    /// 13 - 15        Residue name    resName        Residue name used in this entry.
    /// 17             Character       chainID        Chain identifier.
    /// 19 - 22        Integer         seqNum         Sequence number.
    /// 23             AChar           iCode          Insertion code.
    /// 25 - 27        Residue name    stdRes         Standard residue name.
    /// 30 - 70        String          comment        Description of the residue modification.
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class ModResLine : public Library::ResidueAliasDefinition
    {
        friend class PDB_MODRES;
    public:
        ModResLine();

        bool parseLine( const std::string& line );
        bool matches(const FileParticle& particle) const;
        bool matches(const std::string& resName) const;

    protected:
        char chainID;
        int resNum;
        char iCode;
        std::string chemicalName;
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PDB_MODRES : public Library::AliasMapper
    {
    public:
        PDB_MODRES();

        void parseLine( const std::string& line );
        size_t getCount() const { return m_Lines.size(); }

        /// 'AttemptRename()' is not currently used, but is much more specific than '_LookupAlias_',
        /// being as it will only rename if the residue index and chain also match.
        /// The 'problem' is that '_LookupAlias_' has no provision/need for residue-position specicic name defintions -
        /// the base class assumes than a given alias has a constant mapping to a given residue Type -
        /// This never changes, and I do not believe it should. This function will however be kept as it is the
        /// official PDB definition of a MODRES line.
        bool AttemptRename( FileParticle& particle ) const; 

        virtual bool lookupLongName( char _SourceID, std::string& _DestID ) const; ///< Look up the long name for a given single letter name
        virtual bool lookupShortName( const std::string& _SourceID, char& _DestID ) const; ///< Look up the single letter for a long residue name
        virtual bool lookupAlias( const std::string& _SourceID, std::string& _DestID ) const; ///< Look to see if the source name is an alias, and return the true name

    private:
        std::vector<ModResLine> m_Lines;
    };



    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    /// SSBOND Definition:
    /// 
    /// COLUMNS       DATA TYPE       FIELD          DEFINITION
    /// ----------------------------------------------------------------------------
    ///  1 -  6       Record name     "SSBOND"
    ///  8 - 10       Integer         serNum         Serial number.
    /// 12 - 14       LString(3)      "CYS"          Residue name.
    /// 16            Character       chainID1       Chain identifier.
    /// 18 - 21       Integer         seqNum1        Residue sequence number.
    /// 22            AChar           icode1         Insertion code.
    /// 26 - 28       LString(3)      "CYS"          Residue name.
    /// 30            Character       chainID2       Chain identifier.
    /// 32 - 35       Integer         seqNum2        Residue sequence number.
    /// 36            AChar           icode2         Insertion code.
    /// 60 - 65       SymOP           sym1           Symmetry operator for 1st residue.
    /// 67 - 72       SymOP           sym2           Symmetry operator for 2nd residue.
    ///
    /// CONECT Definition:
    /// 
    /// COLUMNS         DATA TYPE        FIELD           DEFINITION
    /// ---------------------------------------------------------------------------------
    ///  1 -  6         Record name      "CONECT"
    ///  7 - 11         Integer          serial          Atom serial number
    /// 12 - 16         Integer          serial          Serial number of bonded atom
    /// 17 - 21         Integer          serial          Serial number of bonded atom
    /// 22 - 26         Integer          serial          Serial number of bonded atom
    /// 27 - 31         Integer          serial          Serial number of bonded atom
    /// 32 - 36         Integer          serial          Serial number of hydrogen bonded atom
    /// 37 - 41         Integer          serial          Serial number of hydrogen bonded atom
    /// 42 - 46         Integer          serial          Serial number of salt bridged atom
    /// 47 - 51         Integer          serial          Serial number of hydrogen bonded atom
    /// 52 - 56         Integer          serial          Serial number of hydrogen bonded atom
    /// 57 - 61         Integer          serial          Serial number of salt bridged atom
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PDB_Connections : public FileCovalency
    {
    public:
        PDB_Connections();
        void parseCONECT( const std::string& line );
        void parseSSBOND( const std::string& line );
    };


    /// \brief  
    /// Defines a series of PDB atom lines from a single molecule. 
    ///
    /// \details 
    /// Derives from FileMoleculeMaker, meaning that these lines can be converted into
    /// A File_Molecule structure for use by FileImportBase
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PDB_Model : public FileMoleculeMaker
    {
    public:
        friend class PDB_ImportBase;

        PDB_Model(const PDB_ImportBase& _PDBParent, const FileImportBase& _FileParent);
        ~PDB_Model();
        //PDB_Model( const PDB_Model &_clone );
        //const PDB_Model& PDB_Model::operator= (const PDB_Model &_clone);

        PDB_Model( const PDB_Model &_clone ) :
        pdbLines(_clone.pdbLines),
            m_PDBParent(_clone.m_PDBParent),
            FileMoleculeMaker(_clone.getFileParent())
        {
        }
#ifndef SWIG
        const PDB_Model& operator= (const PDB_Model &_clone)
        {
            pdbLines = _clone.pdbLines;
            m_PDBParent = _clone.m_PDBParent;
            FileMoleculeMaker::operator=(_clone);
            return (*this);
        }
#endif
        virtual File_Molecule make( char chainID );
        virtual File_Molecule make( char chainID, const Sequence::BioSequence* buildSequenceOverride ); /// 'buildSequenceOverride' must be a pointer because we need to have a possible unallocated state
        virtual std::vector<char> findChainIDs() const;
        virtual bool HasChainID(char _RequestChain ) const;

        const PDB_ImportBase& getPDBParent() const;

    protected:
        // Line import and adding
        inline void appendLine( PDBAtomLine &_Line ) { pdbLines.push_back(_Line); }
        inline int getLineCount() { return (int) pdbLines.size(); }

        bool useSEQRES() const; ///< Flags whether this class will utilise SEQRES to to make FileMolecules

        std::vector<PDBAtomLine> pdbLines; // The lines from the PDB file corresponding to this model
        const PDB_ImportBase* m_PDBParent;
    };


    //-------------------------------------------------
    //
    /// \brief  Provides base functionality to all PDB importing classes (i.e. PDB_In and PDB_Tools)
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API PDB_ImportBase
    {
    public:
        // Constructor and Destructor. Library::AliasMapper is required for residue name assignment for sequence alignment.
        PDB_ImportBase();

        // Public const member accessors
        const PDB_SEQRES& getSEQRES() const { return m_SEQRES; }
        const PDB_MODRES& getMODRES() const { return m_MODRES; }
        const PDB_Connections& getConnections() const { return m_Connections; }

        // Public variables
        int MaxModelImport; ///< The maximum number of different protein models to import. 0 to import none, -ve numbers disable
        HeaderTagFilter HeaderFilterMode; ///< Affects the initial scan of the file. Has no effect after ScanFile is called
        bool UseSEQRES; ///< Should load() functions use SEQRES definitions as a sequence data source?

    protected:
        // Protected Verbosity::Type settings
        Verbosity::Type getVerbose() const { return m_Verbose; }
        void setVerbosity( Verbosity::Type _BeVerbose ) { m_Verbose = _BeVerbose; }

        // Protected member accessors
        PDB_SEQRES& getSEQRES() { return m_SEQRES; }
        PDB_MODRES& getMODRES() { return m_MODRES; }
        PDB_Connections& getConnections() { return m_Connections; }

        // Filename and import options
        void ScanFile( const FileImportBase& _DerivedImporter ); /// fill our arrays from the PDB file lines

        // Sequence Models and Header
        std::vector< std::string > m_PDBHeader; ///< Stores lines from the PDB file header as defined by 'HeaderFilterMode'
        std::vector< PDB_Model > m_Models; ///< Stores models imported from the PDB file.

        // Biological Sequence Container
        PDB_SEQRES m_SEQRES; ///< Class representing the details held in the PDB files SEQRES statements. These effectively hold the biological sequenec of the protein.
        PDB_MODRES m_MODRES; ///< Use MODRES definitions to use standard residue types when the complex Type is not defined by the forcefield
        PDB_Connections m_Connections; ///< Provision for additional CYX-CYX bonds. i.e. SSBONDS

        // File Attributes
        float m_Resolution; ///< The resolution of the experimental data in the file. 
        PDB_ExpMethod m_ExpMethod;
        std::string m_NMRMethod;
        int m_UnrecognisedLines; ///< Counter for any Type of line "not-yet-coded-for" (and either used or explicitly ignored)

        // Helper functions
        void parseExpData(const std::string &line);
        void printExpdata() const;
        void printExpdata( std::ofstream &_file ) const;

    private:
        Verbosity::Type m_Verbose;
    };


    //-------------------------------------------------
    //
    /// \brief  PDB Babel does name-set conversions on whole files
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API PDB_Babel : public FileBabelBase
    {
    public:
        PDB_Babel();
        ~PDB_Babel();
        void RenameFile( const std::string& _InFilename, const std::string& _OutFilename );
    };


    //-------------------------------------------------
    //
    /// \brief  
    /// PDB Tools is analagous to Tra_Tools, it should be invoked for geometry calculations on a
    /// PDB file which do not require the use of a named forcefield.
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API PDB_Tools: public FileImportBase, public PDB_ImportBase
    {
    public:
        PDB_Tools( const std::string& _filename );
        void ensureScan();
    };


    //-------------------------------------------------
    //
    /// \brief  
    /// PDB_In is analagous to Tra_In. It derives from System and therefore can be used to
    /// import a system into a sysspec.
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    /// On Importing polymer sub-ranges...
    ///
    /// void loadFrom( unsigned int modelIndex, char chainID, int startIndex, int length, ResidueClass _Filter = All );
    /// void loadFrom1stModel( char chainID, int startIndex, int length, ResidueClass _Filter = All );
    ///
    /// These is not currently implemented.
    /// Do we use sequential start and length - or PDB start and length??
    /// It is probably better and more generalised to have a
    /// Molecule::makeSubMolecule(int startIndex, int length) instead
    /// This should be facilitated by calls made to the parent class instance of 'System'
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API PDB_In: public FileInBase, public PDB_ImportBase
    {
    public:
        PDB_In(const FFParamSet &_ffps, const std::string& _filename);

        void setVerbosity( Verbosity::Type _BeVerbose ) 
        { 
            FileInBase::setVerbosity(_BeVerbose); 
            PDB_ImportBase::setVerbosity(_BeVerbose); 
        }

        void ensureScanned(); ///< Calls PDB_ImportBase::scan() if it hasn't already been called. This is automatically called, but a public call won't hurt in any way.

        // Full Load Call
        void loadAll(); ///< Load everything from all models in the context of the current filter (See SetFilter())

        // Normal load calls
        void loadModel( unsigned int modelIndex ); ///< Load everything a given model in the context of the current filter (See SetFilter())
        void load(); ///< Load everything from the first model in the context of the current filter (See SetFilter())
        void load( char chainID ); ///< Load everything a given chain in the context of the current filter (See SetFilter())
        void load( char chainID, unsigned int modelIndex ); ///< Load everything a given chain and model in the context of the current filter (See SetFilter())

        // Here we use the PDB file as a structural model from which we allocate all particles we can to an arbritary 
        // sequence via an internal aligment of the data.
        void loadExplicit( char chainID, Library::ResidueClass _class, unsigned int modelIndex = 0);
        void loadExplicit( char chainID, Library::ResidueClass _class, const Sequence::BioSequence& _SequenceOverride, unsigned int modelIndex = 0 );
        void loadExplicit( char chainID, const std::string& resName, int residueNum, char iCode, unsigned int modelIndex = 0 );


    private:
        void loadExplicitInner( char chainID, Library::ResidueClass _class, const Sequence::BioSequence* _SequenceOverride, unsigned int modelIndex);
        bool m_ScanPerformed; ///< Once load is called, the file will be scanned and this will be flagged
    };


    //-------------------------------------------------
    //
    /// \brief  
    /// Class containing the RAW calls to create PDB files of different types
    /// There are occasions when you may want to have access to the raw calls.
    /// For a better end user inteface, use the derived class PDB_Writer.
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PDB_RawWriter : public StreamWriter
    {
    public:
        PDB_RawWriter();

        bool UseBioInfo;

        /// Reset the internal model counter
        void resetModelCount() { m_ModelIndex = 0; }

        // Expose the underlying scream session control
        void beginStream() { StreamWriter::beginStream(); }
        void endStream()   { StreamWriter::endStream(); }

        /// REMARK statements about the origin of the file
        void rawRemarks() const;

        // Calls to start and end a given mode
        void rawBeginModel() const; ///< Prints the MODEL start tag
        void rawEndModel() const; ///< Prints the ENDMDL end tag

        void rawEnd() const; ///< Prints the END tag

        void rawWriteFullHeader( const MoleculeBase& _Mol ) const;
        void rawWriteFullHeader( const System& _Sys ) const;

        // Write data to the underlying stream
        int rawWrite( const MoleculeBase& _Mol, int _AtomOffset = 0 ) const;
        int rawWrite( const System& _Sys, int _AtomOffset = 0 ) const;
        int rawWrite( const std::vector<PDBAtomLine> &_lines, char chainID = ' ', int _AtomOffset = 0 ) const;
        int rawWrite( const std::vector<Maths::dvector> &_Coordinates, char chainID = ' ', int _AtomOffset = 0 )const;

        void rawWrite( 
            const Maths::dvector &_Coordinate, 
            const std::string& atomName, 
            const std::string& resName, 
            int& atomnum, 
            int& resnum 
            ) const;

        void rawWrite( 
            const std::vector<Maths::dvector> &_Coordinates, 
            int& atomnum, 
            int& resnum 
            ) const;

    protected:
        mutable int m_ModelIndex; ///< The current model index to be printed. Needs to be mutable otherwise you cant call rawBeginModel from a const pointer.
    };


    //-------------------------------------------------
    //
    /// \brief  BRIEF DESCRIPTION
    /// PDB writing class with nice user friendly calls to write a 
    /// System, Individal Molecule or Workspace 
    ///
    /// \details DETAILED USER'S DESCRIPTION
    ///
    /// \author  Jon Rea 
    ///
    /// \todo STATE OF DEVELOPMENT
    ///
    /// \bug BUGS?
    ///
    class PD_API PDB_Writer : private PDB_RawWriter
    {
    public:		
        PDB_Writer( const std::string& _filename, bool _MultiModel = false ); ///< Writes to a filename
        PDB_Writer( const char* _filename, bool _MultiModel = false ); ///< Writes to a filename
        PDB_Writer( std::ostream& _stream, bool _MultiModel = false ); ///< Writes to an arbritrary stream (can be std::cout)
        PDB_Writer( bool _MultiModel ); ///< Writes to the screen

        void write( const MoleculeBase& _Mol ); ///< Write our molecule
        void write( const System& _Sys ); ///< Write an entire system

    protected:
        bool m_MultiModel;
        bool m_DoneHeader;
    };


    //-------------------------------------------------
    //
    /// \brief Creates a makeshift trajectory using the PDB format and saving each structure as a seperate model. 
    ///
    /// \details 
    ///
    /// \author  Jon Rea 
    ///
    /// \todo RENAME to OutTra_PDB_MultiModel 
    ///
    /// \bug BUGS?
    ///
    class PD_API OutTra_PDB: public OutputTrajectoryFile, public PDB_RawWriter
    {
    public:
        OutTra_PDB( const std::string &_filestem, WorkSpace &_ps, bool _CallCreate = false );
        virtual OutTra_PDB* clone() const { return new OutTra_PDB(*this); }

        virtual int create();
        virtual int append();
    };


    //-------------------------------------------------
    //
    /// \brief Creates a single-entry PDB file when the trajector(y/ies) are created. This is useful when creating a
    ///        typical 2-file trajectory (e.g. PDB/XTC or PDB/DCD)
    ///
    /// \details 
    ///        It is a direct decendant of OutTra_PDB_MultiModel but has the append function essentailly blanked out.
    ///
    /// \author Mike Tyka 
    ///
    class PD_API OutTra_PDB_Single : public OutTra_PDB 
    {
    public:
        OutTra_PDB_Single( const std::string &_filestem, WorkSpace &_ps ):
          OutTra_PDB( _filestem, _ps, false )
          {};
          virtual OutTra_PDB_Single* clone() const { return new OutTra_PDB_Single(*this); }

          virtual int append();
    };


    //-------------------------------------------------
    //
    /// \brief Saves a WorkSpace, System or Molecule to a PDB file
    ///       
    /// \details This class is used to create PDB files from a WorkSpace, System or Molecule.
    ///          Note that this only creates single files, not trajectories. For Trajectories use
    ///          the class hierarchy deriving from OutputTrajectory.
    ///          FileOutut_PDB can be called in two different ways. Either a system/workspace is passed as an argument
    ///          or the "save" function of a system/workspace is used. Both lead to the same result:
    ///          THe constructor requires a filename or a filestem (the extension is added automatically if its missing). 
                 /*! 
                 \code

                 // examples
                 OutputFile_PDB  mypdb("test.pdb");
                 mypdb.save( mysystem );
                 mypdb.save( myworkspace );
                 mypdb.save( mymolecule );

                 mysystem.save( OutputFile_PDB("test.pdb") );
                 myworkspace.save( OutputFile_PDB("test.pdb") );
                 mymolecule.save( OutputFile_PDB("test.pdb") );

                 \endcode
                 */
    /// \author Mike Tyka 
    ///
    class PD_API OutputFile_PDB : public OutputFile
    {
    public:
        OutputFile_PDB(const std::string& _filestem): OutputFile( _filestem ) {}
        virtual OutputFile* clone() const { return new OutputFile_PDB( *this ); } 

        /// Saves a WorkSpace to a PDB file
        virtual void save( WorkSpace &_wspace );
        virtual void save( System &_system );
    };
} // namespace IO

#endif

