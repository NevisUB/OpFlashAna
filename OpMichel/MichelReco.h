/**
 * \file MichelReco.h
 *
 * \ingroup OpMichel
 * 
 * \brief Class def header for a class MichelReco
 *
 * @author david caratelli
 */

/** \addtogroup OpMichel

    @{*/

#ifndef LARLITE_MICHELRECO_H
#define LARLITE_MICHELRECO_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class MichelReco
     User custom analysis class made by SHELL_USER_NAME
   */
  class MichelReco : public ana_base{
  
  public:

    /// Default constructor
    MichelReco() :
      _michel_tree(nullptr)
    { _name="MichelReco"; _fout=0; _PEmin = 0;}
    
    /// Default destructor
    virtual ~MichelReco(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    /**
     * @brief Set minimum number of PEs for flash to be considered
     */
    void SetPEMin(double pe) { _PEmin = pe; }

  protected:

    // minimum PEs for a flash to be considered
    double _PEmin;

    // map linking time -> flash index
    std::map<double,size_t> _flashMap;

    // tree
    TTree* _michel_tree;
    int _evt;      /// event number
    double _dt;    /// separation between consecutive flashes
    double _T;     /// time of candidate michel
    double _Tpre;  /// time of candidate muon
    double _PEpre; /// PEs for previous flash
    double _PE;    /// PEs for this (possible michel) flash
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
