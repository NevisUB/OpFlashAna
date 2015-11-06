/**
 * \file LateLightComponent.h
 *
 * \ingroup FlashStudies
 * 
 * \brief Class def header for a class LateLightComponent
 *
 * @author david
 */

/** \addtogroup FlashStudies

    @{*/

#ifndef LARLITE_LATELIGHTCOMPONENT_H
#define LARLITE_LATELIGHTCOMPONENT_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class LateLightComponent
     User custom analysis class made by SHELL_USER_NAME
   */
  class LateLightComponent : public ana_base{
  
  public:

    /// Default constructor
    LateLightComponent();

    /// Default destructor
    virtual ~LateLightComponent(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setMinPE(double pe) { _min_pe = pe; }

  protected:

    TTree* _tree;
    // tree branches
    int _n_flash;
    double _t_flash;
    double _pe_flash;
    double _q_tot;
    double _q_max; // hit in event with largest amount of charge
    double _t_start;
    double _dt;
    std::vector<double> _hit_t_v;
    std::vector<double> _hit_q_v;

    // min pe
    double _min_pe;
    
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
