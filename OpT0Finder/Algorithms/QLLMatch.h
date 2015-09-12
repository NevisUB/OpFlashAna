/**
 * \file QLLMatch.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class QLLMatch
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithms

    @{*/
#ifndef QLLMATCH_H
#define QLLMATCH_H

#include <iostream>
#include "Base/BaseFlashMatch.h"
#include <TMinuit.h>
namespace flashana {
  /**
     \class QLLMatch
     User defined class QLLMatch ... these comments are used to generate
     doxygen documentation!
  */
  class QLLMatch : public BaseFlashMatch {
    
  public:
    
    /// Default constructor
    QLLMatch();
    
    /// Default destructor
    ~QLLMatch(){}

    static QLLMatch& GetME() {

      if(!_me) _me = new QLLMatch;
      return *_me;
    }

    FlashMatch_t Match(const QCluster_t&, const Flash_t&);

    void SetOpDetPositions( const std::vector<double>&,
			    const std::vector<double>&,
			    const std::vector<double>&);

    double ChargeHypothesis(double x);
  private:

    void CallMinuit();

    static QLLMatch* _me;
    
    std::vector<double> _pmt_x_v;
    std::vector<double> _pmt_y_v;
    std::vector<double> _pmt_z_v;

    std::vector<double> _trk_x_v;
    std::vector<double> _trk_y_v;
    std::vector<double> _trk_z_v;
    std::vector<double> _trk_q_v;

    std::vector<double> _qll_hypothesis_v;
    std::vector<double> _flash_pe_v;
    
    double _reco_x;
    double _reco_x_err;
    double _qll;

    TMinuit* _minuit_ptr;
  };
}

#endif
/** @} */ // end of doxygen group 

