/**
 * \file MichelReco_OpDigit.h
 *
 * \ingroup OpMichel
 * 
 * \brief Class def header for a class MichelReco_OpDigit
 *
 * @author david
 */

/** \addtogroup OpMichel

    @{*/

#ifndef LARLITE_MICHELRECO_OPDIGIT_H
#define LARLITE_MICHELRECO_OPDIGIT_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class MichelReco_OpDigit
     User custom analysis class made by SHELL_USER_NAME
   */
  class MichelReco_OpDigit : public ana_base{
  
  public:

    /// Default constructor
    MichelReco_OpDigit();

    /// Default destructor
    virtual ~MichelReco_OpDigit(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setUseMC(bool on) { _use_mc = on; }

  protected:

    void fillTree();

    void cleanTree();

    void add(const short& pmt,
	     const std::vector<short>& ADCmax,
	     const std::vector<short>& TDCmax,
	     const std::vector<short>& wf);

    int findMaxima(const std::vector<short>& wf,
		   const size_t& tick_min,
		   const size_t& tick_max,
		   std::vector<short>& ADCmax,
		   std::vector<short>& TDCmax);

    void findMuonPeak(const std::vector<short>& wf,
		      size_t& max_tdc,
		      short& max_adc);

    short _baseline;

    bool _use_mc;

    TTree* _tree;
    int _event;
    int _michel;
    double _E;
    double _t;
    double _x;
    double _y;
    double _z;
    std::vector<short> _tdc00;
    std::vector<short> _tdc01;
    std::vector<short> _tdc02;
    std::vector<short> _tdc03;
    std::vector<short> _tdc04;
    std::vector<short> _tdc05;
    std::vector<short> _tdc06;
    std::vector<short> _tdc07;
    std::vector<short> _tdc08;
    std::vector<short> _tdc09;
    std::vector<short> _tdc10;
    std::vector<short> _tdc11;
    std::vector<short> _tdc12;
    std::vector<short> _tdc13;
    std::vector<short> _tdc14;
    std::vector<short> _tdc15;
    std::vector<short> _tdc16;
    std::vector<short> _tdc17;
    std::vector<short> _tdc18;
    std::vector<short> _tdc19;
    std::vector<short> _tdc20;
    std::vector<short> _tdc21;
    std::vector<short> _tdc22;
    std::vector<short> _tdc23;
    std::vector<short> _tdc24;
    std::vector<short> _tdc25;
    std::vector<short> _tdc26;
    std::vector<short> _tdc27;
    std::vector<short> _tdc28;
    std::vector<short> _tdc29;
    std::vector<short> _tdc30;
    std::vector<short> _tdc31;
    std::vector<short> _adc00;
    std::vector<short> _adc01;
    std::vector<short> _adc02;
    std::vector<short> _adc03;
    std::vector<short> _adc04;
    std::vector<short> _adc05;
    std::vector<short> _adc06;
    std::vector<short> _adc07;
    std::vector<short> _adc08;
    std::vector<short> _adc09;
    std::vector<short> _adc10;
    std::vector<short> _adc11;
    std::vector<short> _adc12;
    std::vector<short> _adc13;
    std::vector<short> _adc14;
    std::vector<short> _adc15;
    std::vector<short> _adc16;
    std::vector<short> _adc17;
    std::vector<short> _adc18;
    std::vector<short> _adc19;
    std::vector<short> _adc20;
    std::vector<short> _adc21;
    std::vector<short> _adc22;
    std::vector<short> _adc23;
    std::vector<short> _adc24;
    std::vector<short> _adc25;
    std::vector<short> _adc26;
    std::vector<short> _adc27;
    std::vector<short> _adc28;
    std::vector<short> _adc29;
    std::vector<short> _adc30;
    std::vector<short> _adc31;
    std::vector<short> _wf00;
    std::vector<short> _wf01;
    std::vector<short> _wf02;
    std::vector<short> _wf03;
    std::vector<short> _wf04;
    std::vector<short> _wf05;
    std::vector<short> _wf06;
    std::vector<short> _wf07;
    std::vector<short> _wf08;
    std::vector<short> _wf09;
    std::vector<short> _wf10;
    std::vector<short> _wf11;
    std::vector<short> _wf12;
    std::vector<short> _wf13;
    std::vector<short> _wf14;
    std::vector<short> _wf15;
    std::vector<short> _wf16;
    std::vector<short> _wf17;
    std::vector<short> _wf18;
    std::vector<short> _wf19;
    std::vector<short> _wf20;
    std::vector<short> _wf21;
    std::vector<short> _wf22;
    std::vector<short> _wf23;
    std::vector<short> _wf24;
    std::vector<short> _wf25;
    std::vector<short> _wf26;
    std::vector<short> _wf27;
    std::vector<short> _wf28;
    std::vector<short> _wf29;
    std::vector<short> _wf30;
    std::vector<short> _wf31;
    
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
