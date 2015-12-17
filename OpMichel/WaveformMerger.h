/**
 * \file WaveformMerger.h
 *
 * \ingroup OpMichel
 * 
 * \brief Class def header for a class WaveformMerger
 *
 * @author david caratelli
 */

/** \addtogroup OpMichel

    @{*/

#ifndef LARLITE_WAVEFORMMERGER_H
#define LARLITE_WAVEFORMMERGER_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class WaveformMerger
     User custom analysis class made by SHELL_USER_NAME
   */
  class WaveformMerger : public ana_base{
  
  public:

    /// Default constructor
    WaveformMerger();

    /// Default destructor
    virtual ~WaveformMerger(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    // set producers here
    void setPMTProducer(std::string s)  { _PMTproducer  = s; }
    void setTrigProducer(std::string s) { _TRIGproducer = s; }
    // set start time for the output waveform [us]
    void setStartTime(double t) { _t_min = t; }
    // set end time for the output waveform [us]
    void setEndTime(double t) { _t_max = t; }
    // set the baseline value to assign to empty tick [ADCs]
    void setBaseline(short b) { _baseline = b; }
    // set the verbosity
    void setVerbose(bool on) { _verbose = on; }

  protected:

    // verbosity flag
    bool _verbose;

    // range within which to search for waveforms to be merged
    double _t_min;
    double _t_max;

    // baseline value to assign empty ticks?
    short _baseline;

    // sampling [ us ]
    double _sampling;

    // producer names
    std::string _PMTproducer;
    std::string _TRIGproducer;
    
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
