/**
 * \file MergeWF_Paddles.h
 *
 * \ingroup OpMichel
 * 
 * \brief Class def header for a class OverlayWF_Paddles
 *
 * @author david caratelli
 */

/** \addtogroup OpMichel

    @{*/

#ifndef LARLITE_OVERLAYWF_PADDLES_H
#define LARLITE_OVERLAYWF_PADDLES_H

#include "Analysis/ana_base.h"
#include "SignalProcessing.h"

namespace larlite {
  /**
     \class OverlayWF_Paddles
     User custom analysis class made by SHELL_USER_NAME
   */
  class OverlayWF_Paddles : public ana_base{
  
  public:

    /// Default constructor
    OverlayWF_Paddles();

    /// Default destructor
    virtual ~OverlayWF_Paddles(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    // set producers here
    void setPMTProducer(std::string s)  { _PMTproducer  = s; }
    void setTrigProducer(std::string s) { _TRIGproducer = s; }
    // use MC info
    void useMC(bool on) { _useMC = on; }
    // PE threshold to claim a muon pulse
    void setMuonPEThresh(double pe) { _muon_PE_thresh = pe; }
    // set the baseline PE (what to expect at any given moment)
    void setNoisePE(double pe) { _noise_PE = pe; }
    // set whether to require a muon peak
    void setRequireMuonPeak(bool on) { _require_muon_peak = on; }
    // set maximum time for muon peak
    void setMaximumMuonTime(double t) { _max_muon_time = t; }
    // set maximum number of muons in event's beam-gate RO window
    void setMaximumMuonNumber(size_t n) { _max_muon_number = n; }
    // set the late-light amplitude factor
    void setLateLightAmplitude(double a) { _ll_amp = a; }
    // set the late-light time-constant
    void setLateLightTimeConstant(double t) { _ll_tau = t; }
    // set verobisty flag
    void setVerbose(bool on) { _verbose = on; }
    // set the signal processing module instance
    void setSignalProcessor(signalana::SignalProcessing s) { _signalProcessor = s; }

    // get the signifiance of the value o coming from a 
    // poisson distribution with expectation e
    // measured in sigma
    double LateLightProb(const double& e, const double& o);

  protected:

    void fillTree();

    void cleanTree();

    // verbosity boolean
    bool _verbose;

    // require muon peak?
    bool _require_muon_peak;
    // maximum muon time (since stitched wf starts) [ us ]
    double _max_muon_time;
    // maximum number of muons allowed in event's beam-gate RO window
    size_t _max_muon_number;

    double _muon_PE_thresh;

    // resize waveforms for all 32 PMTs
    void resizeWaveforms(const size_t& nticks);


    // exponential function with which to model the late-light
    double lateLight(const double& time,
		     const double& A,
		     const double& tau);
    // late-light amplitude factor
    double _ll_amp;
    // late-light time-constant
    double _ll_tau;

    // instance of the SignalProcessing class
    signalana::SignalProcessing _signalProcessor;


    // function that returns the poisson expectation for
    // the value k for a poisson distribution of mean lambda
    double Poisson(const int& lambda, const int& k);

    // Poisson CDF: prob of getting a value < k
    double PoissonCDF(const int& lambda, const int& k);

    // approximation for the inverse error function
    double InvERF(const double& x);

    // given a time-tick search the maximum ADC in a small region
    // around this tick for all PMTs
    // and calcuulate a "flash position" by making
    // a weighted average of the PEs seen by each PMT
    std::pair<double,double> getPulseYZ(const size_t& tick);


    // baseline PE to be used in developing PE expectation 
    double _noise_PE;

    // producer names
    std::string _PMTproducer;
    std::string _TRIGproducer;

    // use MC info
    bool _useMC;

    TTree* _tree;
    int _event;
    std::vector<double> _merged_wf;
    std::vector<double> _theory_wf;
    std::vector<double> _hit_significance;
    std::vector<double> _hit_time;
    std::vector<double> _all_muon_time;
    std::vector<double> _all_muon_ampl;
    std::vector<double> _all_muon_Y;
    std::vector<double> _all_muon_Z;
    std::vector<double> _all_hit_time;
    std::vector<double> _all_hit_ampl;
    std::vector<int>    _all_hit_pmtn;
    std::vector<int>    _sig_hit_idx; // index in _all_hit_*
    std::vector<double> _sig_hit_Y;   // Y position of this hit
    std::vector<double> _sig_hit_Z;   // Z position of this hit
    // mc michel info
    int _michel;
    double _michel_E;
    double _michel_t;
    double _michel_x;
    double _michel_y;
    double _michel_z;

    // per-PMT waveforms
    std::vector< std::vector<double> > _pmt_wfs;

    // per-PMT positions (YZ)
    std::vector< std::pair<double,double> > _pmt_pos;

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
