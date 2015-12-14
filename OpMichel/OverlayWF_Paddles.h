/**
 * \file MergeWF_Paddles.h
 *
 * \ingroup OpMichel
 * 
 * \brief Class def header for a class OverlayWF_Paddles
 *
 * @author david
 */

/** \addtogroup OpMichel

    @{*/

#ifndef LARLITE_OVERLAYWF_PADDLES_H
#define LARLITE_OVERLAYWF_PADDLES_H

#include "Analysis/ana_base.h"

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
    // PE threshold to claim a hit
    void setHitPEDifferentialThresh(double pe) { _hit_PE_differential_thresh = pe; }
    // set dead-time for muons
    void setDeadTime(double t) { _deadTime = t; }
    // set the baseline PE (what to expect at any given moment)
    void setBaselinePE(double pe) { _baseline_PE = pe; }
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

    // resize waveforms for all 32 PMTs
    void resizeWaveforms(const size_t& nticks);

    // find the time of a pulse based on the differential
    // of the waveform
    // search for peaks followed by zero-crossings
    // the zero-crossing is the time of the pulse
    void findPulseTimes(const std::vector<double>& wfdiff,
			std::vector<size_t>& pulseTicks);
    // PE threshold for hit amplitude
    double _hit_PE_differential_thresh;

    
    // for each hit tick time, search in a neightborhood
    // of that tick for the local maximum
    // that will be the corrected hit time
    // how far out to search depends on shaping time
    void findLocalMaxima(const std::vector<double>& wf,
			 std::vector<size_t>& hitTicks);

    // given a vector of values and time-ticks
    // return a vector of times (in us)
    // and a vector of values for those points
    void getPoints(const std::vector<double>& wf,
		   const std::vector<size_t>& pts,
		   std::vector<double>& times,
		   std::vector<double>& vals);

    // given the waveform and the tick-values at which
    // the pulses are to be found return a vector of 
    // peak-times [usec] and ADCs [baseline subtracted]
    // for each point, scan neighbors to find the true maximum
    void getPeaks(const std::vector<double>& wf,
		  const std::vector<size_t>& ticks,
		  std::vector<double>& times,
		  std::vector<double>& PEs);

    // exponential function with which to model the late-light
    double lateLight(const double& time,
		     const double& A,
		     const double& tau);
    // late-light amplitude factor
    double _ll_amp;
    // late-light time-constant
    double _ll_tau;

    
    // function to get rms and baseline
    void getRMS(const std::vector<double>& wf,
		double& base, double& rms);


    // function to get the rolling baseline from a vector
    void getBaseline(const std::vector<short>& wf,
		     std::vector<double>& baseline);


    // get the differential waveform
    void getDifferential(const std::vector<double>& wf,
			 std::vector<double>& wfdiff);

    // function to find all muon peaks in a waveform
    void getMuonPeaksTicks(const std::vector<double>& wf,
			   std::vector<size_t>& muonTicks);
    // threshold for muon in PE
    double _muon_PE_thresh;

    // function that applies a dead-time for the muon peaks
    void applyMuonDeadTime(std::vector<double>& peakTimes,
			   std::vector<double>& peakAmps);
    // dead-time for muons
    double _deadTime;

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
    double _baseline_PE;

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
