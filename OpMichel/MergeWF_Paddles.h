/**
 * \file MergeWF_Paddles.h
 *
 * \ingroup OpMichel
 * 
 * \brief Class def header for a class MergeWF_Paddles
 *
 * @author david
 */

/** \addtogroup OpMichel

    @{*/

#ifndef LARLITE_MERGEWF_PADDLES_H
#define LARLITE_MERGEWF_PADDLES_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class MergeWF_Paddles
     User custom analysis class made by SHELL_USER_NAME
   */
  class MergeWF_Paddles : public ana_base{
  
  public:

    /// Default constructor
    MergeWF_Paddles();

    /// Default destructor
    virtual ~MergeWF_Paddles(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    // set producers here
    void setPMTProducer(std::string s)  { _PMTproducer  = s; }
    void setTrigProducer(std::string s) { _TRIGproducer = s; }
    // use MC info
    void useMC(bool on) { _useMC = on; }

  protected:

    void fillTree();

    void cleanTree();

    void add(const short& pmt,
	     const std::vector<short>& wf,
	     const std::vector<double>& th,
	     const std::vector<double>& pmt_wf,
	     const std::vector<double>& wfdiff);

    // find the time of a pulse based on the differential
    // of the waveform
    // search for peaks followed by zero-crossings
    // the zero-crossing is the time of the pulse
    void findPulseTimes(const std::vector<double>& wfdiff,
			std::vector<size_t>& pulseTicks);

    
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


    // function that returns the poisson expectation for
    // the value k for a poisson distribution of mean lambda
    double Poisson(const int& lambda, const int& k);

    // Poisson CDF: prob of getting a value < k
    double PoissonCDF(const int& lambda, const int& k);

    // approximation for the inverse error function
    double InvERF(const double& x);

    // get the signifiance of the value o coming from a 
    // poisson distribution with expectation e
    // measured in sigma
    double LateLightProb(const double& e, const double& o);

    // producer names
    std::string _PMTproducer;
    std::string _TRIGproducer;

    // use MC info
    bool _useMC;

    TTree* _tree;
    int _event;
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
    std::vector<double> _wfdiff00;
    std::vector<double> _wfdiff01;
    std::vector<double> _wfdiff02;
    std::vector<double> _wfdiff03;
    std::vector<double> _wfdiff04;
    std::vector<double> _wfdiff05;
    std::vector<double> _wfdiff06;
    std::vector<double> _wfdiff07;
    std::vector<double> _wfdiff08;
    std::vector<double> _wfdiff09;
    std::vector<double> _wfdiff10;
    std::vector<double> _wfdiff11;
    std::vector<double> _wfdiff12;
    std::vector<double> _wfdiff13;
    std::vector<double> _wfdiff14;
    std::vector<double> _wfdiff15;
    std::vector<double> _wfdiff16;
    std::vector<double> _wfdiff17;
    std::vector<double> _wfdiff18;
    std::vector<double> _wfdiff19;
    std::vector<double> _wfdiff20;
    std::vector<double> _wfdiff21;
    std::vector<double> _wfdiff22;
    std::vector<double> _wfdiff23;
    std::vector<double> _wfdiff24;
    std::vector<double> _wfdiff25;
    std::vector<double> _wfdiff26;
    std::vector<double> _wfdiff27;
    std::vector<double> _wfdiff28;
    std::vector<double> _wfdiff29;
    std::vector<double> _wfdiff30;
    std::vector<double> _wfdiff31;
    std::vector<double> _th00;
    std::vector<double> _th01;
    std::vector<double> _th02;
    std::vector<double> _th03;
    std::vector<double> _th04;
    std::vector<double> _th05;
    std::vector<double> _th06;
    std::vector<double> _th07;
    std::vector<double> _th08;
    std::vector<double> _th09;
    std::vector<double> _th10;
    std::vector<double> _th11;
    std::vector<double> _th12;
    std::vector<double> _th13;
    std::vector<double> _th14;
    std::vector<double> _th15;
    std::vector<double> _th16;
    std::vector<double> _th17;
    std::vector<double> _th18;
    std::vector<double> _th19;
    std::vector<double> _th20;
    std::vector<double> _th21;
    std::vector<double> _th22;
    std::vector<double> _th23;
    std::vector<double> _th24;
    std::vector<double> _th25;
    std::vector<double> _th26;
    std::vector<double> _th27;
    std::vector<double> _th28;
    std::vector<double> _th29;
    std::vector<double> _th30;
    std::vector<double> _th31;
    std::vector<double> _pmt_wf00;
    std::vector<double> _pmt_wf01;
    std::vector<double> _pmt_wf02;
    std::vector<double> _pmt_wf03;
    std::vector<double> _pmt_wf04;
    std::vector<double> _pmt_wf05;
    std::vector<double> _pmt_wf06;
    std::vector<double> _pmt_wf07;
    std::vector<double> _pmt_wf08;
    std::vector<double> _pmt_wf09;
    std::vector<double> _pmt_wf10;
    std::vector<double> _pmt_wf11;
    std::vector<double> _pmt_wf12;
    std::vector<double> _pmt_wf13;
    std::vector<double> _pmt_wf14;
    std::vector<double> _pmt_wf15;
    std::vector<double> _pmt_wf16;
    std::vector<double> _pmt_wf17;
    std::vector<double> _pmt_wf18;
    std::vector<double> _pmt_wf19;
    std::vector<double> _pmt_wf20;
    std::vector<double> _pmt_wf21;
    std::vector<double> _pmt_wf22;
    std::vector<double> _pmt_wf23;
    std::vector<double> _pmt_wf24;
    std::vector<double> _pmt_wf25;
    std::vector<double> _pmt_wf26;
    std::vector<double> _pmt_wf27;
    std::vector<double> _pmt_wf28;
    std::vector<double> _pmt_wf29;
    std::vector<double> _pmt_wf30;
    std::vector<double> _pmt_wf31;
    std::vector<double> _hit_significance;
    std::vector<double> _hit_time;
    std::vector<double> _all_hit_time;
    std::vector<double> _all_hit_ampl;
    std::vector<int>    _all_hit_pmtn;
    // mc michel info
    int _michel;
    double _michel_E;
    double _michel_t;
    double _michel_x;
    double _michel_y;
    double _michel_z;
    
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
