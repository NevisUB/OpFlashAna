/**
 * \file SignalProcessing.h
 *
 * \ingroup OpMichel
 * 
 * \brief Class def header for a class SignalProcessing
 *
 * @author davidc1
 */

/** \addtogroup OpMichel

    @{*/
#ifndef SIGNALPROCESSING_H
#define SIGNALPROCESSING_H

#include <iostream>
#include <vector>
#include <math.h>

namespace signalana {

  /**
     \class SignalProcessing
     User defined class SignalProcessing ... these comments are used to generate
     doxygen documentation!
  */
  class SignalProcessing{
    
  public:
    
    /// Default constructor
    SignalProcessing();
    
    /// Default destructor
    ~SignalProcessing(){}


    // **** SETTER FUNCTIONS ***

    // PE threshold to claim a hit
    void setHitPEDifferentialThresh(double pe) { _hit_PE_differential_thresh = pe; }
    // set dead-time for muons
    void setDeadTime(double t) { _deadTime = t; }
    // set muon PE threshold
    void setMuonPEThresh(double pe) { _muon_PE_thresh = pe; }

    // **** ALGORITHMS ****

    // function to get rms and baseline
    void getRMS(const std::vector<double>& wf,
		double& base, double& rms);


    // function to get the rolling baseline from a vector
    void getBaseline(const std::vector<short>& wf,
		     std::vector<double>& baseline);


    // get the differential waveform
    void getDifferential(const std::vector<double>& wf,
			 std::vector<double>& wfdiff);


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
    

    // function to find all muon peaks in a waveform
    void getMuonPeaksTicks(const std::vector<double>& wf,
			   std::vector<size_t>& muonTicks);


    // function that applies a dead-time for the muon peaks
    void applyMuonDeadTime(std::vector<double>& peakTimes,
			   std::vector<double>& peakAmps);
    // dead-time for muons
    double _deadTime;

  private:

    double _hit_PE_differential_thresh;
    double _deadtime;
    double _muon_PE_thresh;

  };

}

#endif
/** @} */ // end of doxygen group 

