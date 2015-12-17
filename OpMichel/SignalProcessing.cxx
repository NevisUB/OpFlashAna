#ifndef SIGNALPROCESSING_CXX
#define SIGNALPROCESSING_CXX

#include "SignalProcessing.h"

namespace signalana {

  SignalProcessing::SignalProcessing()
  {
    _hit_PE_differential_thresh = 1;
    _deadtime = 3; // usec
    _muon_PE_thresh = 50;
  }


  void SignalProcessing::getRMS(const std::vector<double>& wf,
			       double& base, double& rms)
  {

    base = 0.;
    rms  = 0.;

    for (auto const& t : wf)
      base += t;
    base /= wf.size();
    
    for (auto const& t : wf)
      rms += (t-base)*(t-base);
    rms = sqrt(rms/(wf.size()-1));

  }

  void SignalProcessing::getBaseline(const std::vector<short>& wf,
				    std::vector<double>& baseline)
  {

    // number of steps to look before-after
    size_t N = 5;

    //std::cout << "wf size : " << wf.size() << std::endl;
    
    baseline.clear();
    baseline = std::vector<double>(wf.size(),0.);

    double currentBase = 0;
    double firstBase   = 0;

    double b,rms;
    
    for (size_t n=N; n < wf.size()-N; n++){

      std::vector<double> localWF(wf.begin()+n-N,wf.begin()+n+N);
      getRMS(localWF,b,rms);

      if (rms < 2.){
	baseline[n] = b;
	currentBase = b;
	if (firstBase == 0)
	  firstBase = b;
      }
      else
	baseline[n] = currentBase;
    }// for all ticks

    // pad the last few entries that may not have a baseline value
    for (size_t n=0; n < N; n++)
      baseline[wf.size()-N+n] = currentBase;
    // pad the first few entries...
    size_t t=0;
    while (baseline[t] == 0){
      baseline[t] = firstBase;
      t += 1;
    }

    return;
  }


  void SignalProcessing::getDifferential(const std::vector<double>& wf,
					std::vector<double>& wfdiff)
  {

    wfdiff.clear();
    
    for (size_t i=0; i < wf.size()-1; i++)
      wfdiff.push_back( wf[i+1] - wf[i] );

    return;
  }



  void SignalProcessing::findPulseTimes(const std::vector<double>& wfdiff,
				       std::vector<size_t>& pulseTicks)
  {

    pulseTicks.clear();
    
    if (wfdiff.size() == 0){
      std::cout << "differential waveform has size 0!" << std::endl;
      return;
    }
    
    for (size_t i=1; i < wfdiff.size()-1; i++){

      double t1 = wfdiff[i-1];
      double t2 = wfdiff[i];
      double t3 = wfdiff[i+1];

      // first find a maximum above a certain threshold
      if ( (t2 > t1) && (t2 > t3) && (t2 > _hit_PE_differential_thresh) ){
	// now open a window for several ticks
	// to find the zero-crossing
	//std::cout << "found a maximum @ time " << i*15.625/1000. << std::endl;
	for (size_t j=0; j < 20; j++){
	  // don't go out of bounds
	  if (i+j+1 >= wfdiff.size()){
	    //std::cout << "reached end of wf " << std::endl;
	    break;
	  }
	  auto const& z1 = wfdiff[i+j];
	  auto const& z2 = wfdiff[i+j+1];
	  // if one above and the other below baseline
	  if ( (z1 > 0) && (z2 <= 0) ){
	    pulseTicks.push_back(i+j);
	    //std::cout << "saving a maximum @ time " << i*15.625/1000. << std::endl;
	    break;
	  }
	}
	// if this far -> no crossing found
	//std::cout << "no zero-crossing was found @ max " << i*15.625/1000. << std::endl;
      }
    }
    return;
  }


  void SignalProcessing::findLocalMaxima(const std::vector<double>& wf,
					std::vector<size_t>& hitTicks)
  {

    for (size_t i=0; i < hitTicks.size(); i++){

      size_t hitT = hitTicks[i];
      size_t minTick = 0;
      if (int(hitT) - 3 > 0) minTick = hitT-3;
      size_t maxTick = hitT+4;
      if (maxTick >= wf.size()) maxTick = wf.size()-1;
      size_t bestTick = hitT;
      double bestADC  = wf[hitT];
      for (size_t j = minTick; j < maxTick; j++){
	if (wf[j] > bestADC){
	  bestTick = j;
	  bestADC  = wf[j];
	}
      }
      // update what the best-tick is
      hitTicks[i] = bestTick;
    }// for all hit ticks found
      
    return;
  }

  


  void SignalProcessing::getPoints(const std::vector<double>& wf,
				  const std::vector<size_t>& pts,
				  std::vector<double>& times,
				  std::vector<double>& vals)
  {

    times.clear();
    vals.clear();
    
    for (auto const& pt : pts){
      if (pt >= wf.size())
	std::cout << "ERROR -> pt out out bounds!" << std::endl;
      times.push_back(pt*15.6/1000.);
      vals.push_back(wf[pt]);
    }
    return;
  }

  
  void SignalProcessing::getPeaks(const std::vector<double>& wf,
				 const std::vector<size_t>& ticks,
				 std::vector<double>& times,
				 std::vector<double>& PEs)
  {
    
    times.clear();
    PEs.clear();

    for (size_t i=0; i < ticks.size(); i++){

      // scan neighboring points for the true maximum
      size_t maxTDC = ticks[i];
      size_t TDC = maxTDC;
      double maxADC = wf[maxTDC];

      for (size_t j=-2; j<3; j++){
	// make sure we don't go out of bounds
	if ( (j >= TDC) or ((TDC+j) > wf.size()) )
	  continue;
	if ( wf[TDC+j] > maxADC){
	  maxADC = wf[TDC+j];
	  maxTDC = TDC+j;
	}
      }
      times.push_back((maxTDC)*15.625/1000.);
      PEs.push_back(maxADC);
    }// for all peaks

    return;
  }

  void SignalProcessing::getMuonPeaksTicks(const std::vector<double>& wf,
					  std::vector<size_t>& muonTicks)
  {

    // dead-time for successive peaks (in ticks)
    size_t deadtime = (int)_deadTime*1000./15.625;

    // ticks for maxima:
    muonTicks.clear();
    
    double currentMaxADC = 0;
    double currentMaxTDC = -(double)(deadtime)-10;

    for (size_t idx=1; idx < wf.size()-1; idx++){
      
      double adc1 = wf[idx-1];
      double adc2 = wf[idx];
      double adc3 = wf[idx+1];
      
      // are we above threshold?
      if (adc2 < _muon_PE_thresh)
	continue;

      // are we at a maximum?
      if ( (adc2 > adc1) && (adc2 <= adc3) ){

	// if we have passed the dead-time
	// create a new maximum!
	if ( (idx-currentMaxTDC) > deadtime ){
	  currentMaxADC = adc2;
	  currentMaxTDC = idx;
	  muonTicks.push_back(currentMaxTDC);
	}

	// otherwise, if we have not passed the dead-time
	else{
	  // is it a new maximum?
	  // update the previous muon peak
	  if (adc2 > currentMaxADC){
	    currentMaxADC = adc2;
	    currentMaxTDC = idx;
	    muonTicks.at(muonTicks.size()-1) = currentMaxTDC;
	  }// if it is a new maximum
	}// if we have not passed the dead-time

      }// if we are at a maximum

    }// for all ticks in wf

    return;
  }


  void SignalProcessing::applyMuonDeadTime(std::vector<double>& peakTimes,
					    std::vector<double>& peakAmps)
  {

    std::vector<double> peakTimesAfter;
    std::vector<double> peakAmpsAfter;
    
    // unless we find a larger pulse -> apply a dead-time
    double maxTime = 0;
    double maxAmp  = 0;
    for (size_t i=0; i < peakTimes.size(); i++){

      if ( ( ((peakTimes[i]-maxTime)*15.625/1000.) > _deadTime) or
	   ( peakAmps[i] > maxAmp) ){
	maxAmp = peakAmps[i];
	maxTime = peakTimes[i];
	peakTimesAfter.push_back(maxAmp);
	peakAmpsAfter.push_back(maxTime);
      }// if it passes the dead-time cut
    }// for all peaks found
    
    peakTimes = peakTimesAfter;
    peakAmps  = peakAmpsAfter;

    return;
  }


}

#endif
