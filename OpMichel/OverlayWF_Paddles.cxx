#ifndef LARLITE_OVERLAYWF_PADDLES_CXX
#define LARLITE_OVERLAYWF_PADDLES_CXX

#include "OverlayWF_Paddles.h"
#include "DataFormat/opdetwaveform.h"
#include "DataFormat/trigger.h"
#include "DataFormat/mcshower.h"
#include <cmath>

namespace larlite {

  OverlayWF_Paddles::OverlayWF_Paddles()
    : _tree(nullptr)
  {
    _name = "OverlayWF_Paddles";
    _PMTproducer  = "saturation";
    _TRIGproducer = "daq";
    _useMC = false;
    _fout = 0;
    _muon_PE_thresh = 1000;
    _deadTime = 3;
  }


  bool OverlayWF_Paddles::initialize() {

   fillTree();

    return true;
  }
  
  bool OverlayWF_Paddles::analyze(storage_manager* storage) {

    cleanTree();


    // get the MC info, if requested
    if (_useMC){

      auto ev_mcsh = storage->get_data<event_mcshower>("mcreco");
      
      if ( (!ev_mcsh) or (ev_mcsh->size() == 0) ){
	std::cout << "No mcshower found..." << std::endl;
	return false;
      }
      else{
	
	for (size_t n=0; n < ev_mcsh->size(); n++){
	  
	  auto const& shr = ev_mcsh->at(n);
	  
	  if (shr.Process() == "muMinusCaptureAtRest"){
	    //std::cout << "found a michel" << std::endl;
	    _michel = 1;
	    _michel_E = shr.DetProfile().E();
	    _michel_t = shr.DetProfile().T();
	    _michel_x = shr.DetProfile().X();
	    _michel_y = shr.DetProfile().Y();
	    _michel_z = shr.DetProfile().Z();
	    
	  }// if a michel
	}// for all mcshowers
      }// if we have mcshowers
    }// if use MC

    if ( (_useMC) && (_michel == 0) )
      return true;

    // get the trigger info
    auto trig = storage->get_data<trigger>(_TRIGproducer);

    if ( !trig ){
      std::cout << "No trigger found..." << std::endl;
      return false;
    }

    // get the trigger information
    double trig_time = trig->TriggerTime();

    // get the waveform information
    auto ev_opdetwaveform = storage->get_data<event_opdetwaveform>(_PMTproducer);

    if ( (!ev_opdetwaveform) or (ev_opdetwaveform->size() == 0) ){
      std::cout << "No waveforms found..." << std::endl;
      return false;
    }

    _event = storage->get_data<event_opdetwaveform>(_PMTproducer)->event_id();


    // for each channel, prepare a map that connects
    // 1) pmt -> ev_opdetwaveform index w/ beam-gate
    std::map<int,size_t> beamgateMap;
    // 2) pmt -> ev_opdetwaveform w/ pre-beam-gate (if it exists)
    std::map<int,size_t> prebeamMap;

    // find cosmic-discriminator waveforms
    // in the ~ usec before the beam-gate
    double t_min = 10;

    // keep track of the beam-gate size
    size_t beamGateTicks = 0;

    //std::cout << "fill beamgate map info" << std::endl;
    for (size_t i=0; i < ev_opdetwaveform->size(); i++){

      auto const& wf = ev_opdetwaveform->at(i);

      short pmt = wf.ChannelNumber();

      if (pmt >= 32)
	continue;

      if (wf.size() > 800){
	if (wf.size() > 1511){
	  std::cout << "Found beam gate that is too large!" << std::endl;
	  return true;
	}
	if (wf.size() > beamGateTicks) { beamGateTicks = wf.size(); }
	beamgateMap[pmt] = i;
      }

      // get the time relative to the trigger time
      double wf_t = wf.TimeStamp() - trig_time;

      if ( (wf_t < 0) && (wf_t > -2) ){

	//std::cout << "Found pre-beam pulse!" << std::endl;
	prebeamMap[pmt] = i;

	if (wf_t < t_min)
	  t_min = wf_t;
      }// if a couple us before beam-gate

    }// for all wfs    
    //std::cout << "done" << std::endl;
  
    //std::cout << "the minimum time for this event is : " << t_min << std::endl;
    size_t preticks = 0;
    if (t_min < 0)
      preticks =  (size_t)(-1000*t_min/15.625);
    //std::cout << "pre-ticks for this event: " << preticks << std::endl;
    // if t_min is still 10 -> ignore event (does not have a muon)
    if ( (t_min > 0) and (_useMC == false) )
      return true;

    // keep track of the significance of all hits found along
    // in this event
    _hit_significance.clear();
    _hit_time.clear();
    
    // prepare vector to hold overlay of all WFs
    std::vector<double> overlay_wf(preticks+beamGateTicks,0);

    // do baseline subtraction, then add wf of each PMT to the overlay
    for (size_t pmt=0; pmt < 32; pmt++){

      std::vector<short> padded_wf(preticks+beamGateTicks,2048);
      // get the pre-beam-gate, if it exists
      if ( prebeamMap.find(pmt) != prebeamMap.end() ){
	//	std::cout << "found pre-beam stuff" << std::endl;
	auto const& wf = ev_opdetwaveform->at(prebeamMap[pmt]);
	int wf_tick = (int)(-1000*(wf.TimeStamp() - trig_time)/15.625);
	for (size_t idx=0; idx< wf.size(); idx++){
	  size_t this_tick = preticks-wf_tick+idx;
	  padded_wf[this_tick] = wf[idx];
	}// for all ticks
      }// if, for this pmt, we found a pre-beamgate window

      // get the in-beamgate window
      if ( beamgateMap.find(pmt) != beamgateMap.end() ){
	auto const& wf = ev_opdetwaveform->at(beamgateMap[pmt]);
	for (size_t idx=0; idx< 1500; idx++)
	  padded_wf[preticks+idx] = wf[idx];
      }// if we find a beam-gate window
      // if no beam-gate window -> something went wrong
      else
	std::cout << "did not find a beam-gate for this PMT!!! WHAT!!!" << std::endl;
      // add this waveform to the tree
      // calculate and subtract baseline for this waveform
      // additionally change scale to PE (for now /20)
      std::vector<double> baseline;
      std::vector<double> pmt_wf;//(padded_wf.size(),0.);
      pmt_wf.resize(preticks+beamGateTicks);
      getBaseline(padded_wf,baseline);
      for (size_t idx=0; idx < (preticks+beamGateTicks); idx++){
	double pe = (double)(padded_wf[idx]) - baseline[idx];
	overlay_wf[idx] += pe;
      }
    }// for all PMTs

    // now calculate the hits on the overlay wf

    // now find the peaks associated to a muon
    std::vector<size_t> muonTicks;
    getMuonPeaksTicks(overlay_wf,muonTicks);
    // get the times and amplitudes of the muon peaks
    std::vector<double> muonPeakT;
    std::vector<double> muonPeakA;
    getPoints(overlay_wf,muonTicks,muonPeakT,muonPeakA);
    // implement a dead-time for muon pulses
    //applyMuonDeadTime(muonPeakT,muonPeakA);
    // get the differential waveform
    std::vector<double> wfdiff;
    getDifferential(overlay_wf,wfdiff);
    // get the pulse-times from the differential vector
    std::vector<size_t> hit_ticks;
    std::vector<double> hit_times;
    std::vector<double> hit_PEs;
    findPulseTimes(wfdiff,hit_ticks);
    // replace each hit_tick with the local maximum from the wf
    findLocalMaxima(overlay_wf,hit_ticks);
    getPoints(overlay_wf,hit_ticks,hit_times,hit_PEs);

    for (size_t hh=0; hh < hit_PEs.size(); hh++){
      _all_hit_time.push_back(hit_times[hh]);
      _all_hit_ampl.push_back(hit_PEs[hh]);
    }
    
    // prepare a vector that carries the expected
    // PE to be seen at each tick based on the muons
    // found
    std::vector<double> expectation(overlay_wf.size(),30);
    //std::cout << "find expectation" << std::endl;
    // loop through all ticks
    for (size_t this_tick = 0; this_tick < expectation.size(); this_tick++){
      // for all muons
      for (size_t m=0; m < muonPeakT.size(); m++){
	double muonamp  = muonPeakA[m];
	double muon_time = muonPeakT[m] - 0.05; // account for shaping time
	// are we passed this time in the vector?
	// if not -> this exponential does not contribute
	size_t muon_tick = (size_t)(muon_time*1000./15.625);
	double this_time = this_tick*15.625/1000.;
	if ( (this_time > muon_time) ){
	  // current time based on tick
	  double ll = lateLight(this_time-muon_time,0.3*muonamp,1.5);
	  expectation[this_tick] += ll;
	}// if we've passed this muon along the vector
      }// for all muons
    }// for all time-ticks

    // finally, loop through all hits and measure their significance
    for (size_t idx = 0; idx < hit_times.size(); idx++){
      double hit_t = hit_times[idx];
      double hit_a = hit_PEs[idx];
      // if PE is larger than what counts as a muon -> don't calculate significance for this)
      if (hit_a > 3000)
	continue;
      size_t tick = hit_t*1000./15.625;
      if (tick > expectation.size())
	continue;
      double prediction = expectation[tick];
      // if the observed PEs are above the expected calculate the significance
      if (hit_a > prediction){
	double significance = LateLightProb(prediction, hit_a);
	// add to the significance tally for this event
	_hit_significance.push_back(significance);
	_hit_time.push_back(hit_t);
      }
    }// for all hits

    //std::cout << "overlay size : " << overlay_wf.size() << std::endl;
    
    _merged_wf = overlay_wf;
    _theory_wf = expectation;

    _tree->Fill();

    //std::cout << "done filling" << std::endl;
      
    return true;
  }

  bool OverlayWF_Paddles::finalize() {

    if (_fout){
      if (_tree) _tree->Write();
    }

    return true;
  }


  void OverlayWF_Paddles::findLocalMaxima(const std::vector<double>& wf,
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

  
  void OverlayWF_Paddles::findPulseTimes(const std::vector<double>& wfdiff,
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
      if ( (t2 > t1) && (t2 > t3) && (t2 > 20.15) ){
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


  void OverlayWF_Paddles::getPoints(const std::vector<double>& wf,
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

  
  void OverlayWF_Paddles::getPeaks(const std::vector<double>& wf,
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

  
  double OverlayWF_Paddles::lateLight(const double& time,
				    const double& A,
				    const double& tau)
  { return A*exp(-time/tau); }

  

  void OverlayWF_Paddles::getRMS(const std::vector<double>& wf,
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

  void OverlayWF_Paddles::getBaseline(const std::vector<short>& wf,
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


  void OverlayWF_Paddles::getDifferential(const std::vector<double>& wf,
					std::vector<double>& wfdiff)
  {

    wfdiff.clear();
    
    for (size_t i=0; i < wf.size()-1; i++)
      wfdiff.push_back( wf[i+1] - wf[i] );

    return;
  }

  
  void OverlayWF_Paddles::getMuonPeaksTicks(const std::vector<double>& wf,
					  std::vector<size_t>& muonTicks)
  {

    // dead-time for successive peaks
    size_t deadtime = 200;

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


  void OverlayWF_Paddles::applyMuonDeadTime(std::vector<double>& peakTimes,
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

  
  double OverlayWF_Paddles::Poisson(const int& lambda, const int& k)
  {

    double fact = 1;
    if (k > 0){
      for (int j=1; j < k+1; j++)
	fact *= j;
    }

    double poiss = pow(lambda,k)*exp(-lambda)/(double)(fact);

    return poiss;
  }

  double OverlayWF_Paddles::PoissonCDF(const int& lambda, const int& k)
  {

    double prob = 0;
    
    for (int i=0; i < k; i++)
      prob += Poisson(lambda,i);

    return prob;
  }

  double OverlayWF_Paddles::InvERF(const double& x)
  {

    double a  = 0.147;
    double pi = 3.14159265359;
    double b  = 2/(a*pi) + log(1-x*x)/2.;
    
    double erfinv = sqrt( sqrt( b*b - log(1-x*x)/a ) - b );

    if (x < 0)
      erfinv *= -1;

    return erfinv;
  }

  double OverlayWF_Paddles::LateLightProb(const double& e, const double& o)
  {

    // round the expectation and observed values to
    // the nearest integer
    int ex = round(e);
    int ob = round(o);
    if (ex < 1)
      ex = 1;

    //std::cout << "expecation : " << e << "\t observed : " << o << std::endl;;

    double prob = 0.;

    // if the expectation is very large -> approximate the poisson
    // with a gaussian of mean and sigma = lambda
    if (ex > 50){
      if (ob < ex) { std::cout << "ERROR : we should not call when observed is < expectation" << std::endl; }
      double xprime = (o-e)/sqrt(2*e);
      prob = erf(xprime);
    }
    else
      prob = PoissonCDF(ex,ob);
    
    //std::cout << "\tCDF prob : " << prob << std::endl;
    
    // if the probability is 1 -> return a significance of 20
    if (prob >= 1)
      return 20;
    
    // otherwise use error function to calculate significance
    double sig = sqrt(2)*InvERF(prob);
    
    //std::cout << "\t significance : " << sig << std::endl;

    return sig;
  }

  void OverlayWF_Paddles::fillTree(){

    if (_tree) delete _tree;
    
    _tree = new TTree("_tree","michel tree");
    _tree->Branch("_event",&_event,"event/I");

    _tree->Branch("hit_significance","std::vector<double>",&_hit_significance);
    _tree->Branch("hit_time","std::vector<double>",&_hit_time);
    _tree->Branch("merged_wf","std::vector<double>",&_merged_wf);
    _tree->Branch("theory_wf","std::vector<double>",&_theory_wf);
    _tree->Branch("all_hit_time","std::vector<double>",&_all_hit_time);
    _tree->Branch("all_hit_ampl","std::vector<double>",&_all_hit_ampl);
    _tree->Branch("all_hit_pmtn","std::vector<int>",&_all_hit_pmtn);
    _tree->Branch("_michel",&_michel,"michel/I");
    _tree->Branch("_michel_E",&_michel_E,"michel_E/D");
    _tree->Branch("_michel_t",&_michel_t,"michel_t/D");
    _tree->Branch("_michel_x",&_michel_x,"michel_x/D");
    _tree->Branch("_michel_y",&_michel_y,"michel_y/D");
    _tree->Branch("_michel_z",&_michel_z,"michel_z/D");
    
    return;
  }

  void OverlayWF_Paddles::cleanTree(){

    _hit_significance.clear();
    _hit_time.clear();
    _merged_wf.clear();
    _theory_wf.clear();
    _all_hit_time.clear();
    _all_hit_ampl.clear();
    _all_hit_pmtn.clear();
    _michel = 0;
    _michel_E = -1;
    _michel_t = -1;
    _michel_x = -100;
    _michel_y = -200;
    _michel_z = -100;

    return;
  }
    
}
#endif
