#ifndef LARLITE_MERGEWF_PADDLES_CXX
#define LARLITE_MERGEWF_PADDLES_CXX

#include "MergeWF_Paddles.h"
#include "DataFormat/opdetwaveform.h"
#include "DataFormat/trigger.h"
#include <cmath>

namespace larlite {

  MergeWF_Paddles::MergeWF_Paddles()
    : _tree(nullptr)
  {
    _name = "MergeWF_Paddles";
    _fout = 0;
  }


  bool MergeWF_Paddles::initialize() {

   fillTree();

    return true;
  }
  
  bool MergeWF_Paddles::analyze(storage_manager* storage) {

    cleanTree();

    // get the trigger info
    auto trig = storage->get_data<trigger>("daq");

    if ( !trig ){
      std::cout << "No trigger found..." << std::endl;
      return false;
    }

    // get the trigger information
    double trig_time = trig->TriggerTime();

    // get the waveform information
    auto ev_opdetwaveform = storage->get_data<event_opdetwaveform>(_producer);

    if ( (!ev_opdetwaveform) or (ev_opdetwaveform->size() == 0) ){
      std::cout << "No waveforms found..." << std::endl;
      return false;
    }

    _event = storage->get_data<event_opdetwaveform>(_producer)->event_id();


    // for each channel, prepare a map that connects
    // 1) pmt -> ev_opdetwaveform index w/ beam-gate
    std::map<int,size_t> beamgateMap;
    // 2) pmt -> ev_opdetwaveform w/ pre-beam-gate (if it exists)
    std::map<int,size_t> prebeamMap;

    // find cosmic-discriminator waveforms
    // in the ~ usec before the beam-gate
    double t_min = 10;

    //std::cout << "fill beamgate map info" << std::endl;
    for (size_t i=0; i < ev_opdetwaveform->size(); i++){

      auto const& wf = ev_opdetwaveform->at(i);

      short pmt = wf.ChannelNumber();

      if (pmt >= 32)
	continue;

      if (wf.size() > 800){
	if (wf.size() > 1501) return true;
	beamgateMap[pmt] = i;
      }

      // get the time relative to the trigger time
      double wf_t = wf.TimeStamp() - trig_time;

      if ( (wf_t < 0) && (wf_t > -2) ){

	prebeamMap[pmt] = i;

	if (wf_t < t_min)
	  t_min = wf_t;
      }// if a couple us before beam-gate

    }// for all wfs    
    //std::cout << "done" << std::endl;
  
    //std::cout << "the minimum time for this event is : " << t_min << std::endl;
    size_t preticks =  (size_t)(-1000*t_min/15.625);
    //std::cout << "pre-ticks for this event: " << preticks << std::endl;
    // if t_min is still 10 -> ignore event (does not have a muon)
    if (t_min > 0)
      return true;

    // keep track of the significance of all hits found along
    // in this event
    _hit_significance.clear();
    _hit_time.clear();

    for (size_t pmt=0; pmt < 32; pmt++){
      /*
      std::cout << "****************************" << std::endl
		<< "    PMT # " << pmt << "      " << std::endl
		<< "****************************" << std::endl;
      */
      // for each PMT create a std::vector<short> the length of the beam-gate + pre-samples
      // due to cosmic discriminator
      //      std::cout << "pad wf" << std::endl;
      std::vector<short> padded_wf(preticks+1501,2048);
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
	//	std::cout << "wf size " << wf.size() << std::endl;
	for (size_t idx=0; idx< 1500; idx++)
	  padded_wf[preticks+idx] = wf[idx];
      }// if we find a beam-gate window
      // if no beam-gate window -> something went wrong
      else
	std::cout << "did not find a beam-gate for this PMT!!! WHAT!!!" << std::endl;
      // add this waveform to the tree
	//      std::cout << "done" << std::endl;
      
      // calculate and subtract baseline for this waveform
      // additionally change scale to PE (for now /20)
      //std::cout << "calculate baseline" << std::endl;
      //std::cout << "padded wf size " << padded_wf.size() << std::endl;
      std::vector<double> baseline;
      //std::cout << "padded wf size " << padded_wf.size() << std::endl;
      std::vector<double> pmt_wf;//(padded_wf.size(),0.);
      //std::cout << "padded wf size " << padded_wf.size() << std::endl;
      pmt_wf.resize(preticks+1501);
      //std::cout << "padded wf size " << padded_wf.size() << std::endl;
      getBaseline(padded_wf,baseline);
      //std::cout << "done1" << std::endl;
      for (size_t idx=0; idx < pmt_wf.size(); idx++){
	pmt_wf[idx] = (double)(padded_wf[idx]) - baseline[idx];
	pmt_wf[idx] /= 20.;
      }
      //std::cout << "done" << std::endl;

      // now find the peaks associated to a muon
      //std::cout << "find muon peaks" << std::endl;
      std::vector<size_t> muonTicks;
      getMuonPeaksTicks(pmt_wf,muonTicks);
      //std::cout << "done" << std::endl;
      // get the times and amplitudes of the muon peaks
      //std::cout << "get peak points" << std::endl;
      std::vector<double> muonPeakT;
      std::vector<double> muonPeakA;
      getPoints(pmt_wf,muonTicks,muonPeakT,muonPeakA);
      //std::cout << "done" << std::endl;
      // get the differential waveform
      //std::cout << "calculate differential vector" << std::endl;
      std::vector<double> wfdiff;
      getDifferential(pmt_wf,wfdiff);
      //std::cout << "done" << std::endl;
      // get the pulse-times from the differential vector
      //std::cout << "get hit ticks" << std::endl;
      std::vector<size_t> hit_ticks;
      std::vector<double> hit_times;
      std::vector<double> hit_PEs;
      findPulseTimes(wfdiff,hit_ticks);
      //std::cout << "done" << std::endl;
      // replace each hit_tick with the local maximum from the wf
      //std::cout << "find local maxima" << std::endl;
      findLocalMaxima(pmt_wf,hit_ticks);
      //std::cout << "done" << std::endl;
      //std::cout << "find points" << std::endl;
      getPoints(pmt_wf,hit_ticks,hit_times,hit_PEs);
      //std::cout << "done" << std::endl;

      for (size_t hh=0; hh < hit_PEs.size(); hh++){
	_all_hit_time.push_back(hit_times[hh]);
	_all_hit_ampl.push_back(hit_PEs[hh]);
	_all_hit_pmtn.push_back(pmt);
      }

      // prepare a vector that carries the expected
      // PE to be seen at each tick based on the muons
      // found
      std::vector<double> expectation(pmt_wf.size(),0);
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
	    /*
	    if (ll > 100){
	      std::cout << "muon time : " << muon_time << std::endl;
	      std::cout << "this tick : " << this_tick << std::endl;
	      std::cout << "this time : " << this_time << std::endl;
	      std::cout << "time : " << this_time-muon_time << std::endl;
	      std::cout << "A = " << 0.3*muonamp << std::endl;
	      std::cout << "ll is " << ll << std::endl;
	    }
	    */
	    expectation[this_tick] += ll;
	    if (expectation[this_tick] > 100)
	      //std::cout << "@ tick " << this_tick << " expectation is " << expectation[this_tick] << std::endl;
	    // if the late-light expectation goes below 1 -> ignore this track from now on
	    if (ll < 1)
	      continue;
	  }// if we've passed this muon along the vector
	}// for all muons
      }// for all time-ticks
      //std::cout << "done" << std::endl;

      //std::cout << "get significance" << std::endl;
      // finally, loop through all hits and measure their significance
      for (size_t idx = 0; idx < hit_times.size(); idx++){
	double hit_t = hit_times[idx];
	double hit_a = hit_PEs[idx];
	size_t tick = hit_t*1000./15.625;
	if (tick > expectation.size())
	  continue;
	double prediction = expectation[tick];
	// if the observed PEs are above the expected calculate the significance
	if (hit_a > prediction){
	  double significance = LateLightProb(prediction, hit_a);
	  // add to the significance tally for this event
	  //std::cout << "new significance: " << significance << std::endl;
	  _hit_significance.push_back(significance);
	  _hit_time.push_back(hit_t);
	}
      }// for all hits
      //std::cout << "done" << std::endl;

      add(pmt,padded_wf,expectation,pmt_wf,wfdiff);      
      //std::cout << "done w/ this pmt" << std::endl << std::endl;

    }// for all PMTs
    
    _tree->Fill();

    //std::cout << "done filling" << std::endl;
      
    return true;
  }

  bool MergeWF_Paddles::finalize() {

    if (_fout){
      if (_tree) _tree->Write();
    }

    return true;
  }


  void MergeWF_Paddles::findLocalMaxima(const std::vector<double>& wf,
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

  
  void MergeWF_Paddles::findPulseTimes(const std::vector<double>& wfdiff,
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
      if ( (t2 > t1) && (t2 > t3) && (t2 > 0.15) ){
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


  void MergeWF_Paddles::getPoints(const std::vector<double>& wf,
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

  
  void MergeWF_Paddles::getPeaks(const std::vector<double>& wf,
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

  
  double MergeWF_Paddles::lateLight(const double& time,
				    const double& A,
				    const double& tau)
  { return A*exp(-time/tau); }

  

  void MergeWF_Paddles::getRMS(const std::vector<double>& wf,
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

  void MergeWF_Paddles::getBaseline(const std::vector<short>& wf,
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


  void MergeWF_Paddles::getDifferential(const std::vector<double>& wf,
					std::vector<double>& wfdiff)
  {

    wfdiff.clear();
    
    for (size_t i=0; i < wf.size()-1; i++)
      wfdiff.push_back( wf[i+1] - wf[i] );

    return;
  }

  
  void MergeWF_Paddles::getMuonPeaksTicks(const std::vector<double>& wf,
					  std::vector<size_t>& muonTicks)
  {

    // dead-time for successive peaks
    size_t deadtime = 20;
    // threshold for new muon (in PE)
    double thresh = 6;

    // ticks for maxima:
    muonTicks.clear();
    
    double currentMaxADC = 0;
    double currentMaxTDC = -(double)(deadtime)-10;

    for (size_t idx=1; idx < wf.size()-1; idx++){
      
      double adc1 = wf[idx-1];
      double adc2 = wf[idx];
      double adc3 = wf[idx+1];
      
      // are we above threshold?
      if (adc2 < thresh)
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

  
  double MergeWF_Paddles::Poisson(const int& lambda, const int& k)
  {

    double fact = 1;
    if (k > 0){
      for (int j=1; j < k+1; j++)
	fact *= j;
    }

    /*
    std::cout << "\t\tlambda^k     = " << pow(lambda,k) << std::endl;
    std::cout << "\t\texp(-lambda) = " << exp(-lambda) << std::endl;
    std::cout << "\t\tfact(k)      = " << fact << std::endl;
    */
    double poiss = pow(lambda,k)*exp(-lambda)/(double)(fact);
    
    //std::cout << "\t\tpoisson " << poiss << std::endl;

    return poiss;
  }

  
  double MergeWF_Paddles::PoissonCDF(const int& lambda, const int& k)
  {

    double prob = 0;
    
    for (int i=0; i < k; i++)
      prob += Poisson(lambda,i);

    //std::cout << "\tpoisson CDF (lambda,i) : " << prob << std::endl;

    return prob;
  }

  double MergeWF_Paddles::InvERF(const double& x)
  {

    double a  = 0.147;
    double pi = 3.14159265359;
    double b  = 2/(a*pi) + log(1-x*x)/2.;
    
    double erfinv = sqrt( sqrt( b*b - log(1-x*x)/a ) - b );

    if (x < 0)
      erfinv *= -1;

    return erfinv;
  }

  double MergeWF_Paddles::LateLightProb(const double& e, const double& o)
  {

    //std::cout << "expectation : " << e << std::endl;
    //std::cout << "observation : " << o << std::endl;

    // round the expectation and observed values to
    // the nearest integer
    int ex = round(e);
    int ob = round(o);
    if (ex < 1)
      ex = 1;

    double prob = PoissonCDF(ex,ob);

    //std::cout << "prob " << prob << std::endl;

    // if the probability is 1 -> return a significance of 20
    if (prob >= 1)
      return 20;

    // otherwise use error function to calculate significance
    double sig = sqrt(2)*InvERF(prob);

    //std::cout << "significance " << sig << std::endl;

    return sig;
  }

  void MergeWF_Paddles::fillTree(){

    if (_tree) delete _tree;
    
    _tree = new TTree("_tree","michel tree");
    _tree->Branch("_event",&_event,"event/I");
    /*
    _tree->Branch("wf00","std::vector<short>",&_wf00);
    _tree->Branch("wf01","std::vector<short>",&_wf01);
    _tree->Branch("wf02","std::vector<short>",&_wf02);
    _tree->Branch("wf03","std::vector<short>",&_wf03);
    _tree->Branch("wf04","std::vector<short>",&_wf04);
    _tree->Branch("wf05","std::vector<short>",&_wf05);
    _tree->Branch("wf06","std::vector<short>",&_wf06);
    _tree->Branch("wf07","std::vector<short>",&_wf07);
    _tree->Branch("wf08","std::vector<short>",&_wf08);
    _tree->Branch("wf09","std::vector<short>",&_wf09);
    _tree->Branch("wf10","std::vector<short>",&_wf10);
    _tree->Branch("wf11","std::vector<short>",&_wf11);
    _tree->Branch("wf12","std::vector<short>",&_wf12);
    _tree->Branch("wf13","std::vector<short>",&_wf13);
    _tree->Branch("wf14","std::vector<short>",&_wf14);
    _tree->Branch("wf15","std::vector<short>",&_wf15);
    _tree->Branch("wf16","std::vector<short>",&_wf16);
    _tree->Branch("wf17","std::vector<short>",&_wf17);
    _tree->Branch("wf18","std::vector<short>",&_wf18);
    _tree->Branch("wf19","std::vector<short>",&_wf19);
    _tree->Branch("wf20","std::vector<short>",&_wf20);
    _tree->Branch("wf21","std::vector<short>",&_wf21);
    _tree->Branch("wf22","std::vector<short>",&_wf22);
    _tree->Branch("wf23","std::vector<short>",&_wf23);
    _tree->Branch("wf24","std::vector<short>",&_wf24);
    _tree->Branch("wf25","std::vector<short>",&_wf25);
    _tree->Branch("wf26","std::vector<short>",&_wf26);
    _tree->Branch("wf27","std::vector<short>",&_wf27);
    _tree->Branch("wf28","std::vector<short>",&_wf28);
    _tree->Branch("wf29","std::vector<short>",&_wf29);
    _tree->Branch("wf30","std::vector<short>",&_wf30);
    _tree->Branch("wf31","std::vector<short>",&_wf31);
    _tree->Branch("wfdiff00","std::vector<double>",&_wfdiff00);
    _tree->Branch("wfdiff01","std::vector<double>",&_wfdiff01);
    _tree->Branch("wfdiff02","std::vector<double>",&_wfdiff02);
    _tree->Branch("wfdiff03","std::vector<double>",&_wfdiff03);
    _tree->Branch("wfdiff04","std::vector<double>",&_wfdiff04);
    _tree->Branch("wfdiff05","std::vector<double>",&_wfdiff05);
    _tree->Branch("wfdiff06","std::vector<double>",&_wfdiff06);
    _tree->Branch("wfdiff07","std::vector<double>",&_wfdiff07);
    _tree->Branch("wfdiff08","std::vector<double>",&_wfdiff08);
    _tree->Branch("wfdiff09","std::vector<double>",&_wfdiff09);
    _tree->Branch("wfdiff10","std::vector<double>",&_wfdiff10);
    _tree->Branch("wfdiff11","std::vector<double>",&_wfdiff11);
    _tree->Branch("wfdiff12","std::vector<double>",&_wfdiff12);
    _tree->Branch("wfdiff13","std::vector<double>",&_wfdiff13);
    _tree->Branch("wfdiff14","std::vector<double>",&_wfdiff14);
    _tree->Branch("wfdiff15","std::vector<double>",&_wfdiff15);
    _tree->Branch("wfdiff16","std::vector<double>",&_wfdiff16);
    _tree->Branch("wfdiff17","std::vector<double>",&_wfdiff17);
    _tree->Branch("wfdiff18","std::vector<double>",&_wfdiff18);
    _tree->Branch("wfdiff19","std::vector<double>",&_wfdiff19);
    _tree->Branch("wfdiff20","std::vector<double>",&_wfdiff20);
    _tree->Branch("wfdiff21","std::vector<double>",&_wfdiff21);
    _tree->Branch("wfdiff22","std::vector<double>",&_wfdiff22);
    _tree->Branch("wfdiff23","std::vector<double>",&_wfdiff23);
    _tree->Branch("wfdiff24","std::vector<double>",&_wfdiff24);
    _tree->Branch("wfdiff25","std::vector<double>",&_wfdiff25);
    _tree->Branch("wfdiff26","std::vector<double>",&_wfdiff26);
    _tree->Branch("wfdiff27","std::vector<double>",&_wfdiff27);
    _tree->Branch("wfdiff28","std::vector<double>",&_wfdiff28);
    _tree->Branch("wfdiff29","std::vector<double>",&_wfdiff29);
    _tree->Branch("wfdiff30","std::vector<double>",&_wfdiff30);
    _tree->Branch("wfdiff31","std::vector<double>",&_wfdiff31);
    _tree->Branch("th00","std::vector<double>",&_th00);
    _tree->Branch("th01","std::vector<double>",&_th01);
    _tree->Branch("th02","std::vector<double>",&_th02);
    _tree->Branch("th03","std::vector<double>",&_th03);
    _tree->Branch("th04","std::vector<double>",&_th04);
    _tree->Branch("th05","std::vector<double>",&_th05);
    _tree->Branch("th06","std::vector<double>",&_th06);
    _tree->Branch("th07","std::vector<double>",&_th07);
    _tree->Branch("th08","std::vector<double>",&_th08);
    _tree->Branch("th09","std::vector<double>",&_th09);
    _tree->Branch("th10","std::vector<double>",&_th10);
    _tree->Branch("th11","std::vector<double>",&_th11);
    _tree->Branch("th12","std::vector<double>",&_th12);
    _tree->Branch("th13","std::vector<double>",&_th13);
    _tree->Branch("th14","std::vector<double>",&_th14);
    _tree->Branch("th15","std::vector<double>",&_th15);
    _tree->Branch("th16","std::vector<double>",&_th16);
    _tree->Branch("th17","std::vector<double>",&_th17);
    _tree->Branch("th18","std::vector<double>",&_th18);
    _tree->Branch("th19","std::vector<double>",&_th19);
    _tree->Branch("th20","std::vector<double>",&_th20);
    _tree->Branch("th21","std::vector<double>",&_th21);
    _tree->Branch("th22","std::vector<double>",&_th22);
    _tree->Branch("th23","std::vector<double>",&_th23);
    _tree->Branch("th24","std::vector<double>",&_th24);
    _tree->Branch("th25","std::vector<double>",&_th25);
    _tree->Branch("th26","std::vector<double>",&_th26);
    _tree->Branch("th27","std::vector<double>",&_th27);
    _tree->Branch("th28","std::vector<double>",&_th28);
    _tree->Branch("th29","std::vector<double>",&_th29);
    _tree->Branch("th30","std::vector<double>",&_th30);
    _tree->Branch("th31","std::vector<double>",&_th31);
    _tree->Branch("pmt_wf00","std::vector<double>",&_pmt_wf00);
    _tree->Branch("pmt_wf01","std::vector<double>",&_pmt_wf01);
    _tree->Branch("pmt_wf02","std::vector<double>",&_pmt_wf02);
    _tree->Branch("pmt_wf03","std::vector<double>",&_pmt_wf03);
    _tree->Branch("pmt_wf04","std::vector<double>",&_pmt_wf04);
    _tree->Branch("pmt_wf05","std::vector<double>",&_pmt_wf05);
    _tree->Branch("pmt_wf06","std::vector<double>",&_pmt_wf06);
    _tree->Branch("pmt_wf07","std::vector<double>",&_pmt_wf07);
    _tree->Branch("pmt_wf08","std::vector<double>",&_pmt_wf08);
    _tree->Branch("pmt_wf09","std::vector<double>",&_pmt_wf09);
    _tree->Branch("pmt_wf10","std::vector<double>",&_pmt_wf10);
    _tree->Branch("pmt_wf11","std::vector<double>",&_pmt_wf11);
    _tree->Branch("pmt_wf12","std::vector<double>",&_pmt_wf12);
    _tree->Branch("pmt_wf13","std::vector<double>",&_pmt_wf13);
    _tree->Branch("pmt_wf14","std::vector<double>",&_pmt_wf14);
    _tree->Branch("pmt_wf15","std::vector<double>",&_pmt_wf15);
    _tree->Branch("pmt_wf16","std::vector<double>",&_pmt_wf16);
    _tree->Branch("pmt_wf17","std::vector<double>",&_pmt_wf17);
    _tree->Branch("pmt_wf18","std::vector<double>",&_pmt_wf18);
    _tree->Branch("pmt_wf19","std::vector<double>",&_pmt_wf19);
    _tree->Branch("pmt_wf20","std::vector<double>",&_pmt_wf20);
    _tree->Branch("pmt_wf21","std::vector<double>",&_pmt_wf21);
    _tree->Branch("pmt_wf22","std::vector<double>",&_pmt_wf22);
    _tree->Branch("pmt_wf23","std::vector<double>",&_pmt_wf23);
    _tree->Branch("pmt_wf24","std::vector<double>",&_pmt_wf24);
    _tree->Branch("pmt_wf25","std::vector<double>",&_pmt_wf25);
    _tree->Branch("pmt_wf26","std::vector<double>",&_pmt_wf26);
    _tree->Branch("pmt_wf27","std::vector<double>",&_pmt_wf27);
    _tree->Branch("pmt_wf28","std::vector<double>",&_pmt_wf28);
    _tree->Branch("pmt_wf29","std::vector<double>",&_pmt_wf29);
    _tree->Branch("pmt_wf30","std::vector<double>",&_pmt_wf30);
    _tree->Branch("pmt_wf31","std::vector<double>",&_pmt_wf31);
    */
    _tree->Branch("hit_significance","std::vector<double>",&_hit_significance);
    _tree->Branch("hit_time","std::vector<double>",&_hit_time);
    _tree->Branch("all_hit_time","std::vector<double>",&_all_hit_time);
    _tree->Branch("all_hit_ampl","std::vector<double>",&_all_hit_ampl);
    _tree->Branch("all_hit_pmtn","std::vector<int>",&_all_hit_pmtn);
    
    return;
  }

  void MergeWF_Paddles::add(const short& pmt,
			    const std::vector<short>& wf,
			    const std::vector<double>& th,
			    const std::vector<double>& pmt_wf,
			    const std::vector<double>& wfdiff)
  {

    if ( (pmt >= 32) or (pmt < 0) )
      return;

    if (pmt ==  0){ _wf00 = wf; _th00 = th; _pmt_wf00 = pmt_wf; _wfdiff00 = wfdiff; }
    if (pmt ==  1){ _wf01 = wf; _th01 = th; _pmt_wf01 = pmt_wf; _wfdiff01 = wfdiff; }
    if (pmt ==  2){ _wf02 = wf; _th02 = th; _pmt_wf02 = pmt_wf; _wfdiff02 = wfdiff; }
    if (pmt ==  3){ _wf03 = wf; _th03 = th; _pmt_wf03 = pmt_wf; _wfdiff03 = wfdiff; }
    if (pmt ==  4){ _wf04 = wf; _th04 = th; _pmt_wf04 = pmt_wf; _wfdiff04 = wfdiff; }
    if (pmt ==  5){ _wf05 = wf; _th05 = th; _pmt_wf05 = pmt_wf; _wfdiff05 = wfdiff; }
    if (pmt ==  6){ _wf06 = wf; _th06 = th; _pmt_wf06 = pmt_wf; _wfdiff06 = wfdiff; }
    if (pmt ==  7){ _wf07 = wf; _th07 = th; _pmt_wf07 = pmt_wf; _wfdiff07 = wfdiff; }
    if (pmt ==  8){ _wf08 = wf; _th08 = th; _pmt_wf08 = pmt_wf; _wfdiff08 = wfdiff; }
    if (pmt ==  9){ _wf09 = wf; _th09 = th; _pmt_wf09 = pmt_wf; _wfdiff09 = wfdiff; }
    if (pmt == 10){ _wf10 = wf; _th10 = th; _pmt_wf10 = pmt_wf; _wfdiff10 = wfdiff; }
    if (pmt == 11){ _wf11 = wf; _th11 = th; _pmt_wf11 = pmt_wf; _wfdiff11 = wfdiff; }
    if (pmt == 12){ _wf12 = wf; _th12 = th; _pmt_wf12 = pmt_wf; _wfdiff12 = wfdiff; }
    if (pmt == 13){ _wf13 = wf; _th13 = th; _pmt_wf13 = pmt_wf; _wfdiff13 = wfdiff; }
    if (pmt == 14){ _wf14 = wf; _th14 = th; _pmt_wf14 = pmt_wf; _wfdiff14 = wfdiff; }
    if (pmt == 15){ _wf15 = wf; _th15 = th; _pmt_wf15 = pmt_wf; _wfdiff15 = wfdiff; }
    if (pmt == 16){ _wf16 = wf; _th16 = th; _pmt_wf16 = pmt_wf; _wfdiff16 = wfdiff; }
    if (pmt == 17){ _wf17 = wf; _th17 = th; _pmt_wf17 = pmt_wf; _wfdiff17 = wfdiff; }
    if (pmt == 18){ _wf18 = wf; _th18 = th; _pmt_wf18 = pmt_wf; _wfdiff18 = wfdiff; }
    if (pmt == 19){ _wf19 = wf; _th19 = th; _pmt_wf19 = pmt_wf; _wfdiff19 = wfdiff; }
    if (pmt == 20){ _wf20 = wf; _th20 = th; _pmt_wf20 = pmt_wf; _wfdiff20 = wfdiff; }
    if (pmt == 21){ _wf21 = wf; _th21 = th; _pmt_wf21 = pmt_wf; _wfdiff21 = wfdiff; }
    if (pmt == 22){ _wf22 = wf; _th22 = th; _pmt_wf22 = pmt_wf; _wfdiff22 = wfdiff; }
    if (pmt == 23){ _wf23 = wf; _th23 = th; _pmt_wf23 = pmt_wf; _wfdiff23 = wfdiff; }
    if (pmt == 24){ _wf24 = wf; _th24 = th; _pmt_wf24 = pmt_wf; _wfdiff24 = wfdiff; }
    if (pmt == 25){ _wf25 = wf; _th25 = th; _pmt_wf25 = pmt_wf; _wfdiff25 = wfdiff; }
    if (pmt == 26){ _wf26 = wf; _th26 = th; _pmt_wf26 = pmt_wf; _wfdiff26 = wfdiff; }
    if (pmt == 27){ _wf27 = wf; _th27 = th; _pmt_wf27 = pmt_wf; _wfdiff27 = wfdiff; }
    if (pmt == 28){ _wf28 = wf; _th28 = th; _pmt_wf28 = pmt_wf; _wfdiff28 = wfdiff; }
    if (pmt == 29){ _wf29 = wf; _th29 = th; _pmt_wf29 = pmt_wf; _wfdiff29 = wfdiff; }
    if (pmt == 30){ _wf30 = wf; _th30 = th; _pmt_wf30 = pmt_wf; _wfdiff30 = wfdiff; }
    if (pmt == 31){ _wf31 = wf; _th31 = th; _pmt_wf31 = pmt_wf; _wfdiff31 = wfdiff; }

    return;
  }

  void MergeWF_Paddles::cleanTree(){

    _wf00.clear();
    _wf01.clear();
    _wf02.clear();
    _wf03.clear();
    _wf04.clear();
    _wf05.clear();
    _wf06.clear();
    _wf07.clear();
    _wf08.clear();
    _wf09.clear();
    _wf10.clear();
    _wf11.clear();
    _wf12.clear();
    _wf13.clear();
    _wf14.clear();
    _wf15.clear();
    _wf16.clear();
    _wf17.clear();
    _wf18.clear();
    _wf19.clear();
    _wf20.clear();
    _wf21.clear();
    _wf22.clear();
    _wf23.clear();
    _wf24.clear();
    _wf25.clear();
    _wf26.clear();
    _wf27.clear();
    _wf28.clear();
    _wf29.clear();
    _wf30.clear();
    _wf31.clear();
    _wfdiff00.clear();
    _wfdiff01.clear();
    _wfdiff02.clear();
    _wfdiff03.clear();
    _wfdiff04.clear();
    _wfdiff05.clear();
    _wfdiff06.clear();
    _wfdiff07.clear();
    _wfdiff08.clear();
    _wfdiff09.clear();
    _wfdiff10.clear();
    _wfdiff11.clear();
    _wfdiff12.clear();
    _wfdiff13.clear();
    _wfdiff14.clear();
    _wfdiff15.clear();
    _wfdiff16.clear();
    _wfdiff17.clear();
    _wfdiff18.clear();
    _wfdiff19.clear();
    _wfdiff20.clear();
    _wfdiff21.clear();
    _wfdiff22.clear();
    _wfdiff23.clear();
    _wfdiff24.clear();
    _wfdiff25.clear();
    _wfdiff26.clear();
    _wfdiff27.clear();
    _wfdiff28.clear();
    _wfdiff29.clear();
    _wfdiff30.clear();
    _wfdiff31.clear();
    _th00.clear();
    _th01.clear();
    _th02.clear();
    _th03.clear();
    _th04.clear();
    _th05.clear();
    _th06.clear();
    _th07.clear();
    _th08.clear();
    _th09.clear();
    _th10.clear();
    _th11.clear();
    _th12.clear();
    _th13.clear();
    _th14.clear();
    _th15.clear();
    _th16.clear();
    _th17.clear();
    _th18.clear();
    _th19.clear();
    _th20.clear();
    _th21.clear();
    _th22.clear();
    _th23.clear();
    _th24.clear();
    _th25.clear();
    _th26.clear();
    _th27.clear();
    _th28.clear();
    _th29.clear();
    _th30.clear();
    _th31.clear();
    _pmt_wf00.clear();
    _pmt_wf01.clear();
    _pmt_wf02.clear();
    _pmt_wf03.clear();
    _pmt_wf04.clear();
    _pmt_wf05.clear();
    _pmt_wf06.clear();
    _pmt_wf07.clear();
    _pmt_wf08.clear();
    _pmt_wf09.clear();
    _pmt_wf10.clear();
    _pmt_wf11.clear();
    _pmt_wf12.clear();
    _pmt_wf13.clear();
    _pmt_wf14.clear();
    _pmt_wf15.clear();
    _pmt_wf16.clear();
    _pmt_wf17.clear();
    _pmt_wf18.clear();
    _pmt_wf19.clear();
    _pmt_wf20.clear();
    _pmt_wf21.clear();
    _pmt_wf22.clear();
    _pmt_wf23.clear();
    _pmt_wf24.clear();
    _pmt_wf25.clear();
    _pmt_wf26.clear();
    _pmt_wf27.clear();
    _pmt_wf28.clear();
    _pmt_wf29.clear();
    _pmt_wf30.clear();
    _pmt_wf31.clear();
    _hit_significance.clear();
    _hit_time.clear();
    _all_hit_time.clear();
    _all_hit_ampl.clear();
    _all_hit_pmtn.clear();

    return;
  }
    
}
#endif
