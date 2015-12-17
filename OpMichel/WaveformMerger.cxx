#ifndef LARLITE_WAVEFORMMERGER_CXX
#define LARLITE_WAVEFORMMERGER_CXX

#include "WaveformMerger.h"
#include "DataFormat/opdetwaveform.h"
#include "DataFormat/trigger.h"
#include <map>

namespace larlite {

  WaveformMerger::WaveformMerger()
  {
    _name = "WaveformMerger";
    _fout = 0;
    _sampling = 0.015625;
    _t_min = 0;
    _t_max = 23.4;
    _baseline = 0;
    _verbose = false;
  }

  bool WaveformMerger::initialize() {

    return true;
  }
  
  bool WaveformMerger::analyze(storage_manager* storage) {


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

    auto event  = ev_opdetwaveform->event_id();
    auto subrun = ev_opdetwaveform->subrun();
    auto run    = ev_opdetwaveform->run();

    if (_verbose)
      std::cout << std::endl << std::endl << "[run,subrun,event] -> [ " << run << ", " << subrun << ", " << event << "]" << std::endl;

    // create a multi-map object to store
    // pmt -> [ idx of wf for that pmt ]
    std::multimap<int,size_t> wfMap;

    // keep track of the smallest and largest times in the event
    // for th waveforms we intend to keep
    // these will be used to determine the size of the merged wf
    double wf_t_min = 4600.; // us
    double wf_t_max =    0.; // us

    // loop through all waveforms and select
    // those in the time-range specified by the user
    for (size_t i=0; i < ev_opdetwaveform->size(); i++){
    
      auto const& wf = ev_opdetwaveform->at(i);

      short pmt = wf.ChannelNumber();

      if (pmt >= 32)
	continue;

      // filter by time and add to multi-map

      // get the wf's start time relative to the trigger time
      double wf_start_t = wf.TimeStamp() - trig_time;
      // get the wf's end time relative to the trigger time
      double wf_end_t   = wf_start_t + wf.size()*_sampling;

      if ( ( (wf_end_t  > _t_min) and (wf_end_t < _t_max) ) or 
	   ( (wf_start_t < _t_max) and (wf_start_t > _t_min) ) ){

      if (_verbose)
	std::cout << "\twf for PMT " << pmt << "\tsize : " << wf.size()
		  << "\tstart : " << wf_start_t << "\tend : " << wf_end_t << std::endl;
	
	wfMap.emplace(pmt,i);

	if (wf_end_t   > wf_t_max) { wf_t_max = wf_end_t; }
	if (wf_start_t < wf_t_min) { wf_t_min = wf_start_t; }

      }// if wf in the correct time-range

    }// for all wfs

    // create output ef_opdetwaveform
    auto out_opdetwaveform = storage->get_data<event_opdetwaveform>("mergedwf");

    storage->set_id(run,subrun,event);

    // have we found any wfs to save? if no, return
    if (wfMap.size() == 0){
      if (_verbose)
	std::cout << "NO WF FOR THIS EVENT -> SKIP " << std::endl;
      return true;
    }

    // now that we know the time-interval needed
    // create a wf of that size
    size_t ticks = (size_t)((wf_t_max-wf_t_min)/_sampling);

    if (_verbose){
      std::cout << "time-range : [" << wf_t_min << ", " << wf_t_max << "]" << std::endl;
      std::cout << "number of ticks : " << ticks << std::endl;
    }

    // loop through all PMTs
    for (int pmt=0; pmt < 32; pmt++){

      larlite::opdetwaveform mergedWf;
      mergedWf.resize(ticks,_baseline);
      if (_verbose)
	std::cout << "merged wf has " << mergedWf.size() << " ticks" << std::endl;

      // run through multimap entries for this PMT and add their ticks
      // at the appropriate times
      
      // find this pmt in the multimap
      if (wfMap.find(pmt) == wfMap.end())
	continue;

      if (_verbose)
	std::cout << " found wfs for pmt " << pmt << std::endl;
      
      std::pair <std::multimap<int,size_t>::iterator, std::multimap<int,size_t>::iterator> ret;
      ret = wfMap.equal_range(pmt);

      for (std::multimap<int,size_t>::iterator it = ret.first; it != ret.second; it++){
	
	// get the wf
	auto const& wf = ev_opdetwaveform->at(it->second);

	// get start tick for this wf
	size_t start_tick = (wf.TimeStamp() - trig_time - wf_t_min) / _sampling;

	if (_verbose){
	  std::cout << "wf start time : " << wf.TimeStamp() - trig_time << std::endl;
	  std::cout << "wf start tick : " << start_tick << "\t end tick : " << start_tick + wf.size() << std::endl;
	}

	size_t t_max = wf.size();
	if ( (wf.size()+start_tick) >= mergedWf.size() )
	  t_max = mergedWf.size() - start_tick - 1;

	for (size_t tick = 0; tick < t_max; tick++)
	  mergedWf[tick+start_tick] = wf[tick];
	
	mergedWf.SetTimeStamp(wf_t_min);
	mergedWf.SetChannelNumber(pmt);

	out_opdetwaveform->push_back(mergedWf);
	
      } // for all found wfs for this pmt
      

    }// for all pmts
	
    return true;
  }

  bool WaveformMerger::finalize() {

    return true;
  }

}
#endif
