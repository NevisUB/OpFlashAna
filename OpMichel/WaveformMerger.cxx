#ifndef LARLITE_WAVEFORMMERGER_CXX
#define LARLITE_WAVEFORMMERGER_CXX

#include "WaveformMerger.h"
#include "DataFormat/opdetwaveform.h"
#include "DataFormat/trigger.h"
#include <map>

namespace larlite {

  WaveformMerger::WaveformMerger()
    :
    _tree(nullptr)
  {
    _name = "WaveformMerger";
    _fout = 0;
    _sampling = 0.015625;
    _t_min = 0;
    _t_max = 23.4;
    _baseline = 0;
    _write = false;
    _verbose = false;
  }

  bool WaveformMerger::initialize() {

    _wfs.resize(32);

    _tree = new TTree("_tree","merged waveform tree");
    _tree->Branch("_event",&_event,"event/I");
    _tree->Branch("_run",&_run,"run/I");
    _tree->Branch("_subrun",&_subrun,"subrun/I");
    _tree->Branch("_timestamp",&_timestamp,"timestamp/D");
    _tree->Branch("wf00","std::vector<short>",&_wfs[0]);
    _tree->Branch("wf01","std::vector<short>",&_wfs[1]);
    _tree->Branch("wf02","std::vector<short>",&_wfs[2]);
    _tree->Branch("wf03","std::vector<short>",&_wfs[3]);
    _tree->Branch("wf04","std::vector<short>",&_wfs[4]);
    _tree->Branch("wf05","std::vector<short>",&_wfs[5]);
    _tree->Branch("wf06","std::vector<short>",&_wfs[6]);
    _tree->Branch("wf07","std::vector<short>",&_wfs[7]);
    _tree->Branch("wf08","std::vector<short>",&_wfs[8]);
    _tree->Branch("wf09","std::vector<short>",&_wfs[9]);
    _tree->Branch("wf10","std::vector<short>",&_wfs[10]);
    _tree->Branch("wf11","std::vector<short>",&_wfs[11]);
    _tree->Branch("wf12","std::vector<short>",&_wfs[12]);
    _tree->Branch("wf13","std::vector<short>",&_wfs[13]);
    _tree->Branch("wf14","std::vector<short>",&_wfs[14]);
    _tree->Branch("wf15","std::vector<short>",&_wfs[15]);
    _tree->Branch("wf16","std::vector<short>",&_wfs[16]);
    _tree->Branch("wf17","std::vector<short>",&_wfs[17]);
    _tree->Branch("wf18","std::vector<short>",&_wfs[18]);
    _tree->Branch("wf19","std::vector<short>",&_wfs[19]);
    _tree->Branch("wf20","std::vector<short>",&_wfs[20]);
    _tree->Branch("wf21","std::vector<short>",&_wfs[21]);
    _tree->Branch("wf22","std::vector<short>",&_wfs[22]);
    _tree->Branch("wf23","std::vector<short>",&_wfs[23]);
    _tree->Branch("wf24","std::vector<short>",&_wfs[24]);
    _tree->Branch("wf25","std::vector<short>",&_wfs[25]);
    _tree->Branch("wf26","std::vector<short>",&_wfs[26]);
    _tree->Branch("wf27","std::vector<short>",&_wfs[27]);
    _tree->Branch("wf28","std::vector<short>",&_wfs[28]);
    _tree->Branch("wf29","std::vector<short>",&_wfs[29]);
    _tree->Branch("wf30","std::vector<short>",&_wfs[30]);
    _tree->Branch("wf31","std::vector<short>",&_wfs[31]);

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

    _event  = ev_opdetwaveform->event_id();
    _subrun = ev_opdetwaveform->subrun();
    _run    = ev_opdetwaveform->run();

    if (_verbose)
      std::cout << std::endl << std::endl << "[run,subrun,event] -> [ " << _run << ", " << _subrun << ", " << _event << "]" << std::endl;

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

    storage->set_id(_run,_subrun,_event);

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
	  std::cout << "wf start time : " << wf.TimeStamp() - trig_time
		    << "\tstart tick : " << start_tick << "\t end tick : " << start_tick + wf.size() << std::endl;
	}

	size_t t_max = wf.size();
	if ( (wf.size()+start_tick) >= mergedWf.size() )
	  t_max = mergedWf.size() - start_tick - 1;
	
	for (size_t tick = 0; tick < t_max; tick++)
	  mergedWf[tick+start_tick] = wf[tick];

      }// for all wfs for a given channel
	
      _wfs[pmt] = (std::vector<short>)mergedWf;
      
      mergedWf.SetTimeStamp(wf_t_min);
      mergedWf.SetChannelNumber(pmt);
      
      out_opdetwaveform->push_back(mergedWf);
      
    }// for all pmts

    _timestamp = wf_t_min;

    // write tree
    _tree->Fill();
	
    return true;
  }

  bool WaveformMerger::finalize() {

    if (_fout && _write){
      _fout->cd();
      if (_tree) { _tree->Write(); }
    }

    return true;
  }

}
#endif
