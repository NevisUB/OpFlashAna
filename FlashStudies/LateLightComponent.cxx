#ifndef LARLITE_LATELIGHTCOMPONENT_CXX
#define LARLITE_LATELIGHTCOMPONENT_CXX

#include "LateLightComponent.h"
#include "DataFormat/opflash.h"
#include "DataFormat/ophit.h"

namespace larlite {

  LateLightComponent::LateLightComponent()
    : _tree(nullptr)
  { 
    _name = "LateLightComponent"; 
    _fout = 0;
    _min_pe = 0.;
  }

  bool LateLightComponent::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("tree","_flash_tree");
    _tree->Branch("_n_flash",&_n_flash,"n_flash/I");
    _tree->Branch("_t_flash",&_t_flash,"t_flash/D");
    _tree->Branch("_pe_flash",&_pe_flash,"pe_flash/D");
    _tree->Branch("_q_tot",&_q_tot,"q_tot/D");
    _tree->Branch("_q_max",&_q_max,"q_max/D");
    _tree->Branch("_t_start",&_t_start,"t_start/D");
    _tree->Branch("_dt",&_dt,"dt/D");
    // vector of hit times
    _tree->Branch("_hit_t_v","std::vector<double>",&_hit_t_v);
    // vector of hit amplitudes
    _tree->Branch("_hit_q_v","std::vector<double>",&_hit_q_v);

    return true;
  }
  
  bool LateLightComponent::analyze(storage_manager* storage) {

    auto ev_opflash = storage->get_data<event_opflash>("opflash");
    auto ev_ophit   = storage->get_data<event_ophit>("opflash");
    
    _n_flash = ev_opflash->size();

    // find the hit with the largest charge -> this will be
    // the beginning of our single muon's flash
    // ignore any hits with time smaller than than
    _q_max   = 0.;
    _q_tot   = 0;
    for (size_t i=0; i < ev_ophit->size(); i++){
      auto h = ev_ophit->at(i);
      _q_tot += h.PE();
      if (h.PE() > _q_max){
	_t_start = h.PeakTime();
	_q_max   = h.PE();
      }
    }
    
    _hit_t_v.clear();
    _hit_q_v.clear();
    for (size_t i=0; i < ev_ophit->size(); i++){
      auto h = ev_ophit->at(i);
      // ignore if the hit-time is before the time
      // for the max charge hit
      double dt = h.PeakTime() - _t_start;
      if (dt < 0) continue;
      _hit_t_v.push_back(h.PeakTime() - _t_start);
      _hit_q_v.push_back(h.PE());
    }

    for (size_t i=0; i < ev_opflash->size(); i++){
      auto flash = ev_opflash->at(i);
      _t_flash  = flash.Time();
      _pe_flash = flash.TotalPE();
      // apply threshold to PE
      if (_pe_flash > _min_pe)
	_tree->Fill();
    }
      
    return true;
  }

  bool LateLightComponent::finalize() {

    if (_fout)
      _tree->Write();
    
    return true;
  }

}
#endif
