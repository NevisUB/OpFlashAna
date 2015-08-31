#ifndef LARLITE_MICHELRECO_CXX
#define LARLITE_MICHELRECO_CXX

#include "MichelReco.h"
#include "DataFormat/opflash.h"

namespace larlite {

  bool MichelReco::initialize() {

    //_PEmin = 0;
    _evt   = 0;
    
    // ttree for michels
    if (_michel_tree) delete _michel_tree;
    
    _michel_tree = new TTree("_michel_tree","michel tree");
    _michel_tree->Branch("_evt",&_evt,"evt/I");
    _michel_tree->Branch("_dt",&_dt,"dt/D");
    _michel_tree->Branch("_T",&_T,"T/D");
    _michel_tree->Branch("_Tpre",&_Tpre,"Tpre/D");
    _michel_tree->Branch("_PEpre",&_PEpre,"PEpre/D");
    _michel_tree->Branch("_PE",&_PE,"PE/D");

    return true;
  }
  
  bool MichelReco::analyze(storage_manager* storage) {

    auto ev_opflash = storage->get_data<event_opflash>("opflash");

    // loop through flashes in event
    // keep track of time-diff between flashes

    // also apply a cut on flashes

    // strategy: add flashes to a map
    // time -> flash id in event
    // this way we can go through the flashes
    // in time-order
    
    // clear map
    _flashMap.clear();

    
    for (size_t i=0; i < ev_opflash->size(); i++){
      
      auto& flash = ev_opflash->at(i);

      // if PE count above some limit add to map
      if (flash.TotalPE() > _PEmin)
	_flashMap[flash.Time()] = i;
    }

    if(!_flashMap.size()) return false;
    
    // vector holding indexes ordered by time
    std::vector<size_t> orderedFlashIndices;
    for (auto it = _flashMap.begin(); it != _flashMap.end(); it++)
      orderedFlashIndices.push_back(it->second);

    // loop through time ordered flashes
    for (size_t n=0; n < orderedFlashIndices.size()-1; n++){
      auto flashPre  = ev_opflash->at(orderedFlashIndices[n]);
      auto flash     = ev_opflash->at(orderedFlashIndices[n+1]);
      _T  = flash.Time()*0.5;
      _Tpre = flashPre.Time()*0.5;
      _dt = _T-_Tpre;
      _PE = flash.TotalPE();
      _PEpre = flashPre.TotalPE();
      
      if(_PE < _PEmin) {
	std::cout << "flash had a PE > _PEmin??" << std::endl;
	throw std::exception();
      }

_michel_tree->Fill();
    }
      
    _evt += 1;
    
    return true;
  }

  bool MichelReco::finalize() {

    if (_fout && _michel_tree)
      _michel_tree->Write();

    return true;
  }

}
#endif
