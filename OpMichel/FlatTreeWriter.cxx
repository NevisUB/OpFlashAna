#ifndef LARLITE_FLATTREEWRITER_CXX
#define LARLITE_FLATTREEWRITER_CXX

#include "FlatTreeWriter.h"

namespace larlite {

  bool FlatTreeWriter::initialize() {
    // ttree for michels
    
    if (_flash_tree) delete _flash_tree;
    
    _flash_tree = new TTree("_flash_tree","flash tree");

    _flash_tree->Branch("_evt"     ,&_evt,"evt/I");
    _flash_tree->Branch("_flashnum",&_flashnum,"flashnum/I");

    _flash_tree->Branch("_Time",      &_Time,      "Time/D");
    _flash_tree->Branch("_TimeWidth", &_TimeWidth, "TimeWidth/D");
    _flash_tree->Branch("_AbsTime",   &_AbsTime,   "AbsTime/D");

    _flash_tree->Branch("_Frame", &_Frame, "Frame/I");

    _flash_tree->Branch("_YCenter", &_YCenter, "YCenter/D");
    _flash_tree->Branch("_YWidth",  &_YWidth,  "YWidth/D" );
    _flash_tree->Branch("_ZCenter", &_ZCenter, "ZCenter/D");
    _flash_tree->Branch("_ZWidth",  &_ZWidth,  "ZWidth/D" );

    _flash_tree->Branch("_FastToTotal", &_FastToTotal, "FastToTotal/D");

    _flash_tree->Branch("_WireCenters", &_WireCenters);
    _flash_tree->Branch("_WireWidths", &_WireWidths);

    _flash_tree->Branch("_InBeamFrame", &_InBeamFrame, "InBeamFrame/O");
    _flash_tree->Branch("_OnBeamTime",  &_OnBeamTime, "OnBeamTime/I");
    _flash_tree->Branch("_TotalPE",     &_TotalPE, "TotalPE/I");
    


    return true;
  }
  
  bool FlatTreeWriter::analyze(storage_manager* storage) {

    clear_all();

    auto ev_opflash = storage->get_data<event_opflash>("opflash");
    
    _evt = storage->event_id();
    
    for (size_t i = 0; i < ev_opflash->size(); i++){
      auto& flash = ev_opflash->at(i);
      
      _flashnum = (int)i;
      
      _Time      = flash.Time();
      _TimeWidth = flash.TimeWidth();
      _AbsTime   = flash.AbsTime();

      _Frame     = flash.Frame();

      _YCenter   = flash.YCenter();
      _YWidth    = flash.YWidth();
      _ZCenter   = flash.ZCenter();
      _ZWidth    = flash.ZWidth();
      
      _FastToTotal = flash.FastToTotal();
      
      _InBeamFrame = flash.InBeamFrame();
      _OnBeamTime  = flash.OnBeamTime();
      _TotalPE     = flash.TotalPE();
   

      _WireCenters = flash.WireCenters();
      _WireWidths  = flash.WireWidths() ;

      _flash_tree->Fill();
    }
    
    return true;
  }

  bool FlatTreeWriter::finalize() {
    
    if (_fout && _flash_tree)
      _flash_tree->Write();
    
    return true;
  }

  void FlatTreeWriter::clear_all(){
    _evt = -1;
    
    _flashnum = -1;
	 
    _Time      = -1;
    _TimeWidth = -1;
    _AbsTime   = -1;
    _Frame     = -1;
    _YCenter   = -1;
    _YWidth    = -1;
    _ZCenter   = -1;
    _ZWidth    = -1;
    _FastToTotal = -1;
    _InBeamFrame = false;
    _OnBeamTime  = -1;
    _TotalPE     = -1;

    _WireCenters.clear();
    _WireWidths.clear();
  }
  
}
#endif
