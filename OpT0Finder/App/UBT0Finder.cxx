#ifndef LARLITE_UBT0FINDER_CXX
#define LARLITE_UBT0FINDER_CXX

#include "UBT0Finder.h"
#include "DataFormat/track.h"
#include "DataFormat/opflash.h"
#include "DataFormat/calorimetry.h"
#include "DataFormat/mctrack.h"
namespace larlite {

  bool UBT0Finder::initialize() {

    _tree = new TTree("flash_tree","");
    _tree->Branch("npe",&_npe,"npe/D");
    _tree->Branch("fy",&_flash_y,"fy/D");
    _tree->Branch("fz",&_flash_z,"fz/D");
    _tree->Branch("ty",&_tpc_y,"ty/D");
    _tree->Branch("tz",&_tpc_z,"tz/D");
    _tree->Branch("ft",&_flash_time,"ft/D");
    _tree->Branch("mct",&_mc_time,"mct/D");
    _tree->Branch("score",&_score,"score/D");
    
    return true;
  }
  
  bool UBT0Finder::analyze(storage_manager* storage) {

    _mgr.Reset();

    auto ev_flash = storage->get_data<event_opflash>("opflash");

    if(!ev_flash || ev_flash->empty()) return false;

    //auto ev_track = storage->get_data<event_track>("pandoraCosmicKHit");
    auto ev_track = storage->get_data<event_track>("trackkalmanhit");
    auto ev_mctrack = storage->get_data<event_mctrack>("mcreco");
    
    if(!_use_mc) {
      if (!ev_track || ev_track->empty()) return false;
      for(auto const& trk : *ev_track) {
	    
	::flashana::QCluster_t tpc_obj;
	tpc_obj.reserve(trk.NumberTrajectoryPoints()-1);
	for(size_t i=0; i < (trk.NumberTrajectoryPoints()-1); ++i) {
	  
	  auto const& pt1 = trk.LocationAtPoint(i);
	  auto const& pt2 = trk.LocationAtPoint(i+1);
	  
	  ::flashana::QPoint_t pt;
	  
	  double dx = pt2[0] - pt1[0];
	  double dy = pt2[1] - pt1[1];
	  double dz = pt2[2] - pt1[2];
	  
	  double dist = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
	  pt.q = (dist * 10. * 2.);
	  pt.x = pt1[0] + dx/2.;
	  pt.y = pt1[1] + dy/2.;
	  pt.z = pt1[2] + dz/2.;
	  
	  tpc_obj.emplace_back(pt);
	}
	_mgr.Emplace(std::move(tpc_obj));
      }
      
    }else{
      if (!ev_mctrack || ev_mctrack->empty()) return false;
      for(auto const& trk : *ev_mctrack) {

	//if(trk.size()<2) continue;
	
	::flashana::QCluster_t tpc_obj;

	if(trk.size()>=2) {
	  tpc_obj.reserve(trk.size()-1);
	  
	  for(size_t i=0; i < (trk.size()-1); ++i) {
	    
	    auto const& pt1 = trk[i].Position();
	    auto const& pt2 = trk[i+1].Position();
	    
	    ::flashana::QPoint_t pt;
	    
	    double dx = pt2[0] - pt1[0];
	    double dy = pt2[1] - pt1[1];
	    double dz = pt2[2] - pt1[2];
	    
	    pt.q = (trk[i].E() - trk[i+1].E());
	    pt.x = pt1[0] + dx/2.;
	    pt.y = pt1[1] + dy/2.;
	    pt.z = pt1[2] + dz/2.;
	    
	    tpc_obj.emplace_back(pt);
	  }
	}
	_mgr.Emplace(std::move(tpc_obj));
      }
    }
    
    for(auto const& flash : *ev_flash) {

      ::flashana::Flash_t f;
      f.x = f.x_err = 0;
      f.y = flash.YCenter();
      f.z = flash.ZCenter();
      f.y_err = flash.YWidth();
      f.z_err = flash.ZWidth();
      f.pe_v.reserve(32);
      for(unsigned int i=0; i<32; i++)
	f.pe_v.push_back(flash.PE(i));
      f.time = flash.Time();

      _mgr.Emplace(std::move(f));
    }

    /*
    for(auto const& calo : *ev_calo) {

      ::flashana::TPCObject_t tpc_obj;

      auto const& dedx = calo.dEdx();
      auto const& dx   = calo.TrkPitchVec();

      std::cout<<dedx.size()<<" : "<<dx.size()<<std::endl;

    }
    */

    auto const res = _mgr.Match();

    for(auto const& match : res) {
      auto const& flash = (*ev_flash)[match.flash_id];
      _flash_y = flash.YCenter();
      _flash_z = flash.ZCenter();
      _tpc_y   = match.tpc_point.y;
      _tpc_z   = match.tpc_point.z;
      _npe     = flash.TotalPE();
      _score   = match.score;
      _flash_time = flash.Time();
      if(_use_mc)
	_mc_time = (*ev_mctrack)[match.tpc_id][0].T() * 1.e-3;

      _tree->Fill();
    }
    return true;
  }

  bool UBT0Finder::finalize() {
    if(_fout) {
      _fout->cd();
      _tree->Write();
    }
    return true;
  }

}
#endif
