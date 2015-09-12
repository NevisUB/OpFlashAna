#ifndef QLLMATCH_CXX
#define QLLMATCH_CXX

#include "QLLMatch.h"
#include "Base/OpT0FinderException.h"
#include <TMinuit.h>
namespace flashana {

  QLLMatch* QLLMatch::_me = nullptr;

  void  MIN_vtx_qll (Int_t &, Double_t *, Double_t &, Double_t *, Int_t);

  QLLMatch::QLLMatch()
    : _pmt_x_v()
    , _pmt_y_v()
    , _pmt_z_v()
    , _minuit_ptr(nullptr)
  {}

  void QLLMatch::SetOpDetPositions( const std::vector<double>& pos_x,
				    const std::vector<double>& pos_y,
				    const std::vector<double>& pos_z )
  {
    if(pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() )
      throw OpT0FinderException("Unmatching optical detector position array length!");
    _pmt_x_v = pos_x;
    _pmt_y_v = pos_y;
    _pmt_z_v = pos_z;
  }

  FlashMatch_t QLLMatch::Match(const QCluster_t& pt_v, const Flash_t& flash)
  {
    if(_qll_hypothesis_v.empty()) {
      if(_pmt_x_v.empty()) throw OpT0FinderException("No optical detector found!");
      _qll_hypothesis_v.resize(_pmt_x_v.size(),0);
      _flash_pe_v.resize(_pmt_x_v.size(),0);
    }else if(_qll_hypothesis_v.size() != _pmt_x_v.size())
      throw OpT0FinderException("PMT dimension has changed!");

    if(flash.pe_v.size() != _qll_hypothesis_v.size())
      throw OpT0FinderException("Flash contians different # of op-det charge info than provided geo!");

    for(auto& v : _qll_hypothesis_v) v=0;
    
    _trk_x_v.resize(pt_v.size());
    _trk_y_v.resize(pt_v.size());
    _trk_z_v.resize(pt_v.size());
    _trk_q_v.resize(pt_v.size());
    double min_x = 1e9;
    for(size_t i=0; i<pt_v.size(); ++i) {
      auto const& pt = pt_v[i];
      _trk_x_v[i] = pt.x;
      _trk_y_v[i] = pt.y;
      _trk_z_v[i] = pt.z;
      _trk_q_v[i] = pt.q;
      if(pt.x < min_x) min_x = pt.x;
    }

    for(auto& v : _trk_x_v) v -= min_x;

    double max_pe = 0;
    for(auto const& v : flash.pe_v) if(v > max_pe) max_pe = v;

    for(size_t i=0; i<flash.pe_v.size(); ++i)

      _flash_pe_v[i] = flash.pe_v[i] / max_pe;

    this->CallMinuit();

    // Estimate position
    FlashMatch_t res;
    if(isnan(_qll)) return res;
    res.tpc_point.x = res.tpc_point.y = res.tpc_point.z = 0;
    res.score = 1./_qll;
    double weight_sum = 0;
    for(auto const& pt : pt_v) {

      double weight = 0;

      for(size_t pmt_index = 0; pmt_index < _pmt_x_v.size(); ++pmt_index) {

	double r2 = ( pow( _pmt_x_v[pmt_index] - (pt.x - min_x + _reco_x), 2) +
		      pow( _pmt_y_v[pmt_index] - pt.y, 2) +
		      pow( _pmt_z_v[pmt_index] - pt.z, 2) );

	weight += pt.q / r2;
      }

      res.tpc_point.x += (pt.x - min_x + _reco_x) * weight;
      res.tpc_point.y += pt.y * weight;
      res.tpc_point.z += pt.z * weight;

      weight_sum += weight;
    }
    res.tpc_point.x /= weight_sum;
    res.tpc_point.y /= weight_sum;
    res.tpc_point.z /= weight_sum;

    return res;
  }

  double QLLMatch::ChargeHypothesis(double x)
  {
    double qmax = 0;
    for(size_t pmt_index=0; pmt_index < _pmt_x_v.size(); ++pmt_index) {

      for(size_t pt_index=0; pt_index < _trk_x_v.size(); ++pt_index) {

	double r2 = ( pow(_pmt_x_v[pmt_index] - (_trk_x_v[pt_index] + x),2) +
		      pow(_pmt_y_v[pmt_index] - _trk_y_v[pt_index], 2) +
		      pow(_pmt_z_v[pmt_index] - _trk_z_v[pt_index], 2) );

	double angle = (_pmt_x_v[pmt_index] - (_trk_x_v[pt_index] + x)) /
	  sqrt( pow(_pmt_y_v[pmt_index] - _trk_y_v[pt_index],2) +
		pow(_pmt_z_v[pmt_index] - _trk_y_v[pt_index],2) );

	if(angle<0) angle *= -1;
	
	_qll_hypothesis_v[pmt_index] += _trk_q_v[pt_index] * angle / r2;
      }

      if(_qll_hypothesis_v[pmt_index] > qmax) qmax = _qll_hypothesis_v[pmt_index];
    }

    for(auto& v : _qll_hypothesis_v) v /= qmax;

    double result = 0;
    double nvalid_pmt = 0;
    for(size_t pmt_index=0; pmt_index < _qll_hypothesis_v.size(); ++pmt_index) {

      if( _flash_pe_v[pmt_index] < 0.01 ) continue;

      nvalid_pmt += 1;
      
      result += std::fabs( _qll_hypothesis_v[pmt_index] - _flash_pe_v[pmt_index] ) * _flash_pe_v[pmt_index];

    }
    return result / nvalid_pmt;
    
  }
  
  void MIN_vtx_qll(Int_t & /*Npar*/,
		   Double_t * /*Grad*/,
		   Double_t & Fval,
		   Double_t * Xval,
		   Int_t /*Flag*/){

    Fval = QLLMatch::GetME().ChargeHypothesis(Xval[0]);

    return;
  }  


  void QLLMatch::CallMinuit() {

    if(!_minuit_ptr) _minuit_ptr = new TMinuit(4);

    double  reco_x = 128.;
    double  reco_x_err = 0;
    double  MinFval;
    int     ierrflag,npari,nparx,istat;
    double  arglist[2],error[4],Fmin,Fedm,Errdef;
    ierrflag = npari = nparx = istat = 0;

    _minuit_ptr->SetPrintLevel(-1);
    arglist[0] = 1.0;
    _minuit_ptr->mnexcm("SET STR",arglist,1,ierrflag);

    _minuit_ptr->SetFCN(MIN_vtx_qll);

    _minuit_ptr->DefineParameter (0,"X", reco_x, 10, 0.0, 256.);

    _minuit_ptr->Command("SET NOW");

    arglist[0]   = 5.0e+2;
    arglist[1]   = 1.0e-6;
    _minuit_ptr->mnexcm ("simplex",arglist,2,ierrflag);

    _minuit_ptr->GetParameter (0,reco_x,reco_x_err);

    _minuit_ptr->mnstat (Fmin,Fedm,Errdef,npari,nparx,istat);
    MinFval = Fmin;
    double *grad=0;
    int nPar=1;
    double fValue[1];
    fValue[0]=reco_x;
    // Transfer the minimization variables:
    MIN_vtx_qll(nPar, grad, Fmin,
		fValue, ierrflag);
    static bool show = true;
    /*
    if(show){
      if(Fmin!=MinFval)std::cout<<"Fmin "<<Fmin<<" not equall to "<<MinFval<<std::endl;
      show=false;
    }
    */

    // Transfer the minimization variables:
    _reco_x = reco_x;
    _reco_x_err = reco_x_err;
    _qll    = MinFval;

    // Clear:
    _minuit_ptr->mnexcm ("clear",arglist,0,ierrflag);

    if (_minuit_ptr) delete _minuit_ptr;
    _minuit_ptr = 0;

    return;
  }

}
#endif
