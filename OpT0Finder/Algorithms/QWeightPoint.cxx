#ifndef OPT0FINDER_QWEIGHTPOINT_CXX
#define OPT0FINDER_QWEIGHTPOINT_CXX

#include "QWeightPoint.h"
#include "Base/OpT0FinderException.h"
#include <cmath>
namespace flashana {

  QWeightPoint::QWeightPoint(const double x_step_size)
    : BaseFlashHypothesis()
    , _pos_x()
    , _pos_y()
    , _pos_z()
  {
    _x_step_size=x_step_size;
  }
    
  QWeightPoint::QWeightPoint(const std::vector<double>& pos_x,
			     const std::vector<double>& pos_y,
			     const std::vector<double>& pos_z,
			     const double x_step_size)
    : _pos_x(pos_x)
    , _pos_y(pos_y)
    , _pos_z(pos_z)
  {
    if(_pos_x.size() != _pos_y.size() || _pos_x.size() != _pos_z.size() )
      throw OpT0FinderException("Unmatching optical detector position array length!");
    _x_step_size = x_step_size;
  }

  QCluster_t QWeightPoint::FlashHypothesis(const QCluster_t& pt_v)
  {
    QCluster_t res;
    if(pt_v.empty()) return res;

    for(double x_step_size=1; x_step_size<250; x_step_size+=_x_step_size) {
      QPoint_t pt;

      pt.x = pt.y = pt.z = pt.q = 0;
      double weight_tot = 0;
      QPoint_t min_pt = pt_v[0];
      for(auto const& tpc_pt : pt_v)
	
	if(min_pt.x > tpc_pt.x) min_pt = tpc_pt;
      
      if(_pos_x.empty()) {
	
	for(auto const& tpc_pt : pt_v) {
	  
	  double weight = 0;
	  weight = tpc_pt.q / pow(tpc_pt.x - min_pt.x + x_step_size,2);
	  //pt.x += (tpc_pt.x - min_pt.x + x_step_size) * weight;
	  pt.y += (tpc_pt.y * weight);
	  pt.z += (tpc_pt.z * weight);
	  weight_tot += weight;
	}

	//pt.x /= weight_tot;
	pt.x = x_step_size;
	pt.y /= weight_tot;
	pt.z /= weight_tot;
	
      }else{
	
	for(auto const& tpc_pt : pt_v) {
	  
	  double weight = 0;
	  for(size_t i=0; i<_pos_x.size(); ++i) {
	    
	    double r2 = 0;
	    r2 += pow(tpc_pt.x - min_pt.x + x_step_size,2);
	    r2 += pow(tpc_pt.y - _pos_y[i],2);
	    r2 += pow(tpc_pt.z - _pos_z[i],2);
	    
	    weight += tpc_pt.q / r2;
	  }

	  pt.x += (tpc_pt.x - min_pt.x + x_step_size) * weight;
	  pt.y += (tpc_pt.y * weight);
	  pt.z += (tpc_pt.z * weight);
	  weight_tot += weight;
	}
	pt.x /= weight_tot;
	pt.y /= weight_tot;
	pt.z /= weight_tot;
      }
      res.push_back(pt);
    }
    return res;
  }


}
#endif
