#ifndef OPT0FINDER_GEOMATCH_CXX
#define OPT0FINDER_GEOMATCH_CXX

#include "GeoMatch.h"
#include <cmath>
namespace flashana {

  GeoMatch::GeoMatch() : BaseFlashMatch()
  {
    _zdiff_max = 50*50;
    //_zdiff_max = 50*50;
  }

  FlashMatch_t GeoMatch::Match(const QCluster_t& pt_v, const Flash_t& flash)
  {
    // Compute the minimum z-position difference
    double min_dist = 1e12;
    ID_t min_id=kINVALID_ID;
    for(ID_t id=0; id<pt_v.size(); ++id) {
      auto const& pt = pt_v[id];
      double dist = pow(pt.z - flash.z,2);
      if(dist < min_dist) {
	min_dist = dist;
	min_id = id;
      }
    }

    // Prepare the return
    FlashMatch_t f;
    // If min-diff is bigger than assigned max, return default match (score<0)
    if(min_dist > _zdiff_max) return f;

    // Assign the score, return
    f.score = 1./min_dist;
    f.tpc_point = pt_v[min_id];

    return f;
  }
}
#endif
