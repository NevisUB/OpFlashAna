/**
 * \file GeoMatch.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class GeoMatch
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_GEOMATCH_H
#define OPT0FINDER_GEOMATCH_H

#include "Base/BaseFlashMatch.h"

namespace flashana {
  /**
     \class GeoMatch
     Implementation of BaseFlashMatch algorithm class. The concept is   \n
     to use a simple geometry-based matching. Currently only implements \n
     a differencei n z-position. The best match is defined by the shortest \n
     difference in z-position, and the distance-squared is assigned as a match score. 
  */
  class GeoMatch : public BaseFlashMatch {
    
  public:
    
    /// Default constructor
    GeoMatch();
    
    /// Default destructor
    ~GeoMatch(){}

    FlashMatch_t Match(const QCluster_t&, const Flash_t&);

    void Configure(double zdiff_max) { _zdiff_max = zdiff_max; }

  private:

    double _zdiff_max; ///< allowed diff in z-direction to be considered as a match
    
  };
}

#endif
/** @} */ // end of doxygen group 

