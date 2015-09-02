#ifndef OPT0FINDER_OPT0FINDERTYPES_H
#define OPT0FINDER_OPT0FINDERTYPES_H

#include <vector>
#include "OpT0FinderConstants.h"

namespace flashana {

  /// Enumerator for different types of algorithm
  enum Algorithm_t {
    kTPCFilter,       ///< Algorithm type to filter out TPC objects from matching candidate list
    kFlashFilter,     ///< Algorithm type to filter out flash from matching candidate list
    kFlashHypothesis, ///< Algorithm type to make flash hypothesis from TPC object info
    kFlashMatch,      ///< Algorithm type to match flash hypothesis and reconstructed flash
    kAlgorithmTypeMax ///< enum flag for algorithm type count & invalid type
  };

  /// Struct to represent an optical flash
  struct Flash_t {
  public:

    std::vector<double> pe_v; ///< PE distribution over photo-detectors
    double x,y,z;             ///< Flash position 
    double x_err,y_err,z_err; ///< Flash timing, a candidate T0
    double time;
    /// Default ctor assigns invalid values
    Flash_t() {
      x = y = z = kINVALID_DOUBLE;
      x_err = y_err = z_err = kINVALID_DOUBLE;
      time = kINVALID_DOUBLE;
    }
  };

  /// Struct to represent an energy deposition point in 3D space
  struct QPoint_t{

    double x,y,z; ///< Spatial position in [cm]
    double q;     ///< Charge in an arbitrary unit
    /// Default ctor assigns invalid values
    QPoint_t()
      : x(kINVALID_DOUBLE)
      , y(kINVALID_DOUBLE)
      , z(kINVALID_DOUBLE)
      , q(kINVALID_DOUBLE)
    {}
    /// Alternative ctor
    QPoint_t(double xvalue,
	     double yvalue,
	     double zvalue,
	     double qvalue)
      : x(xvalue)
      , y(yvalue)
      , z(zvalue)
      , q(qvalue)
    {}
  };

  /// Collection of charge deposition 3D point (cluster)
  typedef std::vector<flashana::QPoint_t> QCluster_t;
  /// Collection of 3D point clusters (one use case is TPC object representation for track(s) and shower(s))
  typedef std::vector<flashana::QCluster_t> QClusterArray_t;
  /// Collection of Flash objects
  typedef std::vector<flashana::Flash_t> FlashArray_t;
  /// Index used to identify Flash_t/QPointCollection_t uniquely in an event
  typedef size_t ID_t;
  /// Index collection
  typedef std::vector<flashana::ID_t> IDArray_t;

  /// Invalid ID
  const ID_t kINVALID_ID = kINVALID_SIZE;

  /// Flash-TPC match info
  struct FlashMatch_t {
    ID_t tpc_id;   ///< matched TPC object ID
    ID_t flash_id; ///< matched Flash ID
    double score;  ///< floating point representing the "goodness" (algorithm dependent) 
    QPoint_t tpc_point; ///< estimated & matched 3D flash hypothesis point from TPC information
    /// Default ctor assigns invalid values
    FlashMatch_t()
    { tpc_id = kINVALID_ID; flash_id = kINVALID_ID; score = -1; }
    /// Alternative ctor
    FlashMatch_t(const ID_t& tpc_id_value,
		 const ID_t& flash_id_value,
		 const double& score_value)
    { tpc_id = tpc_id_value; flash_id = flash_id_value; score = score_value; }
  };
}
#endif
