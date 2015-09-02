/**
 * \file BaseFlashHypothesis.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseFlashHypothesis
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef OPT0FINDER_BASEFLASHHYPOTHESIS_H
#define OPT0FINDER_BASEFLASHHYPOTHESIS_H

#include "BaseAlgorithm.h"

namespace flashana {
  /**
     \class BaseFlashHypothesis
     Algorithm base class for estimating candidate flash points based on TPC information    \n
     information. The relevant function is FlashHypothesis which takes flashana::QCluster_t \n
     (3D TPC object representation) and returns another QCluster_t that is a list of possible \n
     flash point to be used for matching.\n
  */
  class BaseFlashHypothesis : public BaseAlgorithm{
    
  public:
    
    /// Default constructor
    BaseFlashHypothesis() : BaseAlgorithm(kFlashHypothesis)
    {}
    
    /// Default destructor
    virtual ~BaseFlashHypothesis(){}

    /**
       CORE FUNCTION: takes in TPC object information in the form of QCluster_t, \n
       and returns a list of flash hypothesis 3D points in QCluster_t representation.
     */
    virtual QCluster_t FlashHypothesis(const QCluster_t&) = 0;
    
  };
}

#endif
/** @} */ // end of doxygen group 

