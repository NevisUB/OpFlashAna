/**
 * \file FlatTreeWriter.h
 *
 * \ingroup OpMichel
 * 
 * \brief Class def header for a class FlatTreeWriter
 *
 * @author vgenty
 */

/** \addtogroup OpMichel

    @{*/

#ifndef LARLITE_FLATTREEWRITER_H
#define LARLITE_FLATTREEWRITER_H

#include <vector>

#include "Analysis/ana_base.h"

#include "DataFormat/opflash.h"

#include "TTree.h"


namespace larlite {

  class FlatTreeWriter : public ana_base{
  
  public:

    FlatTreeWriter()
    { _name="FlatTreeWriter"; _fout=0; _flash_tree = nullptr;}
    virtual ~FlatTreeWriter(){}
    virtual bool initialize();
    virtual bool analyze(storage_manager* storage);
    virtual bool finalize();

  protected:

    TTree *_flash_tree;
    
    int _evt;
    int _flashnum;
    
    double _Time;
    double _TimeWidth;
    double _AbsTime;

    int _Frame;

    double _YCenter;
    double _YWidth;
    double _ZCenter;
    double _ZWidth;

    double _FastToTotal;
   
    std::vector<double> _WireCenters;
    std::vector<double> _WireWidths;
    
    bool _InBeamFrame;
    int  _OnBeamTime;
    int  _TotalPE;

    void clear_all();
    
  };
}
#endif
