#ifndef LARLITE_MICHELRECO_OPDIGIT_CXX
#define LARLITE_MICHELRECO_OPDIGIT_CXX

#include "MichelReco_OpDigit.h"
#include "DataFormat/opdetwaveform.h"
#include "DataFormat/mcshower.h"

namespace larlite {

  MichelReco_OpDigit::MichelReco_OpDigit()
    : _tree(nullptr)
  {
    _name = "MichelReco_OpDigit";
    _fout = 0;
    _use_mc = false;
  }


  bool MichelReco_OpDigit::initialize() {

   _baseline = 2048;
   
   fillTree();

    return true;
  }
  
  bool MichelReco_OpDigit::analyze(storage_manager* storage) {

    cleanTree();

    if (_use_mc){

      std::cout << "using mc showers" << std::endl;
      auto ev_mcshowers = storage->get_data<event_mcshower>("mcreco");

      for (size_t n=0; n < ev_mcshowers->size(); n++){

	auto const& shr = ev_mcshowers->at(n);

	if (shr.Process() == "muMinusCaptureAtRest"){
	  std::cout << "found a michel" << std::endl;
	  _michel = 1;
	  _E = shr.DetProfile().E();
	  _t = shr.DetProfile().T();
	  _x = shr.DetProfile().X();
	  _y = shr.DetProfile().Y();
	  _z = shr.DetProfile().Z();
	}// if a michel
      }// for all mcshowers
    }// if we use mc

    auto ev_opdetwaveform = storage->get_data<event_opdetwaveform>("pmtreadout");

    if ( (!ev_opdetwaveform) or (ev_opdetwaveform->size() == 0) ){
      std::cout << "No waveforms found..." << std::endl;
      return false;
    }

    // for each waveform, find all the times at which
    // the waveform reaches a maximum, keep track of all these maxima
    std::vector<std::vector<short> > ADCmax(32, std::vector<short>() );
    std::vector<std::vector<short> > TDCmax(32, std::vector<short>() );


    // first find the peak-time for the muons
    std::vector<short> peakADC(32,0.);
    std::vector<size_t> peakTDC(32,0.);
    for (size_t i=0; i < ev_opdetwaveform->size(); i++){
      auto const& wf = ev_opdetwaveform->at(i);
      short pmt = wf.ChannelNumber();
      // if too few samples -> continue
      if (wf.size() < 800)
	continue;
      if (pmt >= 32)
	continue;
      findMuonPeak(wf,peakTDC[pmt],peakADC[pmt]);
    }
    // find the largest peak:
    size_t peakTime = 0;
    short maxADC = 0;
    for (size_t i=0; i < 32; i++){
      if (peakADC[i] > maxADC)
	peakTime = peakTDC[i];
    }

    std::cout << "peak time is : " << peakTime << std::endl;
    

    for (size_t i=0; i < ev_opdetwaveform->size(); i++){

      auto const& wf = ev_opdetwaveform->at(i);

      short pmt = wf.ChannelNumber();

      // if too few samples -> continue
      if (wf.size() < 800)
	continue;

      if (pmt >= 32)
	continue;

      //int max = findMaxima(wf,0,1000,ADCmax[0],TDCmax[0]);
      //std::cout << "PMT " << pmt << " of length " << wf.size() << " has max " << max-2048 << std::endl;      
      
      // loop trhough the wf in search for all the local maxima
      int thismax = findMaxima(wf,0,1000,ADCmax[pmt],TDCmax[pmt]);
      // add info to tree
      add(pmt,ADCmax[pmt],TDCmax[pmt],wf);

    }// for all PMTs

    _tree->Fill();
      
    return true;
  }

  bool MichelReco_OpDigit::finalize() {

    if (_fout){
      if (_tree) _tree->Write();
    }

    return true;
  }

  // find the peak time of the muon
  void MichelReco_OpDigit::findMuonPeak(const std::vector<short>& wf,
					size_t& max_tdc,
					short& max_adc)
  {

    max_tdc = 0;
    max_adc = 0;
    
    for (size_t i=0; i < wf.size(); i++){
      if (wf[i] > max_adc){
	max_tdc = i;
	max_adc = wf[i];
      }
    }
    
    return;
  }

  int MichelReco_OpDigit::findMaxima(const std::vector<short>& wf,
				     const size_t& tick_min,
				     const size_t& tick_max,
				     std::vector<short>& ADCmax,
				     std::vector<short>& TDCmax)
  {

    ADCmax.clear();
    TDCmax.clear();

    int maxmax = 0;

    short currentMaxADC = 0;
    short currentMaxTDC = 0;

    for (size_t i = tick_min; i < tick_max; i++){

      // if out of bounds
      if ( i > wf.size() )
	break;

      if (wf[i] > maxmax)
	maxmax = wf[i];

      //std::cout << "ADC @ tick " << i << " : " << wf[i] << std::endl;
      
      // if we are at a new local maximum
      if (wf[i] > currentMaxADC){
	currentMaxADC = wf[i];
	currentMaxTDC = i;
      }

      if (wf[i] < currentMaxADC){
	// if we previously were at a maximum
	if ( (currentMaxTDC+1) == i){
	  if ( wf[i]-_baseline > 25){
	    //std::cout << "found new maxima!" << std::endl;
	    std::cout << "new max : " << wf[i]-_baseline << std::endl;
	    ADCmax.push_back(wf[i]-_baseline);
	    TDCmax.push_back(i);
	  }
	}
	// if not -> re-set the max ADC info
	else
	  currentMaxADC = 0;
      }// if we are below the maximum
    }// for all TDCs
	
    return maxmax;
  }

  void MichelReco_OpDigit::fillTree(){

    if (_tree) delete _tree;
    
    _tree = new TTree("_tree","michel tree");
    _tree->Branch("_event",&_event,"event/I");
    _tree->Branch("_michel",&_michel,"michel/I");
    _tree->Branch("_E",&_E,"E/D");
    _tree->Branch("_t",&_t,"t/D");
    _tree->Branch("_x",&_x,"x/D");
    _tree->Branch("_y",&_y,"y/D");
    _tree->Branch("_z",&_z,"z/D");
    _tree->Branch("tdc00","std::vector<short>",&_tdc00);
    _tree->Branch("tdc01","std::vector<short>",&_tdc01);
    _tree->Branch("tdc02","std::vector<short>",&_tdc02);
    _tree->Branch("tdc03","std::vector<short>",&_tdc03);
    _tree->Branch("tdc04","std::vector<short>",&_tdc04);
    _tree->Branch("tdc05","std::vector<short>",&_tdc05);
    _tree->Branch("tdc06","std::vector<short>",&_tdc06);
    _tree->Branch("tdc07","std::vector<short>",&_tdc07);
    _tree->Branch("tdc08","std::vector<short>",&_tdc08);
    _tree->Branch("tdc09","std::vector<short>",&_tdc09);
    _tree->Branch("tdc10","std::vector<short>",&_tdc10);
    _tree->Branch("tdc11","std::vector<short>",&_tdc11);
    _tree->Branch("tdc12","std::vector<short>",&_tdc12);
    _tree->Branch("tdc13","std::vector<short>",&_tdc13);
    _tree->Branch("tdc14","std::vector<short>",&_tdc14);
    _tree->Branch("tdc15","std::vector<short>",&_tdc15);
    _tree->Branch("tdc16","std::vector<short>",&_tdc16);
    _tree->Branch("tdc17","std::vector<short>",&_tdc17);
    _tree->Branch("tdc18","std::vector<short>",&_tdc18);
    _tree->Branch("tdc19","std::vector<short>",&_tdc19);
    _tree->Branch("tdc20","std::vector<short>",&_tdc20);
    _tree->Branch("tdc21","std::vector<short>",&_tdc21);
    _tree->Branch("tdc22","std::vector<short>",&_tdc22);
    _tree->Branch("tdc23","std::vector<short>",&_tdc23);
    _tree->Branch("tdc24","std::vector<short>",&_tdc24);
    _tree->Branch("tdc25","std::vector<short>",&_tdc25);
    _tree->Branch("tdc26","std::vector<short>",&_tdc26);
    _tree->Branch("tdc27","std::vector<short>",&_tdc27);
    _tree->Branch("tdc28","std::vector<short>",&_tdc28);
    _tree->Branch("tdc29","std::vector<short>",&_tdc29);
    _tree->Branch("tdc30","std::vector<short>",&_tdc30);
    _tree->Branch("tdc31","std::vector<short>",&_tdc31);
    _tree->Branch("adc00","std::vector<short>",&_adc00);
    _tree->Branch("adc01","std::vector<short>",&_adc01);
    _tree->Branch("adc02","std::vector<short>",&_adc02);
    _tree->Branch("adc03","std::vector<short>",&_adc03);
    _tree->Branch("adc04","std::vector<short>",&_adc04);
    _tree->Branch("adc05","std::vector<short>",&_adc05);
    _tree->Branch("adc06","std::vector<short>",&_adc06);
    _tree->Branch("adc07","std::vector<short>",&_adc07);
    _tree->Branch("adc08","std::vector<short>",&_adc08);
    _tree->Branch("adc09","std::vector<short>",&_adc09);
    _tree->Branch("adc10","std::vector<short>",&_adc10);
    _tree->Branch("adc11","std::vector<short>",&_adc11);
    _tree->Branch("adc12","std::vector<short>",&_adc12);
    _tree->Branch("adc13","std::vector<short>",&_adc13);
    _tree->Branch("adc14","std::vector<short>",&_adc14);
    _tree->Branch("adc15","std::vector<short>",&_adc15);
    _tree->Branch("adc16","std::vector<short>",&_adc16);
    _tree->Branch("adc17","std::vector<short>",&_adc17);
    _tree->Branch("adc18","std::vector<short>",&_adc18);
    _tree->Branch("adc19","std::vector<short>",&_adc19);
    _tree->Branch("adc20","std::vector<short>",&_adc20);
    _tree->Branch("adc21","std::vector<short>",&_adc21);
    _tree->Branch("adc22","std::vector<short>",&_adc22);
    _tree->Branch("adc23","std::vector<short>",&_adc23);
    _tree->Branch("adc24","std::vector<short>",&_adc24);
    _tree->Branch("adc25","std::vector<short>",&_adc25);
    _tree->Branch("adc26","std::vector<short>",&_adc26);
    _tree->Branch("adc27","std::vector<short>",&_adc27);
    _tree->Branch("adc28","std::vector<short>",&_adc28);
    _tree->Branch("adc29","std::vector<short>",&_adc29);
    _tree->Branch("adc30","std::vector<short>",&_adc30);
    _tree->Branch("adc31","std::vector<short>",&_adc31);
    _tree->Branch("wf00","std::vector<short>",&_wf00);
    _tree->Branch("wf01","std::vector<short>",&_wf01);
    _tree->Branch("wf02","std::vector<short>",&_wf02);
    _tree->Branch("wf03","std::vector<short>",&_wf03);
    _tree->Branch("wf04","std::vector<short>",&_wf04);
    _tree->Branch("wf05","std::vector<short>",&_wf05);
    _tree->Branch("wf06","std::vector<short>",&_wf06);
    _tree->Branch("wf07","std::vector<short>",&_wf07);
    _tree->Branch("wf08","std::vector<short>",&_wf08);
    _tree->Branch("wf09","std::vector<short>",&_wf09);
    _tree->Branch("wf10","std::vector<short>",&_wf10);
    _tree->Branch("wf11","std::vector<short>",&_wf11);
    _tree->Branch("wf12","std::vector<short>",&_wf12);
    _tree->Branch("wf13","std::vector<short>",&_wf13);
    _tree->Branch("wf14","std::vector<short>",&_wf14);
    _tree->Branch("wf15","std::vector<short>",&_wf15);
    _tree->Branch("wf16","std::vector<short>",&_wf16);
    _tree->Branch("wf17","std::vector<short>",&_wf17);
    _tree->Branch("wf18","std::vector<short>",&_wf18);
    _tree->Branch("wf19","std::vector<short>",&_wf19);
    _tree->Branch("wf20","std::vector<short>",&_wf20);
    _tree->Branch("wf21","std::vector<short>",&_wf21);
    _tree->Branch("wf22","std::vector<short>",&_wf22);
    _tree->Branch("wf23","std::vector<short>",&_wf23);
    _tree->Branch("wf24","std::vector<short>",&_wf24);
    _tree->Branch("wf25","std::vector<short>",&_wf25);
    _tree->Branch("wf26","std::vector<short>",&_wf26);
    _tree->Branch("wf27","std::vector<short>",&_wf27);
    _tree->Branch("wf28","std::vector<short>",&_wf28);
    _tree->Branch("wf29","std::vector<short>",&_wf29);
    _tree->Branch("wf30","std::vector<short>",&_wf30);
    _tree->Branch("wf31","std::vector<short>",&_wf31);
    
    return;
  }

  void MichelReco_OpDigit::add(const short& pmt,
			       const std::vector<short>& ADCmax,
			       const std::vector<short>& TDCmax,
			       const std::vector<short>& wf)
  {

    if ( (pmt >= 32) or (pmt < 0) )
      return;

    if (pmt ==  0){_tdc00 = TDCmax; _adc00 = ADCmax; _wf00 = wf; }
    if (pmt ==  1){_tdc01 = TDCmax; _adc01 = ADCmax; _wf01 = wf; }
    if (pmt ==  2){_tdc02 = TDCmax; _adc02 = ADCmax; _wf02 = wf; }
    if (pmt ==  3){_tdc03 = TDCmax; _adc03 = ADCmax; _wf03 = wf; }
    if (pmt ==  4){_tdc04 = TDCmax; _adc04 = ADCmax; _wf04 = wf; }
    if (pmt ==  5){_tdc05 = TDCmax; _adc05 = ADCmax; _wf05 = wf; }
    if (pmt ==  6){_tdc06 = TDCmax; _adc06 = ADCmax; _wf06 = wf; }
    if (pmt ==  7){_tdc07 = TDCmax; _adc07 = ADCmax; _wf07 = wf; }
    if (pmt ==  8){_tdc08 = TDCmax; _adc08 = ADCmax; _wf08 = wf; }
    if (pmt ==  9){_tdc09 = TDCmax; _adc09 = ADCmax; _wf09 = wf; }
    if (pmt == 10){_tdc10 = TDCmax; _adc10 = ADCmax; _wf10 = wf; }
    if (pmt == 11){_tdc11 = TDCmax; _adc11 = ADCmax; _wf11 = wf; }
    if (pmt == 12){_tdc12 = TDCmax; _adc12 = ADCmax; _wf12 = wf; }
    if (pmt == 13){_tdc13 = TDCmax; _adc13 = ADCmax; _wf13 = wf; }
    if (pmt == 14){_tdc14 = TDCmax; _adc14 = ADCmax; _wf14 = wf; }
    if (pmt == 15){_tdc15 = TDCmax; _adc15 = ADCmax; _wf15 = wf; }
    if (pmt == 16){_tdc16 = TDCmax; _adc16 = ADCmax; _wf16 = wf; }
    if (pmt == 17){_tdc17 = TDCmax; _adc17 = ADCmax; _wf17 = wf; }
    if (pmt == 18){_tdc18 = TDCmax; _adc18 = ADCmax; _wf18 = wf; }
    if (pmt == 19){_tdc19 = TDCmax; _adc19 = ADCmax; _wf19 = wf; }
    if (pmt == 20){_tdc20 = TDCmax; _adc20 = ADCmax; _wf20 = wf; }
    if (pmt == 21){_tdc21 = TDCmax; _adc21 = ADCmax; _wf21 = wf; }
    if (pmt == 22){_tdc22 = TDCmax; _adc22 = ADCmax; _wf22 = wf; }
    if (pmt == 23){_tdc23 = TDCmax; _adc23 = ADCmax; _wf23 = wf; }
    if (pmt == 24){_tdc24 = TDCmax; _adc24 = ADCmax; _wf24 = wf; }
    if (pmt == 25){_tdc25 = TDCmax; _adc25 = ADCmax; _wf25 = wf; }
    if (pmt == 26){_tdc26 = TDCmax; _adc26 = ADCmax; _wf26 = wf; }
    if (pmt == 27){_tdc27 = TDCmax; _adc27 = ADCmax; _wf27 = wf; }
    if (pmt == 28){_tdc28 = TDCmax; _adc28 = ADCmax; _wf28 = wf; }
    if (pmt == 29){_tdc29 = TDCmax; _adc29 = ADCmax; _wf29 = wf; }
    if (pmt == 30){_tdc30 = TDCmax; _adc30 = ADCmax; _wf30 = wf; }
    if (pmt == 31){_tdc31 = TDCmax; _adc31 = ADCmax; _wf31 = wf; }

    return;
  }

  void MichelReco_OpDigit::cleanTree(){

    _E = 0;
    _t = 0;
    _x = _y = _z = 0;

    _michel = 0;
    
    _tdc00.clear();
    _tdc01.clear();
    _tdc02.clear();
    _tdc03.clear();
    _tdc04.clear();
    _tdc05.clear();
    _tdc06.clear();
    _tdc07.clear();
    _tdc08.clear();
    _tdc09.clear();
    _tdc10.clear();
    _tdc11.clear();
    _tdc12.clear();
    _tdc13.clear();
    _tdc14.clear();
    _tdc15.clear();
    _tdc16.clear();
    _tdc17.clear();
    _tdc18.clear();
    _tdc19.clear();
    _tdc20.clear();
    _tdc21.clear();
    _tdc22.clear();
    _tdc23.clear();
    _tdc24.clear();
    _tdc25.clear();
    _tdc26.clear();
    _tdc27.clear();
    _tdc28.clear();
    _tdc29.clear();
    _tdc30.clear();
    _tdc31.clear();
    _adc00.clear();
    _adc01.clear();
    _adc02.clear();
    _adc03.clear();
    _adc04.clear();
    _adc05.clear();
    _adc06.clear();
    _adc07.clear();
    _adc08.clear();
    _adc09.clear();
    _adc10.clear();
    _adc11.clear();
    _adc12.clear();
    _adc13.clear();
    _adc14.clear();
    _adc15.clear();
    _adc16.clear();
    _adc17.clear();
    _adc18.clear();
    _adc19.clear();
    _adc20.clear();
    _adc21.clear();
    _adc22.clear();
    _adc23.clear();
    _adc24.clear();
    _adc25.clear();
    _adc26.clear();
    _adc27.clear();
    _adc28.clear();
    _adc29.clear();
    _adc30.clear();
    _adc31.clear();
    _wf00.clear();
    _wf01.clear();
    _wf02.clear();
    _wf03.clear();
    _wf04.clear();
    _wf05.clear();
    _wf06.clear();
    _wf07.clear();
    _wf08.clear();
    _wf09.clear();
    _wf10.clear();
    _wf11.clear();
    _wf12.clear();
    _wf13.clear();
    _wf14.clear();
    _wf15.clear();
    _wf16.clear();
    _wf17.clear();
    _wf18.clear();
    _wf19.clear();
    _wf20.clear();
    _wf21.clear();
    _wf22.clear();
    _wf23.clear();
    _wf24.clear();
    _wf25.clear();
    _wf26.clear();
    _wf27.clear();
    _wf28.clear();
    _wf29.clear();
    _wf30.clear();
    _wf31.clear();

    return;
  }
    
}
#endif
