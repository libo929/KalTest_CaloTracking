#ifndef __EXVTXDETECTOR__
#define __EXVTXDETECTOR__

#include "EXVTXVKalDetector.h"
#include "TMath.h"

class EXVTXKalDetector : public EXVTXVKalDetector {
  
 public:
  EXVTXKalDetector(Int_t m = 100);
  ~EXVTXKalDetector();
  
  ClassDef(EXVTXKalDetector,1)   
    
    
  static const Int_t _nLayer = 6;

};

#endif
