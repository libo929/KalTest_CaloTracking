#ifndef __EXVTXEVENTGEN__
#define __EXVTXEVENTGEN__

#include "TKalDetCradle.h"
#include "THelicalTrack.h"

#include "TVector3.h"
#include <list>

class EXVTXEventGen {
public:
   EXVTXEventGen(TKalDetCradle &cradle, TObjArray &kalhits)
             : fCradlePtr(&cradle), fHitBufPtr(&kalhits) {}
   virtual ~EXVTXEventGen() {}

   THelicalTrack GenerateHelix(Double_t pt,
                               Double_t cosmin,
                               Double_t cosmax);
   void          Swim(THelicalTrack &heltrk);

   static void     SetT0(Double_t t0) { fgT0 = t0;   }
   static Double_t GetT0()            { return fgT0; }
   std::vector<TVector3> GetHitVec()  { return fHitVec; }

private:
   TKalDetCradle *fCradlePtr;     // pointer to detector system
   TObjArray     *fHitBufPtr;     // pointer to hit array
   std::vector<TVector3> fHitVec;      // to store 3-D hits

   static Double_t  fgT0;         // t0

   ClassDef(EXVTXEventGen,1)   // Event Generator
};

#endif
