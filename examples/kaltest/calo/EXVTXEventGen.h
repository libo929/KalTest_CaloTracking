#ifndef __EXVTXEVENTGEN__
#define __EXVTXEVENTGEN__

#include "EXVTXEventGen.h"
#include "EXVTXVKalDetector.h"
#include "EXVTXVMeasLayer.h"
//#include "IO/LCReader.h"
//#include  "grid_mnt/opt__exp_soft/ilc/ilcsoft/sl6/v01-17-09/lcio/v02-07/IOIMPL/LCFactory.h"
//#include "IMPL/LCCollectionVec.h"
//#include "EVENT/MCParticle.h"
//#include "EVENT/LCEvent.h"
//#include "EVENT/Track.h"
//#include "EVENT/CalorimeterHit.h"
//#include "EVENT/SimCalorimeterHit.h"
//#include "EVENT/SimTrackerHit.h"
//#include "UTIL/LCTOOLS.h"
//#include "Exceptions.h" 


#include "TPlane.h"
#include "TRandom.h"


#include <vector>
#include <string.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TClassTable.h>
#include <iostream>


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


   void LoadHits();

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
