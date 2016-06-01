#include "EXVTXEventGen.h"
#include "EXVTXVKalDetector.h"
#include "EXVTXVMeasLayer.h"
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "IMPL/LCCollectionVec.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCEvent.h"
#include "EVENT/Track.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/SimTrackerHit.h"
#include "UTIL/LCTOOLS.h"
#include "Exceptions.h" 


#include "TPlane.h"
#include "TRandom.h"


#include <vector>
#include <string.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TClassTable.h>
#include <iostream>


//-----------------------------------
// Track Parameters
//-----------------------------------

#define __DR__     0.
#define __FI0__    0.
#define __DZ__     0.
#define __X0__     0.
#define __Y0__     0.
#define __Z0__     0.

ClassImp(EXVTXEventGen)

Double_t EXVTXEventGen::fgT0 = 14.; // [nsec]

void EXVTXEventGen::LoadHits()
{
	const char* FILEN = "hitsimu.slcio";
	std::string rootFileBaseName( FILEN , strlen(FILEN)-strlen(".slcio") ) ;

	std::string dirname("");
	int dirlen = strlen( dirname.c_str() );
	int baselen = strlen( rootFileBaseName.c_str() );
	std::string tmpname(rootFileBaseName, dirlen, baselen-dirlen);
	
	//IN case this method is executed multiple times
	delete gROOT->GetListOfFiles()->FindObject( FILEN );
	delete gROOT->GetListOfCanvases()->FindObject("c1");

	int nEvents  = 0;
	int maxEvt   = 2000;

	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
	lcReader->open(FILEN);

	const double b = 3.5;			// Hardcoded M-field strength 3.5 T
	EVENT::LCEvent* evt = 0;
	while ( (evt=lcReader->readNextEvent())!=0 && nEvents < maxEvt ){
		nEvents++;

		
	}
	std::cout << "In total, there are : " << nEvents << " events." << std::endl;
}

THelicalTrack EXVTXEventGen::GenerateHelix(Double_t pt,
                                        Double_t cosmin,
                                        Double_t cosmax)
{
   // ---------------------------
   //  Generate a helical track
   // ---------------------------

   std::cout << "cosmin :" << cosmin << ", cosmax :" << cosmax << std::endl;
   Double_t dr  = __DR__;
   Double_t fi0 = __FI0__ + 2*TMath::Pi()*(gRandom->Uniform()-0.5);
   //fi0 = __FI0__;
   Double_t cpa = 1. / pt;
   Double_t dz  = __DZ__;
   Double_t cs  = gRandom->Uniform(cosmin, cosmax);
   Double_t tnl = cs / TMath::Sqrt((1-cs)*(1+cs)); 
   Double_t x0  = __X0__;
   Double_t y0  = __Y0__;
   Double_t z0  = __Z0__;

   Double_t b   = dynamic_cast<const EXVTXVKalDetector &>
                 (dynamic_cast<EXVTXVMeasLayer *>
                 (fCradlePtr->At(0))->GetParent(kFALSE)).GetBfield();

   return THelicalTrack(dr,fi0,cpa,dz,tnl,x0,y0,z0,b);
}

void EXVTXEventGen::Swim(THelicalTrack &heltrk)
{
   // ---------------------------
   //  Swim track and Make hits
   // ---------------------------

   Double_t dfi       = -dynamic_cast<TVSurface *>(fCradlePtr->At(0)) 
	                             ->GetSortingPolicy()
                         / heltrk.GetRho();

   Int_t    nlayers   = fCradlePtr->GetEntries();
   Int_t    dlyr      = 1;
   Double_t dfisum    = 0.;

   fHitVec.clear();

   for (Int_t lyr = 0; lyr >= 0; lyr += dlyr) { // loop over layers
      // change direction if it starts looping back
      if (lyr == nlayers - 1) dlyr = -1;

      EXVTXVMeasLayer &ml = *dynamic_cast<EXVTXVMeasLayer *>(fCradlePtr->At(lyr));
      TVSurface    &ms = *dynamic_cast<TVSurface *>(fCradlePtr->At(lyr));
      TVector3 xx;
      Double_t dfis = dfi;
      if (!ms.CalcXingPointWith(heltrk,xx,dfi,1)
       || TMath::Abs(dfi) > TMath::Pi()
       || TMath::Abs(dfi + dfisum) > TMath::TwoPi()) {
         dfi = dfis;
         continue;
      }
      // should use the material behind the surface since dfi is measured 
      // from the last point to the current surface
      Bool_t   dir    = dlyr < 0 ? kTRUE : kFALSE;

      if (fCradlePtr->IsMSOn()) {
         TKalMatrix Qms(5,5);
         ml.CalcQms(dir, heltrk, dfi, Qms);
         Double_t sgphi  = TMath::Sqrt(Qms(1,1));
         Double_t sgtnl  = TMath::Sqrt(Qms(4,4));
         Double_t delphi = gRandom->Gaus(0.,sgphi);
         Double_t deltnl = gRandom->Gaus(0.,sgtnl);
#if 0
         dfi *= 0.5;
         TVector3 x0ms = heltrk.CalcXAt(dfi);
         heltrk.MoveTo(x0ms,dfi);     // M.S. at mid point

         heltrk.ScatterBy(delphi,deltnl);
         dfis = dfi;
#else
         heltrk.ScatterBy(delphi,deltnl); // multiple scattering
         dfis = 0.;
#endif
         // recalculate crossing point
         if (!ms.CalcXingPointWith(heltrk,xx,dfi,1)
          || TMath::Abs(dfi) > TMath::Pi()
          || TMath::Abs(dfi + dfisum) > TMath::TwoPi()) {
            dfi = dfis;
            continue;
         }
         dfis += dfi;
      }
      dfisum += dfi;

      heltrk.MoveTo(xx,dfi);	// move pivot to current hit

      if (fCradlePtr->IsDEDXOn()) {
         TKalMatrix av(5,1);
         heltrk.PutInto(av);
         av(2,0) += ml.GetEnergyLoss(dir, heltrk, dfis); // energy loss
         heltrk.SetTo(av, heltrk.GetPivot());
      }
      if (ml.IsActive() && dynamic_cast<const EXVTXVKalDetector &>(ml.GetParent(kFALSE)).IsPowerOn()) {
         ml.ProcessHit(xx, *fHitBufPtr); // create hit point
		 fHitVec.push_back(xx);
		 xx.Print();
      }
      if (lyr == nlayers - 1) break;
   }
}
