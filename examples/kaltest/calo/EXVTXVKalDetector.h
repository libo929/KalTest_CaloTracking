#ifndef __EXVTXVDETECTOR__
#define __EXVTXVDETECTOR__

#include "TVector3.h"
#include "TVKalDetector.h"
#include "TAttDrawable.h"
#include "TNode.h"

class TVMeasLayer;

class EXVTXVKalDetector : public TVKalDetector, public TAttDrawable {
public:
   EXVTXVKalDetector(Int_t m = 100);
   virtual ~EXVTXVKalDetector();

   inline virtual Bool_t IsPowerOn() const { return fIsPowerOn;   }
   inline virtual void   PowerOn  ()       { fIsPowerOn = kTRUE;  }
   inline virtual void   PowerOff ()       { fIsPowerOn = kFALSE; }

   static Double_t GetBfield (const TVector3 &xx = TVector3(0., 0., 0.))
                             { return fgBfield; }

   using  TAttDrawable::Draw;
   virtual void Draw(Int_t color, const Char_t *opt = "");

   static void   SetNodePtr(TNode *nodep) { fgNodePtr = nodep; }
   static TNode *GetNodePtr();

private:
   Bool_t  fIsPowerOn;         // power status
   static Double_t fgBfield;   // magnetic field [kG]
   static TNode   *fgNodePtr;  // pointer to TNode

   ClassDef(EXVTXVKalDetector,1)   // Sample hit class
};

#endif
