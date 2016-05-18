#include "EXVTXVKalDetector.h"
#include "EXVTXVMeasLayer.h"
#include "TVKalDetector.h"
#include "TCONE.h" // this is to include TTUBE header
//#include "TTUBE.h"
#include "TNode.h"
#include "TRotMatrix.h"
#include "TVirtualPad.h"


Double_t EXVTXVKalDetector::fgBfield  = 3.5;
TNode   *EXVTXVKalDetector::fgNodePtr = 0;

ClassImp(EXVTXVKalDetector)

EXVTXVKalDetector::EXVTXVKalDetector(Int_t m)
             : TVKalDetector(m),
               fIsPowerOn(kTRUE)
{
}

EXVTXVKalDetector::~EXVTXVKalDetector()
{
}

TNode *EXVTXVKalDetector::GetNodePtr()
{
   if (!fgNodePtr) {
      new TRotMatrix("rotm","rotm", 10.,80.,10.,80.,10.,80.);
      //new TTUBE("Det","Det","void",100.,100.,260.);
      new TTUBE("Det","Det","void",2500.,2500.,2500.);
      fgNodePtr = new TNode("World","World","Det",0.,0.,0.,"rotm");
   }
   return fgNodePtr;
}

void EXVTXVKalDetector::Draw(Int_t color, const Char_t *opt)
{
   if (!gPad) return;
   TNode *nodep = GetNodePtr();
   nodep->cd();
   TIter next(this);
   TObject *objp;
   while ((objp = next())) {
      TAttDrawable *dp = dynamic_cast<TAttDrawable *>(objp);
      if (dp) dp->Draw(color, opt); 
   }
   nodep->Draw("pad same");
}
