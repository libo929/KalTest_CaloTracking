#ifndef __EXVTXVMEASLAYER__
#define __EXVTXVMEASLAYER__
//*************************************************************************
//* ===================
//*  EXVTXVMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by TVTrackHit.
//* (Requires)
//* (Provides)
//*     class EXVTXVMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "TCylinder.h"
#include "TVMeasLayer.h"
#include "TAttDrawable.h"
#include "KalTrackDim.h"
#include "TString.h"

class TVTrackHit;
class TNode;

class EXVTXVMeasLayer : public TVMeasLayer, public TAttDrawable {
public:
   static Bool_t kActive;
   static Bool_t kDummy;

   // Ctors and Dtor

   EXVTXVMeasLayer(TMaterial &min,
                TMaterial &mout,
                Bool_t     type = EXVTXVMeasLayer::kActive,
          const Char_t    *name = "MeasL");
   virtual ~EXVTXVMeasLayer();

   virtual void ProcessHit(const TVector3  &xx,
                                 TObjArray &hits) = 0;

   inline TString GetMLName () const { return fName;    }
   inline TNode  *GetNodePtr() const { return fNodePtr; }

   inline void    SetNodePtr(TNode *nodep) { fNodePtr = nodep; }

private:
   TString  fName;      // layer name
   TNode   *fNodePtr;   // node pointer

   ClassDef(EXVTXVMeasLayer,1)     // Sample measurement layer class
};

#endif
