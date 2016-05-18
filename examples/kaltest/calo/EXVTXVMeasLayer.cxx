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

#include "EXVTXVMeasLayer.h"

Bool_t   EXVTXVMeasLayer::kActive = kTRUE;
Bool_t   EXVTXVMeasLayer::kDummy = kFALSE;

ClassImp(EXVTXVMeasLayer)
                                                                                
EXVTXVMeasLayer::EXVTXVMeasLayer(TMaterial &min,
                           TMaterial &mout,
                           Bool_t     isactive,
                     const Char_t    *name)  
            : TVMeasLayer(min, mout, isactive),
	      fName(name),
	      fNodePtr(0)
{
}

EXVTXVMeasLayer::~EXVTXVMeasLayer()
{
}
