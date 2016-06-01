#include "EXVTXKalDetector.h"
#include "EXVTXMeasLayer.h"
#include "EXVTXHit.h"
#include "TRandom.h"
#include "TMath.h"
#include <sstream>


#include "marlin/Global.h"
#include "gearxml/GearXML.h"
#include "gear/LayerLayout.h"
#include "gear/CalorimeterParameters.h"
#include "gear/GearMgr.h"

using namespace gear;

ClassImp(EXVTXKalDetector)

EXVTXKalDetector::EXVTXKalDetector(Int_t m)
                : EXVTXVKalDetector(m)
{
   Double_t A, Z, density, radlen;
   A       = 14.00674 * 0.7 + 15.9994 * 0.3;
   Z       = 7.3;
   density = 1.205e-3;
   radlen  = 3.42e4;
   TMaterial &air = *new TMaterial("VTXAir", "", A, Z, density, radlen, 0.);

   A       = 28.0855;
   Z       = 14.;
   density = 2.33;   
   radlen  = 9.36;
   TMaterial &si = *new TMaterial("VTXSi", "", A, Z, density, radlen, 0.);
   
   // Dev: Adding a material: Tungsten
   A	   = 183.84;
   Z	   = 74;
   density = 19.3;
   radlen =  0.3504; 
   TMaterial &w = *new TMaterial("VTXW", "", A, Z, density, radlen, 0.);
   // endDev

   Bool_t active = EXVTXMeasLayer::kActive;
   Bool_t dummy  = EXVTXMeasLayer::kDummy;
   Bool_t activeness[2] = { active, dummy};

   GearXML gearXML("ILD_o1_v05.xml");
   GearMgr* gearMgr = gearXML.createGearMgr();

   std::cout << "gearMgr: " << gearMgr << std::endl;

   const CalorimeterParameters& ecalBarrel = gearMgr->getEcalBarrelParameters();
   const LayerLayout& ecalBarrelLayerLayout = ecalBarrel.getLayerLayout();
   const std::vector<double_t>& ecalBarrelExtent = ecalBarrel.getExtent();

   int nLayer = ecalBarrelLayerLayout.getNLayers();
   
   Double_t inner_R = ecalBarrelExtent[0];
   Double_t outer_z = ecalBarrelExtent[3];

//====================================   
//     geometry of 3 doublet vtx 
//====================================

   //InitGeometry();

   //
   //static const Double_t sigmaxi   = 1.44e-4;
   //static const Double_t sigmazeta = 1.44e-4;
   
   //
   //static const Double_t sigmaxi   = 0.15;
   //static const Double_t sigmazeta = 0.15;
   static const Double_t sigmaxi   = 1.5;
   static const Double_t sigmazeta = 1.5;
	
   Double_t yDepth = inner_R;
   Double_t xiWidth = 2 * yDepth * tan(TMath::Pi()/8);			// 8 pieces of calorimeters, each share 40 degrees
   Double_t zetaWidth = outer_z;

   for(int l=0; l<nLayer; ++l) {
	   
	   // 8 different orientations
	   for (int j=0; j<8; j++){									//to be read from GEAR
			Double_t angle=TMath::Pi()*j/4;
			TVector3 normal(sin(angle), cos(angle), 0);
			TVector3 xc(yDepth*sin(angle), yDepth*cos(angle), 0);
			
			Add(new EXVTXMeasLayer(si, w, xc, normal, xc.Perp(), xiWidth, zetaWidth, 0., sigmaxi, sigmazeta));
	   }

	   
	   /** 
	   TVector3 normal(0, 1, 0);
       TVector3 xc(0, -yDepth, 0);
	   TVector3 xc1(0, yDepth, 0);

	   Add(new EXVTXMeasLayer(si, w, xc, normal, xc.Perp(), xiWidth, zetaWidth, 0.,
	   //Add(new EXVTXMeasLayer(air, si, xc, normal, xc.Perp(), 200+l*30, 200+l*30, 0.,
				              sigmaxi, sigmazeta));       

	   Add(new EXVTXMeasLayer(si, w, xc1, normal, xc.Perp(), xiWidth, zetaWidth, 0., sigmaxi, sigmazeta));
	   **/

	   yDepth += ecalBarrelLayerLayout.getThickness(l);
	   xiWidth = 2 * yDepth * tan(TMath::Pi()/8);
	}
   SetOwner();			
}

EXVTXKalDetector::~EXVTXKalDetector()
{
}
