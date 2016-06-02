#include "CaloTrkBuildingProcessor.h"

#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <IMPL/TrackImpl.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
#include <UTIL/Operators.h>

#include <marlin/Global.h>
#include <gear/GEAR.h>

#include "TKalTrackState.h"   // from KalTrackLib
#include "TKalTrackSite.h"    // from KalTrackLib
#include "TKalTrack.h"        // from KalTrackLib

#include "EXVTXHit.h"


#include <iostream>


#include "gear/gearsurf/MeasurementSurfaceStore.h"
#include "gear/gearsurf/MeasurementSurface.h"
#include "gear/gearsurf/ICoordinateSystem.h"
#include "gear/gearsurf/CartesianCoordinateSystem.h"

//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "marlin/ProcessorEventSeeder.h"

#include "UTIL/ILDConf.h"
#include "UTIL/BitSet32.h"

#include "CLHEP/Vector/TwoVector.h"

#include <cmath>
#include <algorithm>
#include <sstream>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

//static const Bool_t gkDir = kIterBackward;
static const Bool_t gkDir = kIterForward;

CaloTrkBuildingProcessor aCaloTrkBuildingProcessor ;


CaloTrkBuildingProcessor::CaloTrkBuildingProcessor() : Processor("CaloTrkBuildingProcessor") {
  
  // modify processor description
  _description = "CaloTrkBuildingProcessor creates TrackerHits from SimTrackerHits, smearing them according to the input parameters." ;
  
  
  // register steering parameters: name, description, class-variable, default value
}


void CaloTrkBuildingProcessor::init() { 
  
  printParameters() ;
  _nRun = 0 ;
  _nEvt = 0 ;

  _cradle.Install(_detector); // install detector into its cradle
  //_cradle.SwitchOffMS();     // switch off multiple scattering
  _cradle.SwitchOnMS();     // switch off multiple scattering
  _cradle.SwitchOffDEDX();     // switch off multiple scattering
}

void CaloTrkBuildingProcessor::processRunHeader( LCRunHeader* run) { 
  ++_nRun ;
} 

void CaloTrkBuildingProcessor::processEvent( LCEvent * evt ) { 

   std::cout << "evt: " << evt->getEventNumber() << std::endl; 
   
   tracking(evt);
	  
  _nEvt ++ ;
}

void CaloTrkBuildingProcessor::kalhitMaking( LCCollection* tv,  TObjArray& hitbuf)
{
	// buffer initialization
	hitbuf.Delete();

	//......
	int iHit = 0;

	cout << "----------->>>>>>>> hit number: " << tv->getNumberOfElements() << endl;

	for(int i=0; i<tv->getNumberOfElements(); ++i) {
		++iHit;

		TrackerHit* caloTrkHit = dynamic_cast<TrackerHit*>(tv->getElementAt(i));

		const double* hitPos = caloTrkHit->getPosition();
		int layer = caloTrkHit->getCellID0(); // layer: by convention
		cout << ">>>>>>>>>>>>> hitPos: " << hitPos[0] << ", " << hitPos[1] << ", " << hitPos[2] 
			 << ", layer: " << layer << std::endl;

		TVector3 xx(hitPos[0], hitPos[1], hitPos[2]);

		// get layer
		EXVTXMeasLayer* measLayer = dynamic_cast<EXVTXMeasLayer*>(_detector[layer]);
		//cout << "measlayer: " <<  measLayer->GetXiwidth() << endl;
		cout << "measlayer center: " <<  endl;
		measLayer->GetXc().Print();
        TKalMatrix h    = measLayer->XvToMv(xx);

        Double_t xi = h(0,0);
        Double_t zeta = h(1,0);
	
        Double_t meas [2];
        Double_t dmeas[2];
        meas [0] = xi;
        meas [1] = zeta;
        dmeas[0] = 15;
        dmeas[1] = 15;

        Double_t b = EXVTXKalDetector::GetBfield();
        hitbuf.Add(new EXVTXHit(*measLayer, meas, dmeas, xx, b));
	}

	hitbuf.SetOwner(kTRUE);
}

void CaloTrkBuildingProcessor::tracking(LCEvent* evt)
{
	cout << "--------------> In CaloTrkBuilding::tracking" << endl;

	// input
	LCCollection* caloHitCol = evt->getCollection("CaloTrackHit");
	LCCollection* cluTrack   = evt->getCollection("ClupatraTracks");
    
	// output
	LCCollectionVec* _trkVec = new LCCollectionVec( LCIO::TRACK );

	LCFlagImpl hitFlag(0) ;
	hitFlag.setBit( LCIO::TRBIT_HITS ) ;
	_trkVec->setFlag( hitFlag.getFlag()  ) ;

	//TODO:
	// we should make a loop for _kalhits
	
	kalhitMaking(caloHitCol, _kalhits);

	bool buildTrk = false;

#if 1
    // ---------------------------
    //  Create a dummy site: sited
    // ---------------------------

    EXVTXHit hitd = *dynamic_cast<EXVTXHit *>(_kalhits.At(0));
	cout << "dummy hit:" << endl;
	hitd.GetExactX().Print();
    hitd(0,1) = 1.e6;   // give a huge error to d
    hitd(1,1) = 1.e6;   // give a huge error to z

    TKalTrackSite &sited = *new TKalTrackSite(hitd);
    sited.SetOwner();   // site owns states

    // ---------------------------
    // Create initial helix
    // ---------------------------
	// maybe we can also start from MarlinTrack ...
    Int_t i1, i2, i3; // (i1,i2,i3) = (1st,mid,last) hit to filter

    if (gkDir == kIterBackward) {
       i3 = 0;
       i1 = _kalhits.GetEntries() - 1;
       i2 = i1 / 2;
    } else {
       i1 = 0;
       i3 = _kalhits.GetEntries() - 1;
       i2 = i3 / 2;
    }

    TVTrackHit &h1 = *dynamic_cast<TVTrackHit *>(_kalhits.At(i1)); // first hit
    TVTrackHit &h2 = *dynamic_cast<TVTrackHit *>(_kalhits.At(i2)); // middle hit
    TVTrackHit &h3 = *dynamic_cast<TVTrackHit *>(_kalhits.At(i3)); // last hit
    TVector3    x1 = h1.GetMeasLayer().HitToXv(h1);
    TVector3    x2 = h2.GetMeasLayer().HitToXv(h2);
    TVector3    x3 = h3.GetMeasLayer().HitToXv(h3);
    THelicalTrack helstart(x1, x2, x3, h1.GetBfield(), gkDir); // initial helix 

    // ---------------------------
    //  Set dummy state to sited
    // ---------------------------

    static TKalMatrix svd(kSdim,1);
    svd(0,0) = 0.;
    svd(1,0) = helstart.GetPhi0();
	//FIXME::for test
    svd(2,0) = helstart.GetKappa();
    svd(3,0) = 0.;
    svd(4,0) = helstart.GetTanLambda();
    if (kSdim == 6) svd(5,0) = 0.;

#if 1
	// track information from tracker
	if(cluTrack->getNumberOfElements()==1) {
		Track* trk = dynamic_cast<Track*>(cluTrack->getElementAt(0));
		double omega = trk->getOmega();
		double cpa   = omega * 1000 *1.e9/(3.5 * TMath::C());
		double tanl  = trk->getTanLambda();
		double phi0  = trk->getPhi() - M_PI/2; 
		
	    while (phi0 < -M_PI) phi0 += 2.0*M_PI;
	    while (phi0 >= M_PI) phi0 -= 2.0*M_PI;

		//cout << "cpa: " << cpa << ", tanl: " << tanl << endl;

		if(cpa>-1&cpa<1) {
			svd(1,0) = phi0;
			svd(2,0) = cpa;
			svd(4,0) = tanl;
		}
	}
#endif

	cout << "cpa.init: " << svd(2,0) << endl;

    static TKalMatrix C(kSdim,kSdim);
    for (Int_t i=0; i<kSdim; i++) {
       //C(i,i) = 1.e4;   // dummy error matrix
       C(i,i) = 1;   // dummy error matrix
    }

#if 1
	// reset the covariance matrix
	C(0,0) = 1;
	C(1,1) = 1;
	C(2,2) = 0.01;//0.001
	C(3,3) = 1;
	C(4,4) = 0.01;//0.01
#endif

    sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
    sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

    // ---------------------------
    //  Add sited to the kaltrack
    // ---------------------------

    TKalTrack kaltrack;    // a track is a kal system
    kaltrack.SetOwner();   // kaltrack owns sites
    kaltrack.Add(&sited);  // add the dummy site to the track

    // ---------------------------
    //  Prepare hit iterrator
    // ---------------------------

    TIter next(&_kalhits, gkDir);   // come in to IP

	cout << "kalhits size: " <<  _kalhits.GetEntries() << endl;

    // ---------------------------
    //  Start Kalman Filter
    // ---------------------------

    EXVTXHit *hitp = 0;
    while ((hitp = dynamic_cast<EXVTXHit *>(next()))) {     // loop over hits
       TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
	   TKalTrackSite& lastSite = dynamic_cast<TKalTrackSite&>(kaltrack.GetCurSite());


	   //cout << "last track hit: " << endl;
	   //lastHit.Print();
	   //cout << "this hit: " << endl;
	   //hitp->GetExactX().Print();
	   
	   TVector3 lastHit = lastSite.GetPivot();
	   TVector3 thisHit = hitp->GetExactX();
	   TVector3 diff    = lastHit - thisHit;

	   // FIXME
	   // very naive!!!
	   double hitDist = sqrt(diff.X()*diff.X() + diff.Z()*diff.Z());
	   cout << "hits dist: " << hitDist <<  endl;

       //only for testing code...


	   //FIXME:: 20 ???
	   if(hitDist>20) { 
		   thisHit.Print();
		   cout << "hit distance : " << hitDist << ", hit rejected!" << endl;
		   cout << endl;
		   continue;
	   }

       if (!kaltrack.AddAndFilter(site)) {               // add and filter this site
          cout << " site discarded!" << endl;           
		  hitp->GetExactX().Print();
          delete &site;                                  // delete this site, if failed
       } 
	   else {
		   cout << "a hit is added to track with DeltaCHI2: " << site.GetDeltaChi2() << endl;
		  hitp->GetExactX().Print();
	      Double_t cpa = kaltrack.GetCurSite().GetCurState()(2, 0);
		  cout << "cpa: " << cpa << endl;
	   }
       cout << endl;           
    }

	const int MINTRACKHITNUM = 25;

	if(kaltrack.GetEntries() > MINTRACKHITNUM) buildTrk = true;

    // ============================================================
    //  Monitor Fit Result
    // ============================================================

    Int_t    ndf  = kaltrack.GetNDF();
    Double_t chi2 = kaltrack.GetChi2();
    Double_t cl   = TMath::Prob(chi2, ndf);

	// get the parameters at last point
	TVKalSite& cursite = kaltrack.GetCurSite();

    Double_t d0  =  cursite.GetCurState()(0, 0);  
    Double_t fi0  = cursite.GetCurState()(1, 0);  
    Double_t cpa  = cursite.GetCurState()(2, 0);  
    Double_t dz   = cursite.GetCurState()(3, 0);  
    Double_t tnl  = cursite.GetCurState()(4, 0);  
   
    Double_t omega = cpa/1000/1.e9*(3.5 * TMath::C());
	TVector3 ref = dynamic_cast<TKalTrackSite&>(cursite).GetPivot();

#if 1
	// try to get parameters at the starting point
	TVTrack& ft = dynamic_cast<TKalTrackState&>(cursite.GetCurState()).CreateTrack();
	TKalTrackSite& firstSite = dynamic_cast<TKalTrackSite&>(*kaltrack.At(0));
	TVector3 refs = dynamic_cast<TKalTrackSite&>(firstSite).GetPivot();
	
	// with (0,0,0)
	//refs = TVector3();

	Double_t angle = 0;
	ft.MoveTo(refs, angle);

	Double_t d0s  = ft.GetDrho();
	Double_t fi0s = ft.GetPhi0();
	Double_t cpas = ft.GetKappa();
	Double_t dzs  = ft.GetDz();
	Double_t tnls = ft.GetTanLambda();
	cout << "ft MoveTo done with a angle: " <<  angle << ", d0s: " << d0s << ", dzs: " << dzs 
		 << ", fi0s: " << fi0s << endl;

    Double_t omegas = cpas/1000/1.e9*(3.5 * TMath::C());

	delete &ft;

	float refPoss[3] = {(float)refs.X(), (float)refs.Y(), (float)refs.Z()};

	if( buildTrk ) {
        TrackImpl* trks = new TrackImpl;
		trks->setReferencePoint(refPoss);
		trks->setD0(d0s);
		trks->setPhi(fi0s + M_PI/2);
		trks->setOmega(omegas);
		trks->setZ0(dzs);
		trks->setTanLambda(tnls);
		cout << "-------omegas: " << omegas << ", d0s: " << d0s << endl;

		// FIXME::
		// only add the hit belong to a track
		for(int i=0; i<caloHitCol->getNumberOfElements(); ++i) {
			TrackerHit* caloTrkHit = dynamic_cast<TrackerHit*>(caloHitCol->getElementAt(i));
			trks->addHit(caloTrkHit);
		}

    	_trkVec->addElement( trks ); 
	}
#endif

	float refPos[3] = {(float)ref.X(), (float)ref.Y(), (float)ref.Z()};

	cout << "ndf: " << ndf << ", chi2: " << chi2 << ", cpa: " << cpa << endl;
#endif

	// Smooth, seems doesn't work ???
	//kaltrack.SmoothBackTo(1);
#if 0
	if( buildTrk ) {
    	// store information of CaloTrkHit into LCIO
        TrackImpl* trk = new TrackImpl;
		trk->setReferencePoint(refPos);
		trk->setD0(d0);
		trk->setPhi(fi0 + M_PI/2);
		trk->setOmega(omega);
		trk->setZ0(dz);
		trk->setTanLambda(tnl);

		cout << "-------omega: " << omega << endl;

		for(int i=0; i<caloHitCol->getNumberOfElements(); ++i) {
			TrackerHit* caloTrkHit = dynamic_cast<TrackerHit*>(caloHitCol->getElementAt(i));
			trk->addHit(caloTrkHit);
		}

        //trk setting...
    
    	_trkVec->addElement( trk ); 
	}
#endif

	evt->addCollection( _trkVec , "CaloTrack" ) ;

	cout << "--------------> End of tracking, " << buildTrk << endl;
}

void CaloTrkBuildingProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void CaloTrkBuildingProcessor::end(){ 

  streamlog_out(MESSAGE) << " end()  " << name() 
  << " processed " << _nEvt << " events in " << _nRun << " runs "
  << std::endl ;
  
}
