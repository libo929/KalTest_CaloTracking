#include "CaloTrkHitClusteringProcessor.h"

#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <IMPL/TrackerHitImpl.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
#include <UTIL/Operators.h>

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>

//#include "TKalDetCradle.h"    // from KalTrackLib
#include "TKalTrackState.h"   // from KalTrackLib
#include "TKalTrackSite.h"    // from KalTrackLib
#include "TKalTrack.h"        // from KalTrackLib

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

static const Bool_t gkDir = kIterBackward;

CaloTrkHitClusteringProcessor aCaloTrkHitClusteringProcessor ;


CaloTrkHitClusteringProcessor::CaloTrkHitClusteringProcessor() : Processor("CaloTrkHitClusteringProcessor") {
  
  // modify processor description
  _description = "" ;
  
}


void CaloTrkHitClusteringProcessor::init() { 
  
  printParameters() ;
  _nRun = 0 ;
  _nEvt = 0 ;
}

void CaloTrkHitClusteringProcessor::processRunHeader( LCRunHeader* run) { 
  ++_nRun ;
} 

void CaloTrkHitClusteringProcessor::processEvent( LCEvent * evt ) { 

   std::cout << "evt: " << evt->getEventNumber() << std::endl; 
   
   clustering(evt);
	  
  _nEvt ++ ;
}

void CaloTrkHitClusteringProcessor::clustering(LCEvent* evt)
{
	cout << "--------------> In preClustering" << endl;

	LCCollection * ecalCol = evt->getCollection("ECALBarrel");// Ecal hits collection
	CellIDDecoder<CalorimeterHit> decoder(ecalCol);
    
	// output
	LCCollectionVec* _trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT );


	//std::vector< std::vector<CalorimeterHit* ecalHit> > allHits; 
	std::vector< CalorimeterHit* > hitArray[8][29]; // 8: staves, 29: layers

	int iHit = 0;

	for(int i=0; i<ecalCol->getNumberOfElements(); ++i) {
		++iHit;

		CalorimeterHit* ecalHit = dynamic_cast<CalorimeterHit*>(ecalCol->getElementAt(i));
		const float* hitPos = ecalHit->getPosition();

		int layerNum = 0, staveNum = 0;

		//FIXME::BitField64 ???
		layerNum = decoder(ecalHit)["K-1"];
		staveNum = decoder(ecalHit)["S-1"];

		std::cout << "CaloHitPos - " << iHit << ": " << hitPos[0] << ", " << hitPos[1] << ", " << hitPos[2] 
			      << ", layer: " << layerNum << ", stave: " << staveNum << std::endl;


		hitArray[staveNum][layerNum].push_back(ecalHit);
		//trkHit->setCovMatrix(covMat);      
	}

	//FIXME
	//we need gear file
	for(int s=0; s<8; ++s) {
		for(int l=0; l<29; ++l) {
			int hitSize = hitArray[s][l].size();

			if(hitSize>0) {
				cout << "stave: " << s << ", layer: " << l << ", hit size: " << hitSize << endl;

				std::vector< CalorimeterHit* >& hits = hitArray[s][l];

				while(!hits.empty()) { // for each hit in hits vector
					// a CaloTrackHit candidate (actually a CaloHit vector)
					std::vector< std::vector< CalorimeterHit* >::iterator > caloTrkHit;

					// initialization by the first hit in the hits vector
					caloTrkHit.push_back(hits.begin());

					std::vector< CalorimeterHit* >::iterator it = hits.begin() + 1;

					// the iterator of caloTrkHit :)
					std::vector< std::vector< CalorimeterHit* >::iterator >::iterator caloTrkHitIt;

					// for each other hit, we are going to decide if it can be put into the
					// caloTrkHit by distance
					for(; it!=hits.end(); ++it) { 
				    	const float* hitPos1 = (*it)->getPosition(); 
						
						caloTrkHitIt = caloTrkHit.begin();

						bool isHit = true;

						// compare the distance to hits which is in caloTrHit already
						for(; caloTrkHitIt!=caloTrkHit.end();++caloTrkHitIt) { // for each hit in the candidate
							const float* hitPos2 = (*(*caloTrkHitIt))->getPosition(); 

				    	    double distance = sqrt( (hitPos1[0]-hitPos2[0])*(hitPos1[0]-hitPos2[0]) +
				    				                (hitPos1[1]-hitPos2[1])*(hitPos1[1]-hitPos2[1]) + 
				    								(hitPos1[2]-hitPos2[2])*(hitPos1[2]-hitPos2[2]) );
		
							std::cout << "hit1: " << hitPos1[0] << ", " << hitPos1[1] << ", " << hitPos1[2] << std::endl;
							std::cout << "hit2: " << hitPos2[0] << ", " << hitPos2[1] << ", " << hitPos2[2] << std::endl;
				    		
				    		// FIXME:: hard coded ...
							// The distance is discrete, so the value of this cut is very safe
				    		if(distance>7.15) { //5.1 * sqrt(2)
								std::cout << "distance: " <<  distance << std::endl;
								isHit = false; break;
							}	
						} 
						
						// this CaloHit is a hit belonging to CaloTrkHit
						if(isHit) caloTrkHit.push_back(it);
					}

					// get CaloTrkHit information
					double caloHitPos[3];
					double hitEnergy = 0.;

					for(caloTrkHitIt=caloTrkHit.begin(); caloTrkHitIt!=caloTrkHit.end(); ++caloTrkHitIt) {
						const float* hp = (**caloTrkHitIt)->getPosition();;
		                std::cout << "Hit components - " << hp[0] << ", " << hp[1] << ", " << hp[2]  << std::endl;

						caloHitPos[0] = hp[0];
						caloHitPos[1] = hp[1];
					    caloHitPos[2] = hp[2];
					}

					// store information of CaloTrkHit into LCIO
                    TrackerHitImpl* trkHit = new TrackerHitImpl ;
                    trkHit->setPosition(caloHitPos);
                    trkHit->setEDep(10.);
					trkHit->setCellID0(l); //layer
					_trkhitVec->addElement( trkHit ); 

					if(caloTrkHit.size()>1) std::cout << "CaloTrackHit has more than two components " << endl;

					// Very tricky: it is only work if we erase from end to begin
					std::vector< std::vector< CalorimeterHit* >::iterator >::reverse_iterator caloTrkHitRit;

					for(caloTrkHitRit=caloTrkHit.rbegin(); caloTrkHitRit!=caloTrkHit.rend(); ++caloTrkHitRit) {
						hits.erase( *caloTrkHitRit );
					}

					std::cout << ">>>>>" << std::endl;
				}
			}
		}// layer
	}// stave

	evt->addCollection( _trkhitVec , "CaloTrackHit" ) ;

	cout << "--------------> End of preClustering" << endl;
}

void CaloTrkHitClusteringProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void CaloTrkHitClusteringProcessor::end(){ 

  streamlog_out(MESSAGE) << " end()  " << name() 
  << " processed " << _nEvt << " events in " << _nRun << " runs "
  << std::endl ;
  
}
