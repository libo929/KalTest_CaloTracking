#ifndef CaloTrkBuildingProcessor_h
#define CaloTrkBuildingProcessor_h 1

#include "marlin/Processor.h"

#include "EXVTXEventGen.h"
#include "EXVTXKalDetector.h"

#include "TNtupleD.h"         // from ROOT
#include "TFile.h"            // from ROOT

#include "TKalDetCradle.h"    // from KalTrackLib
#include "lcio.h"

#include <string>
#include <vector>

#include <IMPL/LCCollectionVec.h>

using namespace lcio ;
using namespace marlin ;

class CaloTrkBuildingProcessor : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new CaloTrkBuildingProcessor ; }
  
  
  CaloTrkBuildingProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;


private:
  void tracking(LCEvent* evt);
  void kalhitMaking( LCCollection* tv,  TObjArray& hitbuf);
  
  
protected:
  
  std::string _inColName ;
  
  std::string _outColName ;
  std::string _outRelColName ;

  TObjArray     _kalhits;    // hit buffer 
  TKalDetCradle _cradle;     // detctor system
  EXVTXKalDetector _detector;   
 
  int _nRun ;
  int _nEvt ;
} ;

#endif
