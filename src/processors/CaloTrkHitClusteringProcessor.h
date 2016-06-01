#ifndef CaloTrkHitClusteringProcessor_h
#define CaloTrkHitClusteringProcessor_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;


/** ======= CaloTrkHitClusteringProcessor ========== <br>
 * 
 */
class CaloTrkHitClusteringProcessor : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new CaloTrkHitClusteringProcessor ; }
  
  
  CaloTrkHitClusteringProcessor() ;
  
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
  void clustering(LCEvent* evt);
  
  
protected:
  
  std::string _inColName ;
  
  std::string _outColName ;
  std::string _outRelColName ;
 
    
  int _nRun ;
  int _nEvt ;
} ;

#endif
