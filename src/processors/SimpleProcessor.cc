#include "SimpleProcessor.h"

//#include <marlin/Global.h>

using namespace marlin ;
using namespace std ;

SimpleProcessor aSimpleProcessor ;


SimpleProcessor::SimpleProcessor() : Processor("SimpleProcessor") {
}


void SimpleProcessor::init() { 
  
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

}

void SimpleProcessor::processRunHeader( LCRunHeader* run) { 

  ++_nRun ;

} 

void SimpleProcessor::processEvent( LCEvent * evt ) { 

  streamlog_out(MESSAGE) << " >>>>>> Event:  " << _nEvt << std::endl;
  _nEvt ++ ;

}



void SimpleProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void SimpleProcessor::end(){ 
  
  streamlog_out(MESSAGE) << " end()  " << name() 
  << " processed " << _nEvt << " events in " << _nRun << " runs "
  << std::endl ;
  
}
