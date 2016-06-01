//#define SAVE_RESIDUAL
#include "TNtupleD.h"
#include "TFile.h"

#include "TKalDetCradle.h"
#include "TKalTrackState.h"
#include "TKalTrackSite.h"
#include "TVTrackHit.h"

#include "EXVTXKalTest.h"
#include "EXVTXKalDetector.h"
#include "EXVTXHit.h"
#include "EXVTXEventGen.h"
#include "EXHYBTrack.h"

#include "TCanvas.h"
#include "TView.h"
#include "TRotMatrix.h"
#include "TNode.h"
#include "TString.h"
#include "TPolyMarker3D.h"

#include <iostream>
#include <iomanip>
#include <sstream>


/**************************************************************************
 * Changes:
 * F.Gaede, 10.11.2010 : added track parameters d0,tnl and errors^2 
 *                       d0err2,fi0err2,cpaerr2,dzerr2,tnlerr2 to ntuple
 **************************************************************************/

//FG: if this is active errors^2 for fi0, tnl and cpa are partly negative !?
//#define SAVE_RESIDUAL

static const Bool_t gkDir = kIterBackward;
//static const Bool_t gkDir = kIterForward;

using namespace std;

int main (Int_t argc, Char_t **argv)
{
   // ===================================================================
   //  Get job parameters from command line arguments, if any
   // ===================================================================

   Int_t    offset  = 0;
   if (argc > 1 && TString(argv[1]) == "-b") {
      offset = 1;
      gROOT->SetBatch();   // batch mode without event display
   }

   Double_t pt      = 10.;   // default Pt [GeV]
   Double_t t0in    = 14.;   // default tp [nsec]
   //Double_t cosmin  = -0.97; // default cos minimum [deg]
   Double_t cosmin  =  0; // default cos minimum [deg]
   Double_t cosmax  =  0; // default cos maximum [deg]
   Int_t    nevents = 1;     // default number of events to generate
   switch (argc-offset) {
      case 6: 
         nevents = atoi(argv[1+offset]);
         pt      = atof(argv[2+offset]);
         t0in    = atof(argv[3+offset]);
         cosmin  = atof(argv[4+offset]);
         cosmax  = atof(argv[5+offset]);
         break;
      case 5: 
         nevents = atoi(argv[1+offset]);
         pt      = atof(argv[2+offset]);
         t0in    = atof(argv[3+offset]);
         cosmin  = cosmax = atof(argv[4+offset]);
         break;
      case 4: 
         nevents = atoi(argv[1+offset]);
         pt      = atof(argv[2+offset]);
         t0in    = atof(argv[3+offset]);
         break;
      case 3: 
         nevents = atoi(argv[1+offset]);
         pt      = atof(argv[2+offset]);
         break;
      case 2:
         nevents = atoi(argv[1+offset]);
         break;
      case 1:
         break;
      default:
         cerr << "Too many command line arguments!" << endl;
         abort();
   }

   // ===================================================================
   //  Create TApplication
   // ===================================================================

   TApplication app("EXVTXKalTest", &argc, argv, 0, 0);

   // ===================================================================
   //  Prepare a detector
   // ===================================================================

   TKalDetCradle    toygld; // toy GLD detector
   EXVTXKalDetector vtxdet; // vertex detector (vtx)

   toygld.Install(vtxdet);  // install vtx into its toygld
   toygld.Close();          // close the cradle
   toygld.Sort();           // sort meas. layers from inside to outside

   //vtxdet.PowerOff();       // power off vtx not to process hit
   toygld.SwitchOffMS();    // switch off multiple scattering
   toygld.SwitchOffDEDX();  // switch off enery loss
   //toygld.SwitchOnMS();    // switch off multiple scattering
   //toygld.SwitchOnDEDX();  // switch off enery loss

   // ===================================================================
   //  Prepare an output n-tuple
   // ===================================================================

   TFile hfile("h.root","RECREATE","KalTest");

   stringstream sout;
   //   sout << "ndf:chi2:cl:fi0:cpa:cs:t0"; // 7 items
   sout << "ndf:chi2:cl:d0:fi0:cpa:dz:tnl:cs:t0:d0err2:fi0err2:cpaerr2:dzerr2:tnlerr2";  // 15 items

#ifdef SAVE_RESIDUAL
   //   Int_t nitems  = 7;
   Int_t nitems0  = 15;
   Int_t nitems  = nitems0 ;
   
   Int_t itemID[2000][3];
   TIter nextlayer(&toygld);
   EXVTXVMeasLayer *mlp;
   while ((mlp = dynamic_cast<EXVTXVMeasLayer *>(nextlayer()))) {
      EXVMeasLayer &ml = *mlp;
      if (ml.IsActive()) {
         Int_t index = ml.GetIndex();
         sout << ":dxin" << setw(3) << setfill('0') << index;
	 itemID[index][0] = nitems++;
         sout << ":dxot" << setw(3) << setfill('0') << index;
	 itemID[index][1] = nitems++;
         sout << ":z" << setw(3) << setfill('0') << index;
	 itemID[index][2] = nitems++;
#if 1
         cerr << "index = " << setw(4) << setfill(' ') << index << " "
              << "name = "  << ml.GetMLName() << endl;
#endif
      }
   }
   Double_t *data = new Double_t [nitems];
#endif
   sout << ends;
   TNtupleD *hTrackMonitor = new TNtupleD("track", "", sout.str().data());

   // ===================================================================
   //  Prepare an Event Generator
   // ===================================================================

   TObjArray kalhits;                // array to store hits
   EXVTXEventGen gen(toygld, kalhits);  // create event generator
   gen.SetT0(t0in);                  // set bunch crossing timing (t0)

   // ===================================================================
   //  Loop over events
   // ===================================================================

   for (Int_t eventno = 0; eventno < nevents; eventno++) { 
      cerr << "------ Event " << eventno << " ------" << endl;

      kalhits.Delete(); // clear hit data

      // ============================================================
      //  Generate a partcle
      // ============================================================

      THelicalTrack hel = gen.GenerateHelix(pt, cosmin, cosmax);

      // ============================================================
      //  Swim the particle in detector
      // ============================================================

      gen.Swim(hel);
	  std::vector<TVector3> hitVec;


	  // ------------------------- Dev ------------------------------

	  std::cout << "Number of Events read from LCIO: " << gen.LoadHits() << std::endl << std::endl;

	  // ----------------------- End Dev ----------------------------

      // ============================================================
      //  Do Kalman Filter
      // ============================================================

      cout << "nhits = " 
           << kalhits.GetEntries() << " >>>>>>>" << endl;

      if (kalhits.GetEntries() < 3) {
         cerr << "<<<<<< Shortage of Hits! nhits = " 
              << kalhits.GetEntries() << " >>>>>>>" << endl;
         continue;
      }

      Int_t i1, i2, i3; // (i1,i2,i3) = (1st,mid,last) hit to filter
      if (gkDir == kIterBackward) {
         i3 = 0;
         i1 = kalhits.GetEntries() - 1;
         i2 = i1 / 2;
      } else {
         i1 = 0;
         i3 = kalhits.GetEntries() - 1;
         i2 = i3 / 2;
      }

      // ---------------------------
      //  Create a dummy site: sited
      // ---------------------------

      TVTrackHit *ht1p = dynamic_cast<TVTrackHit *>(kalhits.At(i1));
      TVTrackHit *htdp = 0;
      if (dynamic_cast<EXVTXHit *>(ht1p)) {
         htdp = new EXVTXHit(*dynamic_cast<EXVTXHit *>(ht1p));
      } 

      TVTrackHit &hitd = *htdp;

      hitd(0,1) = 1.e6;   // give a huge error to d
      hitd(1,1) = 1.e6;   // give a huge error to z

      TKalTrackSite &sited = *new TKalTrackSite(hitd);
      sited.SetHitOwner();// site owns hit
      sited.SetOwner();   // site owns states

      // ---------------------------
      //  Create initial helix
      // ---------------------------

      TVTrackHit &h1 = *dynamic_cast<TVTrackHit *>(kalhits.At(i1)); // first hit
      TVTrackHit &h2 = *dynamic_cast<TVTrackHit *>(kalhits.At(i2)); // middle hit
      TVTrackHit &h3 = *dynamic_cast<TVTrackHit *>(kalhits.At(i3)); // last hit
      TVector3    x1 = h1.GetMeasLayer().HitToXv(h1);
      TVector3    x2 = h2.GetMeasLayer().HitToXv(h2);
      TVector3    x3 = h3.GetMeasLayer().HitToXv(h3);
      THelicalTrack helstart(x1, x2, x3, h1.GetBfield(), gkDir); // initial helix 

	  //cout << "Helix rho: " << helstart.GetRho() << endl;

      // ---------------------------
      //  Set dummy state to sited
      // ---------------------------

      static TKalMatrix svd(kSdim,1);
      svd(0,0) = 0.;                        // dr
      svd(1,0) = helstart.GetPhi0();        // phi0
      svd(2,0) = helstart.GetKappa();       // kappa
      svd(3,0) = 0.;                        // dz
      svd(4,0) = helstart.GetTanLambda();   // tan(lambda)

      if (kSdim == 6) svd(5,0) = 0.;        // t0

      static TKalMatrix C(kSdim,kSdim);

	  for (Int_t i=0; i<kSdim; i++) {
		C(i,i) = 1.e-4;
      }

      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
      sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

      // ---------------------------
      //  Add sited to the kaltrack
      // ---------------------------

      EXHYBTrack kaltrack;   // a track is a kal system
      kaltrack.SetOwner();   // kaltrack owns sites
      kaltrack.Add(&sited);  // add the dummy site to this track

      // ---------------------------
      //  Prepare hit iterrator
      // ---------------------------

      TIter next(&kalhits, gkDir); // come in to IP, if gkDir = kIterBackward

      // ---------------------------
      //  Start Kalman Filter
      // ---------------------------

      TVTrackHit *hitp = 0;
	  Int_t hitIndex = 0;

	  Int_t siteDiscarded = 0;
	  Int_t nsite = -1;
      while ((hitp = dynamic_cast<TVTrackHit *>(next()))) {
		  //cout << endl << "-------------site: " << ++nsite << endl;
	   	  TKalTrackSite  &site = *new TKalTrackSite(*hitp); // new site
		  hitVec.push_back(site.GetPivot());
		  
         if (!kaltrack.AddAndFilter(site)) {               // filter it
            delete &site;                        // delete it if failed
			++siteDiscarded;
			if(siteDiscarded>1) {
				cout << "The discarded site number is 2, end of tracking." << endl;
				break;
			}
			else
				continue;
         }

		 siteDiscarded = 0;

      } // end of Kalman filter

      // ---------------------------
      //  Smooth the track
      // ---------------------------
#ifndef SAVE_RESIDUAL
      TVKalSite &cursite = kaltrack.GetCurSite();
#else
      Int_t isite = 1;
      kaltrack.SmoothBackTo(isite);
      TVKalSite &cursite = static_cast<TVKalSite &>(*kaltrack[isite]);
#endif

      // ============================================================
      //  Monitor Fit Result
      // ============================================================

      Int_t    ndf  = kaltrack.GetNDF();
      Double_t chi2 = kaltrack.GetChi2();
      Double_t cl   = TMath::Prob(chi2, ndf);
      Double_t d0  =  cursite.GetCurState()(0, 0 ); 
      Double_t fi0  = cursite.GetCurState()(1, 0 ); 
      Double_t cpa  = cursite.GetCurState()(2, 0 ); 
      Double_t dz   = cursite.GetCurState()(3, 0 ); 
      Double_t tnl  = cursite.GetCurState()(4, 0 ); 
      Double_t cs   = tnl/TMath::Sqrt(1.+tnl*tnl );
      Double_t t0   = cursite.GetCurState()(5, 0 ); 

      const TKalMatrix& covK = cursite.GetCurState().GetCovMat() ; 

      // errors^2 of track parameters
      double d0err2  = covK( 0 , 0 )   ;
      double fi0err2 = covK( 1 , 1 )   ;
      double cpaerr2 = covK( 2 , 2 )   ;
      double dzerr2  = covK( 3 , 3 )   ;
      double tnlerr2 = covK( 4 , 4 )   ;


#ifndef SAVE_RESIDUAL
      hTrackMonitor->Fill(ndf, chi2, cl, d0, fi0, cpa, dz, tnl , cs, t0 , 
			  d0err2, fi0err2, cpaerr2, dzerr2, tnlerr2 ) ;

#else
      data[0] = ndf ;
      data[1] = chi2 ;
      data[2] = cl ;
      data[3] = d0 ;
      data[4] = fi0 ;
      data[5] = cpa ;
      data[6] = dz ;
      data[7] = tnl ;
      data[8] = cs ;
      data[9] = t0 ;
      data[10] = d0err2 ;
      data[11] = fi0err2 ;
      data[12] = cpaerr2 ;
      data[13] = dzerr2 ;
      data[14] = tnlerr2 ;

      for (Int_t i=nitems0; i<nitems; i++) data[i] = 9999999.;
      TIter nextsite(&kaltrack);
            nextsite(); // skip dummy site
      TKalTrackSite *sitep;
      while ((sitep = static_cast<TKalTrackSite *>(nextsite()))) {
         TKalTrackSite &site = *sitep;
         Int_t index = site.GetHit().GetMeasLayer().GetIndex();
	 data[itemID[index][0]] = site.GetResVec(TVKalSite::kSmoothed)(0,0);
         site.InvFilter();
         data[itemID[index][1]] = site.GetResVec(TVKalSite::kInvFiltered)(0,0);
	 data[itemID[index][2]] = site.GetPivot().Z();
      }
      hTrackMonitor->Fill(data);
#endif

      // ============================================================
      //  Very Primitive Event Display
      // ============================================================

      static TCanvas     *cvp    = 0;
      if (!gROOT->IsBatch()) {
         if (!cvp) {
            cvp = new TCanvas("OED", "Event Display", 400, 400, 610, 610);
         } else {
            cvp->cd();
            cvp->Clear();
         }

         TView   *vwp = TView::CreateView(1,0,0);
		 //const double WVWP = 70;
		 const double WVWP = 2500;
         vwp->SetRange(-WVWP,-WVWP,-WVWP,+WVWP,+WVWP,+WVWP);
         Int_t ierr;
         //vwp->SetView(10.,80.,80.,ierr);
         vwp->SetView(0.,0.,0.,ierr);

         vtxdet.Draw(40);
		 //kaltrack.Print();

		 //Draw all kalhits
		 //cout << "hit number: " << hitVec.size() << endl;
         if (hitVec.size()>0) {
			 //cout << ">>>>>>>>>>>>>>>>>draw all hits: " << hitVec.size() << endl;
            gPad->cd();

            TPolyMarker3D *pm3dp = new TPolyMarker3D(hitVec.size());
            pm3dp->SetBit(kCanDelete);
            pm3dp->SetMarkerColor(40);
            pm3dp->SetMarkerStyle(6);

		    for(int ihit=0; ihit<hitVec.size(); ++ihit) {	
			   TVector3& pos = hitVec[ihit];
               pm3dp->SetPoint(ihit, pos.X(), pos.Y(), pos.Z());
			   //cout << "pos: " << pos.X() << ", " << pos.Y() << ", " << pos.Z() <<  endl;
            }

            pm3dp->Draw();
            gPad->Update();
		 }

         kaltrack.Draw(2,"");         
		 // Dev: shortcut
		 app.Run(kTRUE);
		 // End Dev

         cout << "Next? [yes/no/edit/quit] " << flush;
         static const Int_t kMaxLen = 1024;
         Char_t temp[kMaxLen];
         cin.getline(temp,kMaxLen);
         TString opts(temp);
         opts.ToLower();
         if (!opts.Length()) {
            continue;
         } else if (opts[0] == 'n' || opts[0] == 'q') {
            break;
         } else if (opts[0] == 'e') {
            cout << "Select \"Quit ROOT\" from \"File\" to display next" << endl;
            cout << "\"CNTRL+C\" to really quit" << endl;
            app.Run(kTRUE);
         }
      } // endo fo event display
   } // end of event loop

   // ===================================================================
   //  Write results to file.
   // ===================================================================

   hfile.Write();

   return 0;
}
