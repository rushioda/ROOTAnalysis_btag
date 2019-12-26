#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <MyAnalysisExampleROOT/MyAnalysisAlg.h>
#include <TSystem.h>

// Infrastructure include(s): 
#include "xAODRootAccess/Init.h" 
#include "xAODRootAccess/TEvent.h" 
// ASG status code check 
#include <AsgTools/MessageCheck.h>

#include "xAODEventInfo/EventInfo.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODTracking/TrackParticle.h"

#include "AthLinks/ElementLink.h"
#include "AthLinks/ElementLinkVector.h"
#include <TFile.h>

  
// this is needed to distribute the algorithm to the workers                    
ClassImp(MyAnalysisAlg)

std::pair<unsigned int, unsigned int> res;

MyAnalysisAlg :: MyAnalysisAlg ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



EL::StatusCode MyAnalysisAlg :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  
  job.useXAOD (); 
  ANA_CHECK(xAOD::Init());

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyAnalysisAlg :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  // get the output file, create a new TTree and connect it to that output
  // define what braches will go in that tree
  TFile *outputFile = wk()->getOutputFile (outputName);
  tree = new TTree ("tree", "Ntuple variables collection");
  tree->SetDirectory (outputFile);
  tree->Branch("EventNumber", &m_EventNumber);	

  //Primary vertex
  tree->Branch("VtxType", &m_VtxType);
  tree->Branch("pVtxX", &m_pVtxX);
  tree->Branch("pVtxY", &m_pVtxY);
  tree->Branch("pVtxZ", &m_pVtxZ);
  tree->Branch("pVtxXErr", &m_pVtxXErr);
  tree->Branch("pVtxYErr", &m_pVtxYErr);
  tree->Branch("pVtxZErr", &m_pVtxZErr);

  //Jet
  tree->Branch("nJets", &m_nJets);	
  m_JetPt = new std::vector<double>;  
  tree->Branch("JetPt", &m_JetPt); 
  m_JetEta = new std::vector<double>;  
  tree->Branch("JetEta", &m_JetEta); 
  m_JetPhi = new std::vector<double>;  
  tree->Branch("JetPhi", &m_JetPhi); 
  m_JetType = new std::vector<int>;  
  tree->Branch("JetType", &m_JetType); 
  m_JetOrigin = new std::vector<int>;  
  tree->Branch("JetOrigin", &m_JetOrigin); 
  m_JetTruthID = new std::vector<int>;  
  tree->Branch("JetTruthID", &m_JetTruthID); 
  m_JVFCorr = new std::vector<float>;  
  tree->Branch("JVFCorr", &m_JVFCorr); 
  m_JVT = new std::vector<float>;  
  tree->Branch("JVT", &m_JVT); 
  m_nTrks_SV1 =  new std::vector<int>;
  tree->Branch("nTrks_SV1", &m_nTrks_SV1); 
  m_SV1pu =  new std::vector<double>;  
  tree->Branch("SV1pu", &m_SV1pu);
  m_SV1pc =  new std::vector<double>;  
  tree->Branch("SV1pc", &m_SV1pc);
  m_SV1pb =  new std::vector<double>;  
  tree->Branch("SV1pb", &m_SV1pb);
  m_SV1_mVtx =  new std::vector<float>;  
  tree->Branch("SV1_mVtx", &m_SV1_mVtx);
  m_SV1_EfracVtx =  new std::vector<float>;  
  tree->Branch("SV1_EfracVtx", &m_SV1_EfracVtx);
  m_SV1_nVtx2Trk =  new std::vector<int>;  
  tree->Branch("SV1_nVtx2Trk", &m_SV1_nVtx2Trk);
  m_SV1_dR =  new std::vector<float>;  
  tree->Branch("SV1_dR", &m_SV1_dR);

  m_nTrks_IP2D =  new std::vector<int>;
  tree->Branch("nTrks_IP2D", &m_nTrks_IP2D);
  m_IP2Dpu =  new std::vector<double>;
  tree->Branch("IP2Dpu", &m_IP2Dpu);
  m_IP2Dpc =  new std::vector<double>;
  tree->Branch("IP2Dpc", &m_IP2Dpc);
  m_IP2Dpb =  new std::vector<double>;
  tree->Branch("IP2Dpb", &m_IP2Dpb);
  m_IP2DLR =  new std::vector<double>;
  tree->Branch("IP2DLR", &m_IP2DLR);

  m_nTrks_IP3D =  new std::vector<int>;
  tree->Branch("nTrks_IP3D", &m_nTrks_IP3D);
  m_IP3Dpu =  new std::vector<double>;
  tree->Branch("IP3Dpu", &m_IP3Dpu);
  m_IP3Dpc =  new std::vector<double>;
  tree->Branch("IP3Dpc", &m_IP3Dpc);
  m_IP3Dpb =  new std::vector<double>;
  tree->Branch("IP3Dpb", &m_IP3Dpb);
  m_IP3DLR =  new std::vector<double>;
  tree->Branch("IP3DLR", &m_IP3DLR);

  m_JFpu =  new std::vector<double>;
  tree->Branch("JFpu", &m_JFpu);
  m_JFpc =  new std::vector<double>;
  tree->Branch("JFpc", &m_JFpc);
  m_JFpb =  new std::vector<double>;
  tree->Branch("JFpb", &m_JFpb);
  m_JF_mass =  new std::vector<float>;  
  tree->Branch("JF_mass", &m_JF_mass);
  m_JF_Efrac =  new std::vector<float>;  
  tree->Branch("JF_Efrac", &m_JF_Efrac);
  m_JF_nVtx2Trk =  new std::vector<int>;  
  tree->Branch("JF_nVtx2Trk", &m_SV1_nVtx2Trk);
  m_JF_n1Trk =  new std::vector<int>;  
  tree->Branch("JF_n1Trk", &m_JF_n1Trk);
  m_JF_nVtx =  new std::vector<int>;  
  tree->Branch("JF_nVtx", &m_JF_nVtx);
  m_JF_dR =  new std::vector<float>;  
  tree->Branch("JF_dR", &m_JF_dR);

  m_MV2c10 =  new std::vector<double>;
  tree->Branch("MV2c10", &m_MV2c10);
 
  m_DL1pu =  new std::vector<double>;
  tree->Branch("DL1pu", &m_DL1pu);
  m_DL1pc =  new std::vector<double>;
  tree->Branch("DL1pc", &m_DL1pc);
  m_DL1pb =  new std::vector<double>;
  tree->Branch("DL1pb", &m_DL1pb);

  //Track 
  m_Trackd0 = new std::vector<std::vector<double> >; 
  tree->Branch("Trackd0", &m_Trackd0); 
  m_Trackd0err = new std::vector<std::vector<double> >; 
  tree->Branch("Trackd0err", &m_Trackd0err); 
  m_Trackz0 = new std::vector<std::vector<double> >;  
  tree->Branch("Trackz0", &m_Trackz0);  
  m_Trackz0err = new std::vector<std::vector<double> >;  
  tree->Branch("Trackz0err", &m_Trackz0err);
  m_TrackPt = new std::vector<std::vector<double> >; 
  tree->Branch("TrackPt", &m_TrackPt); 
  m_TrackEta = new std::vector<std::vector<double> >; 
  tree->Branch("TrackEta", &m_TrackEta); 
  m_TrackPhi = new std::vector<std::vector<double> >; 
  tree->Branch("TrackPhi", &m_TrackPhi); 
  m_TrackCharge = new std::vector<std::vector<double> >; 
  tree->Branch("TrackCharge", &m_TrackCharge); 
  m_TrackSignd0 = new std::vector<std::vector<float> >; 
  tree->Branch("TrackSignd0", &m_TrackSignd0); 
  m_TrackSignz0 = new std::vector<std::vector<float> >; 
  tree->Branch("TrackSignz0", &m_TrackSignz0); 
  m_TrackSignd0sig = new std::vector<std::vector<float> >; 
  tree->Branch("TrackSignd0sig", &m_TrackSignd0sig); 
  m_TrackSignz0sig = new std::vector<std::vector<float> >; 
  tree->Branch("TrackSignz0sig", &m_TrackSignz0sig); 
  m_V0Track = new std::vector<std::vector<bool> >; 
  tree->Branch("V0Track", &m_V0Track); 
  m_nPixelHits = new std::vector<std::vector<unsigned char> >; 
  tree->Branch("nPixelHits", &m_nPixelHits); 
  m_nSCTHits = new std::vector<std::vector<unsigned char> >; 
  tree->Branch("nSCTHits", &m_nSCTHits); 
  m_nPixelDead = new std::vector<std::vector<unsigned char> >; 
  tree->Branch("nPixelDead", &m_nPixelDead); 
  m_nSCTDead = new std::vector<std::vector<unsigned char> >; 
  tree->Branch("nSCTDead", &m_nSCTDead); 
  m_nPixelHoles = new std::vector<std::vector<unsigned char> >; 
  tree->Branch("nPixelHoles", &m_nPixelHoles); 
  m_nSCTHoles = new std::vector<std::vector<unsigned char> >; 
  tree->Branch("nSCTHoles", &m_nSCTHoles); 
  m_nPixelShared = new std::vector<std::vector<unsigned char> >; 
  tree->Branch("nPixelShared", &m_nPixelShared); 
  m_nSCTShared = new std::vector<std::vector<unsigned char> >; 
  tree->Branch("nSCTShared", &m_nSCTShared); 
  m_nIBLHits = new std::vector<std::vector<unsigned char> >; 
  tree->Branch("nIBLHits", &m_nIBLHits); 
  m_nBLHits = new std::vector<std::vector<unsigned char> >; 
  tree->Branch("nBLHits", &m_nBLHits); 
 
  //Truth
  m_TruthPt = new std::vector<std::vector<double> >; 
  tree->Branch("TruthPt", &m_TruthPt); 
  m_TruthEta = new std::vector<std::vector<double> >; 
  tree->Branch("TruthEta", &m_TruthEta); 
  m_TruthPhi = new std::vector<std::vector<double> >; 
  tree->Branch("TruthPhi", &m_TruthPhi); 
  m_ID = new std::vector<std::vector<int> >;
  tree->Branch("ID", &m_ID);
  m_Truthd0 = new std::vector<std::vector<float> >; 
  tree->Branch("Truthd0", &m_Truthd0); 
  m_Truthz0 = new std::vector<std::vector<float> >;  
  tree->Branch("Truthz0", &m_Truthz0);  
 
/*
  m_nPixHit = new std::vector<unsigned char>;  
  tree->Branch("nPixHit", &m_nPixHit); 
  m_nSCTHit = new std::vector<unsigned char>;  
  tree->Branch("nSCTHit", &m_nSCTHit); 
  m_nTRTHit = new std::vector<unsigned char>;  
  tree->Branch("nTRTHit", &m_nTRTHit); 
  m_match = new std::vector<int>;  
  tree->Branch("match", &m_match); 

  tree->Branch("nTruth", &m_nTruth);	
  m_ID_truth = new std::vector<int>;  
  tree->Branch("ID_truth", &m_ID_truth); 
*/
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyAnalysisAlg :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyAnalysisAlg :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyAnalysisAlg :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  xAOD::TEvent* event = wk()->xaodEvent(); 
  Info("initialize()", "Number of events = %lli", event->getEntries() );

  // MC truth classifier	
  m_truthClass = new MCTruthClassifier("MCTruthClassifier_MyAnalysis");
  if (m_truthClass->initialize().isFailure())
  ANA_MSG_FATAL("Couldn't initialize MCTruthClassifier");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyAnalysisAlg :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  
  xAOD::TEvent* event = wk()->xaodEvent();
 
  m_JetPt->clear();  
  m_JetEta->clear();  
  m_JetPhi->clear();  
  m_JetType->clear();  
  m_JetOrigin->clear();  
  m_JetTruthID->clear();  
  m_JVFCorr->clear();
  m_JVT->clear();

  m_nTrks_SV1->clear(); 
  m_SV1pu->clear();
  m_SV1pc->clear();
  m_SV1pb->clear();
  m_SV1_mVtx->clear();
  m_SV1_EfracVtx->clear();
  m_SV1_nVtx2Trk->clear();
  m_SV1_dR->clear();

  m_nTrks_IP2D->clear(); 
  m_IP2Dpu->clear();
  m_IP2Dpc->clear();
  m_IP2Dpb->clear();
  m_IP2DLR->clear();

  m_nTrks_IP3D->clear(); 
  m_IP3Dpu->clear();
  m_IP3Dpc->clear();
  m_IP3Dpb->clear();
  m_IP3DLR->clear();

  m_JFpu->clear();
  m_JFpc->clear();
  m_JFpb->clear();
  m_JF_mass->clear();
  m_JF_Efrac->clear();
  m_JF_nVtx2Trk->clear();
  m_JF_n1Trk->clear();
  m_JF_nVtx->clear();
  m_JF_dR->clear();

  m_MV2c10->clear();
 
  m_DL1pu->clear();
  m_DL1pc->clear();
  m_DL1pb->clear();
 
  m_Trackd0->clear();
  m_Trackd0err->clear();
  m_Trackz0->clear();
  m_Trackz0err->clear();
  m_TrackPt->clear();
  m_TrackEta->clear();
  m_TrackPhi->clear();
  m_TrackCharge->clear();
  m_TrackSignd0->clear();
  m_TrackSignz0->clear();
  m_TrackSignd0sig->clear();
  m_TrackSignz0sig->clear();
  m_V0Track->clear();
  m_nPixelHits->clear();
  m_nSCTHits->clear();
  m_nPixelDead->clear();
  m_nSCTDead->clear();
  m_nPixelHoles->clear();
  m_nSCTHoles->clear();
  m_nPixelShared->clear();
  m_nSCTShared->clear();
  m_nIBLHits->clear();
  m_nBLHits->clear();
 
  m_TruthPt->clear();
  m_TruthEta->clear();
  m_TruthPhi->clear();
  m_ID->clear();
  m_Truthd0->clear();
  m_Truthz0->clear();

  // Event information
  const xAOD::EventInfo* eventInfo = 0;
  ANA_CHECK(event->retrieve(eventInfo, "EventInfo"));
  unsigned long runNumber  = eventInfo->runNumber();
  unsigned long eventNumber = eventInfo->eventNumber();
  Info("execute()", "Run = %lu : Event = %lu", runNumber, eventNumber);
  m_EventNumber = eventInfo->eventNumber();

  if(eventInfo){ //event loop
    
    // Vertex
    int numVtx = 0;
/*
    int VtxType;
    float pVtxX = -999.;
    float pVtxY = -999.;
    float pVtxZ = -999.;
    float pVtxXErr = -999.;
    float pVtxYErr = -999.;
    float pVtxZErr = -999.;
*/
    const xAOD::VertexContainer* primaryVertex = 0;
    ANA_CHECK(event->retrieve(primaryVertex,"PrimaryVertices"));
    for (xAOD::VertexContainer::const_iterator vtx_itr=primaryVertex->begin(); vtx_itr!=primaryVertex->end(); vtx_itr++) {
      if ((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx) {
        const Amg::MatrixX& pVtxCov = (*vtx_itr)->covariancePosition();
        m_VtxType = (*vtx_itr)->type();
        m_pVtxX = (*vtx_itr)->x();
        m_pVtxY = (*vtx_itr)->y();
        m_pVtxZ = (*vtx_itr)->z();
        m_pVtxXErr = TMath::Sqrt(pVtxCov(0,0));
        m_pVtxYErr = TMath::Sqrt(pVtxCov(1,1));
        m_pVtxZErr = TMath::Sqrt(pVtxCov(2,2));
//        Info("execute()", " vertex, (x, y, z) = %.2f mm, %.2f mm, %.2f mm ", m_pVtxX, m_pVtxY, m_pVtxZ);
      }
      if ((*vtx_itr)->vertexType()!=0) { numVtx++; }
    }


    // Jet
    int nJets = 0;
    const xAOD::JetContainer* jets = nullptr;
    ANA_CHECK(event->retrieve(jets, "AntiKt4EMTopoJets"));
    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    for ( ; jet_itr != jet_end; ++jet_itr ){ //jet loop
//      Info("execute()", "jet pt = %.2f GeV", ((*jet_itr)->pt()*0.001));
      const xAOD::Jet* thisjet = (*jet_itr);
      res = m_truthClass->particleTruthClassifier(thisjet,true); // for DR matching
//      Info("execute()", "jet type %i/jet origin %i ", res.first, res.second);

      // get b-tag link
      const xAOD::BTagging* btag = thisjet->btagging();	
      if(!btag) continue;
//      Info("execute()", "ID = %i, IP3D_loglikelihoodratio  = %.2f, log(IP3Dpb/IP3Dpu) = %.2f", thisjet->auxdata<int>("ConeTruthLabelID"), btag->IP3D_loglikelihoodratio(), log(btag->IP3D_pb()/btag->IP3D_pu()));

      int nTrack_SV1 = btag->nSV1_TrackParticles();
      int nTrack_IP2D = btag->nIP2D_TrackParticles();
      int nTrack_IP3D = btag->nIP3D_TrackParticles();
     
      //get track links of IP3D
      std::vector<double> trackd0; 
      std::vector<double> trackd0err;
      std::vector<double> trackz0;
      std::vector<double> trackz0err;
      std::vector<double> trackPt;
      std::vector<double> trackEta;
      std::vector<double> trackPhi;
      std::vector<double> trackCharge;
      std::vector<float> trackSignd0;
      std::vector<float> trackSignz0;
      std::vector<float> trackSignd0sig;
      std::vector<float> trackSignz0sig;
      std::vector<bool> V0track;
      std::vector<unsigned char> numPixelHits;
      std::vector<unsigned char> numSCTHits;
      std::vector<unsigned char> numPixelDead;
      std::vector<unsigned char> numSCTDead;
      std::vector<unsigned char> numPixelHoles;
      std::vector<unsigned char> numSCTHoles;
      std::vector<unsigned char> numPixelShared;
      std::vector<unsigned char> numSCTShared;
      std::vector<unsigned char> numIBLHits;
      std::vector<unsigned char> numBLHits;

      std::vector<double> truthPt;
      std::vector<double> truthEta;
      std::vector<double> truthPhi;
      std::vector<int> pdgID;
      std::vector<float> truthd0;
      std::vector<float> truthz0;

      if(!nTrack_IP3D) continue;
      for(int iTrk=0; iTrk<nTrack_IP3D; iTrk++){ //track loop
        const xAOD::TrackParticle* trk_IP3D = btag->IP3D_TrackParticle(iTrk);
        m_truthClass->particleTruthClassifier(trk_IP3D); 
//        Info("execute()", "track d0 = %.2f mm, z0 = %.2f mm", (trk_IP3D->d0()), trk_IP3D->z0());

//-----test----------
        trackd0.push_back(trk_IP3D->d0());
        trackd0err.push_back(trk_IP3D->definingParametersCovMatrixVec().at(0));
        trackz0.push_back(trk_IP3D->z0());
        trackz0err.push_back(trk_IP3D->definingParametersCovMatrixVec().at(2));
        trackPt.push_back(trk_IP3D->pt()*0.001);
        trackEta.push_back(trk_IP3D->eta());
        trackPhi.push_back(trk_IP3D->phi());
        trackCharge.push_back(trk_IP3D->charge());

        trackSignd0 = btag->auxdata<std::vector<float>>("IP3D_valD0wrtPVofTracks");
        trackSignz0 = btag->auxdata<std::vector<float>>("IP3D_valZ0wrtPVofTracks");
        trackSignd0sig = btag->auxdata<std::vector<float>>("IP3D_sigD0wrtPVofTracks");
        trackSignz0sig = btag->auxdata<std::vector<float>>("IP3D_sigZ0wrtPVofTracks");
        V0track = btag->auxdata<std::vector<bool>>("IP3D_flagFromV0ofTracks");
        numPixelHits.push_back(trk_IP3D->auxdata<unsigned char>("numberOfPixelHits"));
        numSCTHits.push_back(trk_IP3D->auxdata<unsigned char>("numberOfSCTHits"));
        numPixelDead.push_back(trk_IP3D->auxdata<unsigned char>("numberOfPixelDeadSensors"));
        numSCTDead.push_back(trk_IP3D->auxdata<unsigned char>("numberOfSCTDeadSensors"));
        numPixelHoles.push_back(trk_IP3D->auxdata<unsigned char>("numberOfPixelHoles"));
        numSCTHoles.push_back(trk_IP3D->auxdata<unsigned char>("numberOfSCTHoles"));
        numPixelShared.push_back(trk_IP3D->auxdata<unsigned char>("numberOfPixelSharedHits"));
        numSCTShared.push_back(trk_IP3D->auxdata<unsigned char>("numberOfSCTSharedHits"));
        numIBLHits.push_back(trk_IP3D->auxdata<unsigned char>("numberOfInnermostPixelLayerHits"));
        numBLHits.push_back(trk_IP3D->auxdata<unsigned char>("numberOfNextToInnermostPixelLayerHits"));

//-----test----------
        
        /* 
        // get vertex
        const xAOD::Vertex* vtx = trk_IP3D->vertex();
//        if(!vtx) continue;
*/
	      // get truth particles
        const xAOD::TruthParticle* tru_IP3D =  m_truthClass->getGenPart();
        if(!tru_IP3D){
          truthPt.push_back(-999);
          truthEta.push_back(-999);
          truthPhi.push_back(-999);
          pdgID.push_back(-999);
          truthd0.push_back(-999);
          truthz0.push_back(-999);
          continue;
        }
//          Info("execute()", "truth pt = %.2f GeV" , tru_IP3D->pt()*0.001);
  
        truthPt.push_back(tru_IP3D->pt()*0.001);
        truthEta.push_back(tru_IP3D->eta());
        truthPhi.push_back(tru_IP3D->phi());
        pdgID.push_back(tru_IP3D->pdgId());
        truthd0.push_back(tru_IP3D->auxdata<float>("d0"));
        truthz0.push_back(tru_IP3D->auxdata<float>("z0"));

      } //end of track loop

      //Param def.
      double jetPt = thisjet->pt()*0.001;
      double jetEta = thisjet->eta();
      double jetPhi = thisjet->phi();
      double jetType = res.first;
      double jetOrigin = res.second;
      int jetTruthID = thisjet->auxdata<int>("ConeTruthLabelID");
      float jvf_corr = thisjet->auxdata<float>("JVFCorr");
      float jvt = thisjet->auxdata<float>("Jvt");

      double SV1_pu = btag->SV1_pu();
      double SV1_pc = btag->SV1_pc();
      double SV1_pb = btag->SV1_pb();
      float SV1_mVtx = btag->auxdata<float>("SV1_masssvx");
      float SV1_EfracVtx = btag->auxdata<float>("SV1_efracsvx");
      int SV1_nVtx2Trk = btag->auxdata<int>("SV1_N2Tpair");
      float SV1_dR = btag->auxdata<float>("SV1_deltaR");

      double IP2D_pu = btag->IP2D_pu();
      double IP2D_pc = btag->IP2D_pc();
      double IP2D_pb = btag->IP2D_pb();
      double IP2D_LR = btag->IP2D_loglikelihoodratio();
      double IP3D_pu = btag->IP3D_pu();
      double IP3D_pc = btag->IP3D_pc();
      double IP3D_pb = btag->IP3D_pb();
      double IP3D_LR = btag->IP3D_loglikelihoodratio();

      double JF_pu = btag->JetFitter_pu();
      double JF_pc = btag->JetFitter_pc();
      double JF_pb = btag->JetFitter_pb();
      float JF_mass = btag->auxdata<float>("JetFitter_mass");
      float JF_Efrac = btag->auxdata<float>("JetFitter_energyFraction");
      int JF_nVtx2Trk = btag->auxdata<int>("JetFitter_N2Tpair");
      int JF_n1Trk = btag->auxdata<int>("JetFitter_nSingleTracks");
      int JF_nVtx = btag->auxdata<int>("JetFitter_nVTX");
      float JF_dR = btag->auxdata<float>("JetFitter_dRFlightDir");

      double MV2c10 = btag->auxdata<double>("MV2c10_discriminant");

      double DL1_pu = btag->auxdata<double>("DL1_pu");
      double DL1_pc = btag->auxdata<double>("DL1_pc");
      double DL1_pb = btag->auxdata<double>("DL1_pb");
 
//      double MV1 = btag->MV1_discriminant();
//      std::cout << "Jet Type = " << jetType << ", MV1 score = " << MV1 << std::endl;

      //Fill
      m_JetPt->push_back(jetPt);
      m_JetEta->push_back(jetEta);
      m_JetPhi->push_back(jetPhi);
      m_JetType->push_back(jetType);
      m_JetOrigin->push_back(jetOrigin);
      m_JetTruthID->push_back(jetTruthID);
      m_JVFCorr->push_back(jvf_corr);
      m_JVT->push_back(jvt);

      m_nTrks_SV1->push_back(nTrack_SV1);
      m_SV1pu->push_back(SV1_pu);
      m_SV1pc->push_back(SV1_pc);
      m_SV1pb->push_back(SV1_pb);
      m_SV1_mVtx->push_back(SV1_mVtx);
      m_SV1_EfracVtx->push_back(SV1_EfracVtx);
      m_SV1_nVtx2Trk->push_back(SV1_nVtx2Trk);
      m_SV1_dR->push_back(SV1_dR);

      m_nTrks_IP2D->push_back(nTrack_IP2D);
      m_IP2Dpu->push_back(IP2D_pu); 
      m_IP2Dpc->push_back(IP2D_pc); 
      m_IP2Dpb->push_back(IP2D_pb); 
      m_IP2DLR->push_back(IP2D_LR); 

      m_nTrks_IP3D->push_back(nTrack_IP3D);
      m_IP3Dpu->push_back(IP3D_pu); 
      m_IP3Dpc->push_back(IP3D_pc); 
      m_IP3Dpb->push_back(IP3D_pb); 
      m_IP3DLR->push_back(IP3D_LR); 

      m_JFpu->push_back(JF_pu); 
      m_JFpc->push_back(JF_pc); 
      m_JFpb->push_back(JF_pb); 
      m_JF_mass->push_back(JF_mass);
      m_JF_Efrac->push_back(JF_Efrac);
      m_JF_nVtx2Trk->push_back(JF_nVtx2Trk);
      m_JF_n1Trk->push_back(JF_n1Trk);
      m_JF_nVtx->push_back(JF_nVtx);
      m_JF_dR->push_back(JF_dR);
 
      m_MV2c10->push_back(MV2c10); 

      m_DL1pu->push_back(DL1_pu); 
      m_DL1pc->push_back(DL1_pc); 
      m_DL1pb->push_back(DL1_pb); 

      m_Trackd0->push_back(trackd0); 
      m_Trackd0err->push_back(trackd0err); 
      m_Trackz0->push_back(trackz0); 
      m_Trackz0err->push_back(trackz0err); 
      m_TrackPt->push_back(trackPt); 
      m_TrackEta->push_back(trackEta); 
      m_TrackPhi->push_back(trackPhi); 
      m_TrackCharge->push_back(trackCharge); 
      m_TrackSignd0->push_back(trackSignd0); 
      m_TrackSignz0->push_back(trackSignz0); 
      m_TrackSignd0sig->push_back(trackSignd0sig); 
      m_TrackSignz0sig->push_back(trackSignz0sig); 
      m_V0Track->push_back(V0track);
      m_nPixelHits->push_back(numPixelHits); 
      m_nSCTHits->push_back(numSCTHits); 
      m_nPixelDead->push_back(numPixelDead); 
      m_nSCTDead->push_back(numSCTDead); 
      m_nPixelHoles->push_back(numPixelHoles); 
      m_nSCTHoles->push_back(numSCTHoles); 
      m_nPixelShared->push_back(numPixelShared); 
      m_nSCTShared->push_back(numSCTShared); 
      m_nIBLHits->push_back(numIBLHits); 
      m_nBLHits->push_back(numBLHits); 

      m_TruthPt->push_back(truthPt); 
      m_TruthEta->push_back(truthEta); 
      m_TruthPhi->push_back(truthPhi); 
      m_ID->push_back(pdgID); 
      m_Truthd0->push_back(truthd0); 
      m_Truthz0->push_back(truthz0); 
 
      trackd0.clear();
      trackd0err.clear();
      trackz0.clear();
      trackz0err.clear();
      trackPt.clear();
      trackEta.clear();
      trackPhi.clear();
      trackCharge.clear();
      trackSignd0.clear(); 
      trackSignz0.clear(); 
      trackSignd0sig.clear(); 
      trackSignz0sig.clear(); 
      V0track.clear();
      numPixelHits.clear(); 
      numSCTHits.clear(); 
      numPixelDead.clear(); 
      numSCTDead.clear(); 
      numPixelHoles.clear(); 
      numSCTHoles.clear(); 
      numPixelShared.clear(); 
      numSCTShared.clear(); 
      numIBLHits.clear(); 
      numBLHits.clear(); 

      truthPt.clear();
      truthEta.clear();
      truthPhi.clear();
      pdgID.clear();
      truthd0.clear();
      truthz0.clear();
 
      nJets++;
    } //end of jet loop 
/*
 	int hoge = 0;
	for (auto trackd0 : *m_Trackd0) {
    std::cout << hoge << ": ";
	  for (auto d0 : trackd0) {
      std::cout << d0 << ", ";
		}
    std::cout << std::endl;
		hoge++;
	}
*/

    m_nJets = static_cast<int>(nJets);	
    tree->Fill();

  } //end of event loop

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyAnalysisAlg :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyAnalysisAlg :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  delete m_JetPt; m_JetPt = 0;
  delete m_JetEta; m_JetEta = 0;
  delete m_JetPhi; m_JetPhi  = 0;
  delete m_JetType; m_JetType = 0;  
  delete m_JetOrigin; m_JetOrigin = 0;  
  delete m_JetTruthID; m_JetTruthID = 0;  
  delete m_JVFCorr; m_JVFCorr = 0;
  delete m_JVT; m_JVT = 0;

  delete m_nTrks_SV1; m_nTrks_SV1 = 0;
  delete m_SV1pu; m_SV1pu = 0;
  delete m_SV1pc; m_SV1pc = 0;
  delete m_SV1pb; m_SV1pb = 0;
  delete m_SV1_mVtx; m_SV1_mVtx = 0;
  delete m_SV1_EfracVtx; m_SV1_EfracVtx = 0;
  delete m_SV1_nVtx2Trk; m_SV1_nVtx2Trk = 0;
  delete m_SV1_dR; m_SV1_dR = 0;

  delete m_nTrks_IP2D; m_nTrks_IP2D = 0;
  delete m_IP2Dpu; m_IP2Dpu = 0; 
  delete m_IP2Dpc; m_IP2Dpc = 0; 
  delete m_IP2Dpb; m_IP2Dpb = 0; 
  delete m_IP2DLR; m_IP2DLR = 0; 

  delete m_nTrks_IP3D; m_nTrks_IP3D = 0;
  delete m_IP3Dpu; m_IP3Dpu = 0; 
  delete m_IP3Dpc; m_IP3Dpc = 0; 
  delete m_IP3Dpb; m_IP3Dpb = 0; 
  delete m_IP3DLR; m_IP3DLR = 0; 

  delete m_JFpu; m_JFpu = 0; 
  delete m_JFpc; m_JFpc = 0; 
  delete m_JFpb; m_JFpb = 0; 
  delete m_JF_mass; m_JF_mass = 0;
  delete m_JF_Efrac; m_JF_Efrac = 0;
  delete m_JF_nVtx2Trk; m_JF_nVtx2Trk = 0;
  delete m_JF_n1Trk; m_JF_n1Trk = 0;
  delete m_JF_nVtx; m_JF_nVtx = 0;
  delete m_JF_dR; m_JF_dR = 0;
 
  delete m_MV2c10; m_MV2c10 = 0; 

  delete m_DL1pu; m_DL1pu = 0; 
  delete m_DL1pc; m_DL1pc = 0; 
  delete m_DL1pb; m_DL1pb = 0; 

  delete m_Trackd0; m_Trackd0 = 0; 
  delete m_Trackd0err; m_Trackd0err = 0; 
  delete m_Trackz0; m_Trackz0 = 0; 
  delete m_Trackz0err; m_Trackz0err = 0; 
  delete m_TrackPt; m_TrackPt = 0; 
  delete m_TrackEta; m_TrackEta = 0; 
  delete m_TrackPhi; m_TrackPhi = 0; 
  delete m_TrackCharge; m_TrackCharge = 0; 
  delete m_TrackSignd0; m_TrackSignd0 = 0; 
  delete m_TrackSignz0; m_TrackSignz0 = 0; 
  delete m_TrackSignd0sig; m_TrackSignd0sig = 0; 
  delete m_TrackSignz0sig; m_TrackSignz0sig = 0;
  delete m_V0Track; m_V0Track = 0;
  delete m_nPixelHits; m_nPixelHits = 0;
  delete m_nSCTHits; m_nSCTHits = 0;
  delete m_nPixelDead; m_nPixelDead = 0;
  delete m_nSCTDead; m_nSCTDead = 0;
  delete m_nPixelHoles; m_nPixelHoles = 0;
  delete m_nSCTHoles; m_nSCTHoles = 0;
  delete m_nPixelShared; m_nPixelShared = 0;
  delete m_nSCTShared; m_nSCTShared = 0;
  delete m_nIBLHits; m_nIBLHits = 0;
  delete m_nBLHits; m_nBLHits = 0;

  delete m_TruthPt; m_TruthPt = 0; 
  delete m_TruthEta; m_TruthEta = 0; 
  delete m_TruthPhi; m_TruthPhi = 0; 
  delete m_ID; m_ID = 0; 
  delete m_Truthd0; m_Truthd0 = 0;
  delete m_Truthz0; m_Truthz0 = 0;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyAnalysisAlg :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}
