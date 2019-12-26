#ifndef MyAnalysisExampleROOT_MyAnalysisAlg_H
#define MyAnalysisExampleROOT_MyAnalysisAlg_H
#include "MCTruthClassifier/MCTruthClassifier.h"
#include <TTree.h>

#include <EventLoop/Algorithm.h>

class MyAnalysisAlg : public EL::Algorithm
{
   MCTruthClassifier* m_truthClass; //!
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)

public:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  TTree *tree; //!
  int m_EventNumber; //!
  int m_VtxType; //!
  double m_pVtxX; //!
  double m_pVtxY; //!
  double m_pVtxZ; //!
  double m_pVtxXErr; //!
  double m_pVtxYErr; //!
  double m_pVtxZErr; //!

  int m_nJets; //!
  std::vector<double> *m_JetPt; //!
  std::vector<double> *m_JetEta; //!
  std::vector<double> *m_JetPhi; //!
  std::vector<int> *m_JetType; //!
  std::vector<int> *m_JetOrigin; //!
  std::vector<int> *m_JetTruthID; //!
  std::vector<float> *m_JVFCorr; //!
  std::vector<float> *m_JVT; //!

  std::vector<int> *m_nTrks_SV1; //!
  std::vector<double> *m_SV1pu; //!
  std::vector<double> *m_SV1pc; //!
  std::vector<double> *m_SV1pb; //!
  std::vector<float> *m_SV1_mVtx; //! 
  std::vector<float> *m_SV1_EfracVtx; //! 
  std::vector<int> *m_SV1_nVtx2Trk; //!
  std::vector<float> *m_SV1_dR; //! 

  std::vector<int> *m_nTrks_IP2D; //!
  std::vector<double> *m_IP2Dpu; //!
  std::vector<double> *m_IP2Dpc; //!
  std::vector<double> *m_IP2Dpb; //!
  std::vector<double> *m_IP2DLR; //!

  std::vector<int> *m_nTrks_IP3D; //!
  std::vector<double> *m_IP3Dpu; //!
  std::vector<double> *m_IP3Dpc; //!
  std::vector<double> *m_IP3Dpb; //!
  std::vector<double> *m_IP3DLR; //!

  std::vector<double> *m_JFpu; //!
  std::vector<double> *m_JFpc; //!
  std::vector<double> *m_JFpb; //!
  std::vector<float> *m_JF_mass; //! 
  std::vector<float> *m_JF_Efrac; //! 
  std::vector<int> *m_JF_nVtx2Trk; //!
  std::vector<int> *m_JF_n1Trk; //!
  std::vector<int> *m_JF_nVtx; //!
  std::vector<float> *m_JF_dR; //! 

  std::vector<double> *m_MV2c10; //!

  std::vector<double> *m_DL1pu; //!
  std::vector<double> *m_DL1pc; //!
  std::vector<double> *m_DL1pb; //!


  std::vector<std::vector<double> > *m_Trackd0; //!
  std::vector<std::vector<double> > *m_Trackd0err; //!
  std::vector<std::vector<double> > *m_Trackz0; //!
  std::vector<std::vector<double> > *m_Trackz0err; //!
  std::vector<std::vector<double> > *m_TrackPt; //!
  std::vector<std::vector<double> > *m_TrackEta; //!
  std::vector<std::vector<double> > *m_TrackPhi; //!
  std::vector<std::vector<double> > *m_TrackCharge; //!
  std::vector<std::vector<float> > *m_TrackSignd0; //!
  std::vector<std::vector<float> > *m_TrackSignz0; //!
  std::vector<std::vector<float> > *m_TrackSignd0sig; //!
  std::vector<std::vector<float> > *m_TrackSignz0sig; //!
  std::vector<std::vector<bool> > *m_V0Track; //!
  std::vector<std::vector<unsigned char> > *m_nPixelHits; //!
  std::vector<std::vector<unsigned char> > *m_nSCTHits; //!
  std::vector<std::vector<unsigned char> > *m_nPixelHoles; //!
  std::vector<std::vector<unsigned char> > *m_nSCTHoles; //!
  std::vector<std::vector<unsigned char> > *m_nPixelDead; //!
  std::vector<std::vector<unsigned char> > *m_nSCTDead; //!
  std::vector<std::vector<unsigned char> > *m_nPixelShared; //!
  std::vector<std::vector<unsigned char> > *m_nSCTShared; //!
  std::vector<std::vector<unsigned char> > *m_nIBLHits; //!
  std::vector<std::vector<unsigned char> > *m_nBLHits; //!
  
  std::vector<std::vector<double> > *m_TruthPt; //!
  std::vector<std::vector<double> > *m_TruthEta; //!
  std::vector<std::vector<double> > *m_TruthPhi; //!
  std::vector<std::vector<int> > *m_ID; //!
  std::vector<std::vector<float> > *m_Truthd0; //!
  std::vector<std::vector<float> > *m_Truthz0; //!

  std::string outputName;

  // this is a standard constructor
  MyAnalysisAlg ();

 
  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers                  
  ClassDef(MyAnalysisAlg, 1);

};

#endif
