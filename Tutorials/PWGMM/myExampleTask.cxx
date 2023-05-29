/// \author Koki Soeda
/// \since 10/05/2023

#include <iostream>
#include <cmath>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "TDatabasePDG.h"
#include "MathUtils/Utils.h"

using namespace std;
using namespace o2;
using namespace o2::framework;

struct MyAnalysisTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsVer{"nBinsVer", 300, "N bins in Collisino Vertex histo"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsXY{"nBinsXY", 10000, "N bins in X and Y axis"};
  Configurable<int> nBinsZ{"nBinsZ", 5000, "N bins in Z axis"};
  Configurable<int> nBinsDCA{"nBinsDCA", 5000, "N bins in DCA"};
  Configurable<int> nBinsDCAT{"nBinsDCAT", 1000, "N bins in DCAT"};
  Configurable<int> nBinsDist{"nBinsDist", 500, ""};
  Configurable<int> fwdTrackType{"fwdTrackType", 0, "N TrackType in fwd"};
  Configurable<int> fwdTrackType2{"fwdTrackType2", 2, "N TrackType in fwd2"};

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisVertex{nBinsVer, -30, +30, "cm"};
    const AxisSpec axistrackType{6, -0.5, 5.5, "TrackType"};
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisChi2(505, -4.5, +500.5, "");
    const AxisSpec axisPDG(19, -0.5, +18.5, "");
    const AxisSpec axisX(nBinsXY, -10, 10, "x(cm)");
    const AxisSpec axisY(nBinsXY, -10, 10, "y(cm)");
    const AxisSpec axisZ(nBinsZ, -50, 50, "z(cm)");
    const AxisSpec axisDist(nBinsDist, -0.1, 50, "z(cm)");
    const AxisSpec axisDCAX(nBinsDCA, -1, 1, "x(cm)");
    const AxisSpec axisDCAY(nBinsDCA, -1, 1, "y(cm)");
    const AxisSpec axisDCAT(nBinsDCAT, -0.01, +1, "cm");
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}"};

    // create histograms
    histos.add("CollisionPos", "CollisionPos", kTH1F, {axisVertex});
    histos.add("TrackType", "TrackType", kTH1F, {axistrackType});
    histos.add("FwdParticleCounter", "FwdParticleCounter", kTH1F, {axisCounter});
    histos.add("FwdMuonCounter", "FwdMuonCounter", kTH1F, {axisCounter});
    histos.add("GMT_pT", "GMT_pT", kTH1F, {axisPt});
    histos.add("GlobalMuonTrackPDG", "GlobalMuonTrackPDG", kTH1F, {axisPDG});
    histos.add("ForwardMuonTrackPDG", "ForwardMuonTrackPDG", kTH1F, {axisPDG});
    histos.add("Chi2GMT", "Chi2GMT", kTH1F, {axisChi2});
    histos.add("muonMomPDG", "muonMomPDG", kTH1F, {axisPDG});
    histos.add("pi2mu_DCA_X", "DCA_X", kTH1F, {axisDCAX});
    histos.add("pi2mu_DCA_Y", "DCA_Y", kTH1F, {axisDCAY});
    histos.add("pi2mu_DCAT", "pi2mu_DCAT", kTH1F, {axisDCAT});
    histos.add("D0DecayX", "D0DecayX", kTH1F, {axisX});
    histos.add("D0DecayY", "D0DecayY", kTH1F, {axisY});
    histos.add("D0DecayZ", "D0DecayZ", kTH1F, {axisZ});
    histos.add("Distance_D0_z", "Distance_D0_z", kTH1F, {axisDist});
    histos.add("D02mu_DCA_X", "DCA_X", kTH1F, {axisDCAX});
    histos.add("D02mu_DCA_Y", "DCA_Y", kTH1F, {axisDCAY});
    histos.add("D02mu_DCAT", "D02mu_DCAT", kTH1F, {axisDCAT});
    histos.add("D0_pT", "D0_pT", kTH1F, {axisPt});
    histos.add("MFT_Kaon_pT", "MFT_Kaon_pT", kTH1F, {axisPt});
    histos.add("KpDecayX", "KpDecayX", kTH1F, {axisX});
    histos.add("KpDecayY", "KpDecayY", kTH1F, {axisY});
    histos.add("KpDecayZ", "KpDecayZ", kTH1F, {axisZ});
    histos.add("Kp2mu_DCA_X", "DCA_X", kTH1F, {axisDCAX});
    histos.add("Kp2mu_DCA_Y", "DCA_Y", kTH1F, {axisDCAY});
    histos.add("Kp2mu_DCAT", "Kp2mu_DCAT", kTH1F, {axisDCAT});
    histos.add("MFTParticleCounter", "MFTParticleCounter", kTH1F, {axisCounter});
    histos.add("MFTMuonCounter", "MFTMuonCounter", kTH1F, {axisCounter});
    histos.add("MFT_muon_pT", "MFT_muon_pT", kTH1F, {axisPt});
    histos.add("counter1", "counter1", kTH1F, {axisCounter});
    histos.add("counter2", "counter2", kTH1F, {axisCounter});

    auto hmpdgcode = histos.get<TH1>(HIST("muonMomPDG"));
    auto* x1 = hmpdgcode->GetXaxis();
    x1->SetBinLabel(1, "Primary"); // 1 - 37
    x1->SetBinLabel(2, "Pion0"); // 111
    x1->SetBinLabel(3, "Pion±"); // 211
    x1->SetBinLabel(4, "pho"); // 113
    x1->SetBinLabel(5, "K_l"); // 130
    x1->SetBinLabel(6, "Eta"); // 221
    x1->SetBinLabel(7, "Omega"); // 223
    x1->SetBinLabel(8, "K_s"); // 310
    x1->SetBinLabel(9, "K*0(892)"); // 313
    x1->SetBinLabel(10, "K±"); // 321
    x1->SetBinLabel(11, "K*±(892)"); // 323
    x1->SetBinLabel(12, "Eta_prim"); // 331
    x1->SetBinLabel(13, "Phi"); // 333
    x1->SetBinLabel(14, "D±"); // 411
    x1->SetBinLabel(15, "D*±"); // 413
    x1->SetBinLabel(16, "D0"); // 421
    x1->SetBinLabel(17, "D_s±"); // 431
    x1->SetBinLabel(18, "Beauty"); // 500-599
    x1->SetBinLabel(19, "Baryon"); // 1000 -
    
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
               soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
               soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
               aod::McParticles const& MCparticle,
               aod::McCollisions const&)
  { 
    if(collision.has_mcCollision()){
      //Informatin of MC collision
      auto mcCol_z = collision.mcCollision().posZ();
      histos.fill(HIST("CollisionPos"), mcCol_z);

      //Information of FwdTracks
      for(auto& fwdtrack : fwdtracks){
        histos.fill(HIST("TrackType"), fwdtrack.trackType());

        if(!fwdtrack.has_mcParticle()) continue;
        auto mcParticle_fwd = fwdtrack.mcParticle();
        histos.fill(HIST("FwdParticleCounter"), 0.5);
        auto fwdTrackPDG = mcParticle_fwd.pdgCode();

        if(fabs(fwdTrackPDG)==13) histos.fill(HIST("FwdMuonCounter"), 0.5);

        //Informaion of D0 meson
        float mudcaX, mudcaY, mumftX, mumftY, mumftZ;
        if(fwdtrack.trackType()==fwdTrackType || fwdtrack.trackType()==fwdTrackType2){//Required MFT-MCH-MID (GlobalMuonTrack)
          auto chi2GMT = fwdtrack.chi2MatchMCHMFT();
          auto GMT_pT = mcParticle_fwd.pt();
          histos.fill(HIST("GlobalMuonTrackPDG"), fwdTrackPDG);
          histos.fill(HIST("Chi2GMT"), chi2GMT);
          histos.fill(HIST("GMT_pT"), GMT_pT);

          if(chi2GMT==-1) continue;
          if(fabs(fwdTrackPDG)==13){
            auto a_ID = mcParticle_fwd.globalIndex();

            if(mcParticle_fwd.has_mothers()){
              auto mcMom = MCparticle.rawIteratorAt(mcParticle_fwd.mothersIds()[0]);
              auto mcMomPDG = mcMom.pdgCode();
              auto Daughters = mcMom.daughters_as<aod::McParticles>();
              
              if(fabs(mcMomPDG) < 38) histos.fill(HIST("muonMomPDG"), 0.0, 1);
              if(fabs(mcMomPDG) == 111) histos.fill(HIST("muonMomPDG"), 1.0, 1);
              if(fabs(mcMomPDG) == 113 ) histos.fill(HIST("muonMomPDG"), 3.0, 1);
              if(fabs(mcMomPDG) == 130 ) histos.fill(HIST("muonMomPDG"), 4.0, 1);
              if(fabs(mcMomPDG) == 221 ) histos.fill(HIST("muonMomPDG"), 5.0, 1);
              if(fabs(mcMomPDG) == 223 ) histos.fill(HIST("muonMomPDG"), 6.0, 1);
              if(fabs(mcMomPDG) == 310 ) histos.fill(HIST("muonMomPDG"), 7.0, 1);
              if(fabs(mcMomPDG) == 313 ) histos.fill(HIST("muonMomPDG"), 8.0, 1);
              if(fabs(mcMomPDG) == 321 ) histos.fill(HIST("muonMomPDG"), 9.0, 1);
              if(fabs(mcMomPDG) == 323 ) histos.fill(HIST("muonMomPDG"), 10.0, 1);
              if(fabs(mcMomPDG) == 333 ) histos.fill(HIST("muonMomPDG"), 12.0, 1);
              if(fabs(mcMomPDG) == 411 ) histos.fill(HIST("muonMomPDG"), 13.0, 1);
              if(fabs(mcMomPDG) == 413 ) histos.fill(HIST("muonMomPDG"), 14.0, 1);
              if(fabs(mcMomPDG) == 431 ) histos.fill(HIST("muonMomPDG"), 16.0, 1);
              if(fabs(mcMomPDG) > 499 && fabs(mcMomPDG) < 600) histos.fill(HIST("muonMomPDG"), 17.0, 1);
              if(fabs(mcMomPDG) > 999 ) histos.fill(HIST("muonMomPDG"), 18.0, 1);
              
              if(fabs(mcMomPDG)==211){
                histos.fill(HIST("muonMomPDG"), 2.0, 1);
                histos.fill(HIST("pi2mu_DCA_X"), fwdtrack.fwdDcaX());
                histos.fill(HIST("pi2mu_DCA_Y"), fwdtrack.fwdDcaY());
                histos.fill(HIST("pi2mu_DCAT"), sqrt(pow(fwdtrack.fwdDcaX(), 2.0)*pow(fwdtrack.fwdDcaY(), 2.0)));
              }
              if(fabs(mcMomPDG)==421){
                //auto momID = mcMom.globalIndex();
                histos.fill(HIST("muonMomPDG"), 15.0, 1);
                histos.fill(HIST("D0DecayX"), mcParticle_fwd.vx());
                histos.fill(HIST("D0DecayY"), mcParticle_fwd.vy());
                histos.fill(HIST("D0DecayZ"), mcParticle_fwd.vz());
                histos.fill(HIST("D0_pT"), mcMom.pt());
                histos.fill(HIST("Distance_D0_z"), fabs(mcCol_z-mcParticle_fwd.vz()));

                /*int d[10] = {};
                int gi[10] = {};
                float eta_k[10] = {};
                int i = 0;
                for(auto Daughter : Daughters){
                  d[i] = Daughter.pdgCode();
                  gi[i] = Daughter.globalIndex();
                  eta_k[i] = Daughter.eta();
                  i++;
                }

                if(d[0]*d[1]!=4173) continue;
                if(eta_k[0]>(-2.45) || eta_k[0]<(-3.6)) continue;
                histos.fill(HIST("counter1"), 0.5);*/

                for(auto Daughter : Daughters){
                  if(fabs(Daughter.pdgCode())==13){//Daughter.globalIndex()==gi[1]){ //muon
                    //if(eta_k[0]<(-2.45) && eta_k[0]>(-3.6))
                    auto mu_ID = Daughter.globalIndex();
                    if(mu_ID==a_ID){
                      mudcaX = fwdtrack.fwdDcaX();
                      mudcaY = fwdtrack.fwdDcaY();
                      mumftX = fwdtrack.x();
                      mumftY = fwdtrack.y();
                      mumftZ = fwdtrack.z();
                      histos.fill(HIST("D02mu_DCA_X"), mudcaX);
                      histos.fill(HIST("D02mu_DCA_Y"), mudcaY);
                      histos.fill(HIST("D02mu_DCAT"), sqrt(pow(mudcaX, 2.0)*pow(mudcaY, 2.0)));
                    }
                  }
                  /*if(Daughter.globalIndex()==gi[0]){ //kaon
                    if(eta_k[0]<(-2.45) && eta_k[0]>(-3.6)){
                      for(auto& mfttrack : mfttracks){
                        if(!mfttrack.has_mcParticle()) continue;
                        auto mcParticle_mft = mfttrack.mcParticle();
                        if(mcParticle_mft.globalIndex()==Daughter.globalIndex()){
                          histos.fill(HIST("counter2"), 0.5);
                          histos.fill(HIST("MFT_Kaon_pT"), mcParticle_mft.pt());
                        }
                      }
                    }
                  }*/
                }
              }
              if(fabs(mcMomPDG)==321){
                histos.fill(HIST("muonMomPDG"), 11.0, 1);
                histos.fill(HIST("KpDecayX"), mcParticle_fwd.vx());
                histos.fill(HIST("KpDecayY"), mcParticle_fwd.vy());
                histos.fill(HIST("KpDecayZ"), mcParticle_fwd.vz());
                histos.fill(HIST("Kp2mu_DCA_X"), fwdtrack.fwdDcaX());
                histos.fill(HIST("Kp2mu_DCA_Y"), fwdtrack.fwdDcaY());
                histos.fill(HIST("Kp2mu_DCAT"), sqrt(pow(fwdtrack.fwdDcaX(), 2.0)*pow(fwdtrack.fwdDcaY(), 2.0)));
              }
            }
          }
        }
      }
      //Information of MFTTracks
      for(auto& mfttrack : mfttracks){
        if(!mfttrack.has_mcParticle()) continue;
        auto mcParticle_mft = mfttrack.mcParticle();
        histos.fill(HIST("MFTParticleCounter"), 0.5);
        if(fabs(mcParticle_mft.pdgCode())==13){
          histos.fill(HIST("MFTMuonCounter"), 0.5);
          auto mft_mupT = mcParticle_mft.pt();
          histos.fill(HIST("MFT_muon_pT"), mft_mupT);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyAnalysisTask>(cfgc)};
}
