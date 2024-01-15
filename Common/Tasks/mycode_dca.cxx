#include <iostream>
#include <fstream>
#include <array>
#include <cmath>
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "CommonConstants/GeomConstants.h"

#include "CCDB/BasicCCDBManager.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "Math/MatrixFunctions.h"
#include "MathUtils/Utils.h"
#include "TMath.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "MFTTracking/Tracker.h"
#include "Framework/ASoAHelpers.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DataFormatsParameters/GRPMagField.h"
#include <math.h>
#include <TLorentzVector.h>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::track;
using namespace evsel;
using o2::field::MagneticField;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

struct DCAAnalysis{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&){
    const AxisSpec axisDCA{10000, -5, 5, "cm"};
    const AxisSpec axisDCAT{3000, 0, 30, "cm"};
    const AxisSpec axisP{500, 0, 50, "p(GeV/c)"};
    const AxisSpec axisD0P{1000, 0, 100, "p(GeV/c)"};
    const AxisSpec axisBeta{5000, 0, 1.1, "p/E (c)"};
    const AxisSpec axispT{500, 0, 50, "p_{T}(GeV/c)"};
    const AxisSpec axisChi2{1000, 0, 100, "chi2"};
    const AxisSpec axisPhi{640, -3.2, 3.2, "rad"};
    const AxisSpec axisEta{100, -5, 5, "Eta"};
    const AxisSpec axisSigma{2000, -10, 10, "#sigma = (DCA_{MC}-DCA_{MFT})/DCA_{MC}"};

    histos.add("AllMFTTrack_DCAT", "AllMFTTrack_DCAT", kTH1F, {axisDCAT});
    histos.add("AllMFTTrack_DCAx", "AllMFTTrack_DCAx", kTH1F, {axisDCA});
    histos.add("AllMFTTrack_DCAy", "AllMFTTrack_DCAy", kTH1F, {axisDCA});

    histos.add("D0_p", "D0_p", kTH1F, {axisD0P});
    histos.add("D0_pT", "D0_pT", kTH1F, {axispT});

    histos.add("Pion_DCAT", "Pion_DCAT", kTH1F, {axisDCAT});
    histos.add("Pion_DCAx", "Pion_DCAx", kTH1F, {axisDCA});
    histos.add("Pion_DCAy", "Pion_DCAy", kTH1F, {axisDCA});
    histos.add("Pion_HF_DCAT", "Pion_HF_DCAT", kTH1F, {axisDCAT});
    histos.add("Pion_HF_p", "Pion_HF_p", kTH1F, {axisP});
    histos.add("Pion_HF_pT", "Pion_HF_pT", kTH1F, {axispT});
    histos.add("Pion_HF_beta", "Pion_HF_beta", kTH1F, {axisBeta});
    histos.add("Pion_HF_DCAx", "Pion_HF_DCAx", kTH1F, {axisDCA});
    histos.add("Pion_HF_DCAy", "Pion_HF_DCAy", kTH1F, {axisDCA});
    histos.add("Pion_HF_DCAreso_p", "Pion_HF_DCAreso_p", kTH2F, {axisSigma, axisP});
    histos.add("Pion_HF_DCAreso_p02", "Pion_HF_DCAreso_p02", kTH1F, {axisSigma});
    histos.add("Pion_HF_DCAreso_p25", "Pion_HF_DCAreso_p25", kTH1F, {axisSigma});
    histos.add("Pion_HF_DCAreso_p510", "Pion_HF_DCAreso_p510", kTH1F, {axisSigma});
    histos.add("Pion_HF_DCAreso_p1015", "Pion_HF_DCAreso_p1015", kTH1F, {axisSigma});
    histos.add("Pion_HF_DCAreso_p1525", "Pion_HF_DCAreso_p1525", kTH1F, {axisSigma});
    histos.add("Pion_HF_DCAreso_p2550", "Pion_HF_DCAreso_2550", kTH1F, {axisSigma});
    histos.add("Pion_HF_Chi2_p02", "Pion_HF_Chi2_p01", kTH1F, {axisChi2});
    histos.add("Pion_HF_Chi2_p25", "Pion_HF_Chi2_p25", kTH1F, {axisChi2});
    histos.add("Pion_HF_Chi2_p510", "Pion_HF_Chi2_p510", kTH1F, {axisChi2});
    histos.add("Pion_HF_Chi2_p1015", "Pion_HF_Chi2_p1015", kTH1F, {axisChi2});
    histos.add("Pion_HF_Chi2_p1525", "Pion_HF_Chi2_p1525", kTH1F, {axisChi2});
    histos.add("Pion_HF_Chi2_p2550", "Pion_HF_Chi2_2550", kTH1F, {axisChi2});
    histos.add("Pion_nonHF_DCAT", "Pion_HF_DCAT", kTH1F, {axisDCAT});
    histos.add("Pion_nonHF_p", "Pion_nonHF_p", kTH1F, {axisP});
    histos.add("Pion_nonHF_DCAx", "Pion_nonHF_DCAx", kTH1F, {axisDCA});
    histos.add("Pion_nonHF_DCAy", "Pion_nonHF_DCAy", kTH1F, {axisDCA});

    histos.add("Kaon_DCAT", "Kaon_DCAT", kTH1F, {axisDCAT});
    histos.add("Kaon_DCAx", "Kaon_DCAx", kTH1F, {axisDCA});
    histos.add("Kaon_DCAy", "Kaon_DCAy", kTH1F, {axisDCA});
    histos.add("Kaon_HF_DCAT", "Kaon_HF_DCAT", kTH1F, {axisDCAT});
    histos.add("Kaon_HF_p", "Kaon_HF_p", kTH1F, {axisP});
    histos.add("Kaon_HF_pT", "Kaon_HF_pT", kTH1F, {axispT});
    histos.add("Kaon_HF_beta", "Kaon_HF_beta", kTH1F, {axisBeta});
    histos.add("Kaon_HF_DCAx", "Kaon_HF_DCAx", kTH1F, {axisDCA});
    histos.add("Kaon_HF_DCAy", "Kaon_HF_DCAy", kTH1F, {axisDCA});
    histos.add("Kaon_HF_DCAreso_p", "Kaon_HF_DCAreso_p", kTH2F, {axisSigma, axisP});
    histos.add("Kaon_HF_DCAreso_p02", "Kaon_HF_DCAreso_p02", kTH1F, {axisSigma});
    histos.add("Kaon_HF_DCAreso_p25", "Kaon_HF_DCAreso_p25", kTH1F, {axisSigma});
    histos.add("Kaon_HF_DCAreso_p510", "Kaon_HF_DCAreso_p510", kTH1F, {axisSigma});
    histos.add("Kaon_HF_DCAreso_p1015", "Kaon_HF_DCAreso_p1015", kTH1F, {axisSigma});
    histos.add("Kaon_HF_DCAreso_p1525", "Kaon_HF_DCAreso_p1525", kTH1F, {axisSigma});
    histos.add("Kaon_HF_DCAreso_p2550", "Kaon_HF_DCAreso_2550", kTH1F, {axisSigma});
    histos.add("Kaon_HF_Chi2_p02", "Kaon_HF_Chi2_p01", kTH1F, {axisChi2});
    histos.add("Kaon_HF_Chi2_p25", "Kaon_HF_Chi2_p25", kTH1F, {axisChi2});
    histos.add("Kaon_HF_Chi2_p510", "Kaon_HF_Chi2_p510", kTH1F, {axisChi2});
    histos.add("Kaon_HF_Chi2_p1015", "Kaon_HF_Chi2_p1015", kTH1F, {axisChi2});
    histos.add("Kaon_HF_Chi2_p1525", "Kaon_HF_Chi2_p1525", kTH1F, {axisChi2});
    histos.add("Kaon_HF_Chi2_p2550", "Kaon_HF_Chi2_2550", kTH1F, {axisChi2});
    histos.add("Kaon_nonHF_DCAT", "Kaon_HF_DCAT", kTH1F, {axisDCAT});
    histos.add("Kaon_nonHF_p", "Kaon_nonHF_p", kTH1F, {axisP});
    histos.add("Kaon_nonHF_DCAx", "Kaon_nonHF_DCAx", kTH1F, {axisDCA});
    histos.add("Kaon_nonHF_DCAy", "Kaon_nonHF_DCAy", kTH1F, {axisDCA});

    histos.add("Muon_DCAT", "Muon_DCAT", kTH1F, {axisDCAT});
    histos.add("Muon_DCAx", "Muon_DCAx", kTH1F, {axisDCA});
    histos.add("Muon_DCAy", "Muon_DCAy", kTH1F, {axisDCA});
    histos.add("Muon_HF_DCAT", "Muon_HF_DCAT", kTH1F, {axisDCAT});
    histos.add("Muon_HF_p", "Muon_HF_p", kTH1F, {axisP});
    histos.add("Muon_HF_pT", "Muon_HF_pT", kTH1F, {axispT});
    histos.add("Muon_HF_beta", "Muon_HF_beta", kTH1F, {axisBeta});
    histos.add("Muon_HF_DCAx", "Muon_HF_DCAx", kTH1F, {axisDCA});
    histos.add("Muon_HF_DCAy", "Muon_HF_DCAy", kTH1F, {axisDCA});
    histos.add("Muon_HF_DCAreso_p", "Muon_HF_DCAreso_p", kTH2F, {axisSigma, axisP});
    histos.add("Muon_HF_DCAreso_p02", "Muon_HF_DCAreso_p02", kTH1F, {axisSigma});
    histos.add("Muon_HF_DCAreso_p25", "Muon_HF_DCAreso_p25", kTH1F, {axisSigma});
    histos.add("Muon_HF_DCAreso_p510", "Muon_HF_DCAreso_p510", kTH1F, {axisSigma});
    histos.add("Muon_HF_DCAreso_p1015", "Muon_HF_DCAreso_p1015", kTH1F, {axisSigma});
    histos.add("Muon_HF_DCAreso_p1525", "Muon_HF_DCAreso_p1525", kTH1F, {axisSigma});
    histos.add("Muon_HF_DCAreso_p2550", "Muon_HF_DCAreso_2550", kTH1F, {axisSigma});
    histos.add("Muon_HF_Chi2_p02", "Muon_HF_Chi2_p01", kTH1F, {axisChi2});
    histos.add("Muon_HF_Chi2_p25", "Muon_HF_Chi2_p25", kTH1F, {axisChi2});
    histos.add("Muon_HF_Chi2_p510", "Muon_HF_Chi2_p510", kTH1F, {axisChi2});
    histos.add("Muon_HF_Chi2_p1015", "Muon_HF_Chi2_p1015", kTH1F, {axisChi2});
    histos.add("Muon_HF_Chi2_p1525", "Muon_HF_Chi2_p1525", kTH1F, {axisChi2});
    histos.add("Muon_HF_Chi2_p2550", "Muon_HF_Chi2_2550", kTH1F, {axisChi2});
    histos.add("Muon_nonHF_DCAT", "Muon_HF_DCAT", kTH1F, {axisDCAT});
    histos.add("Muon_nonHF_p", "Muon_nonHF_p", kTH1F, {axisP});
    histos.add("Muon_nonHF_DCAx", "Muon_nonHF_DCAx", kTH1F, {axisDCA});
    histos.add("Muon_nonHF_DCAy", "Muon_nonHF_DCAy", kTH1F, {axisDCA});

    histos.add("Muon_beta099_DCAx", "Muon_beta099_DCAx", kTH1F, {axisDCA});
    histos.add("Pion_beta099_DCAx", "Pion_beta099_DCAx", kTH1F, {axisDCA});
    histos.add("Kaon_beta099_DCAx", "Kaon_beta099_DCAx", kTH1F, {axisDCA});
    histos.add("Muon_beta099_DCAy", "Muon_beta099_DCAy", kTH1F, {axisDCA});
    histos.add("Pion_beta099_DCAy", "Pion_beta099_DCAy", kTH1F, {axisDCA});
    histos.add("Kaon_beta099_DCAy", "Kaon_beta099_DCAy", kTH1F, {axisDCA});

    histos.add("GlobalFwdTrack_DCAx", "GlobalFwdTrack_DCAx", kTH1F, {axisDCA});
    histos.add("GlobalFwdTrack_DCAy", "GlobalFwdTrack_DCAy", kTH1F, {axisDCA});

    histos.add("Diff_Phi_Eta", "Diff_Phi_Eta", kTH2F, {axisPhi, axisEta});
    histos.add("Diff_PropagateMFT_MFTParameter_perP", "Diff_PropagateMFT_MFTParameter_perP", kTH2F, {axisP, axisDCAT});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions,
               soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA, aod::FwdTracksCov> const& fwdtracks, 
               soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
               aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
               aod::McCollisions const&)
  {
    //const double k_mass = 0.493677;
    //const double mu_mass = 0.105658;
    //const double pi_mass = 0.139571;

    for(auto& fwdtrack : fwdtracks){
      if(fwdtrack.trackType()!=0) continue;
      histos.fill(HIST("GlobalFwdTrack_DCAx"), fwdtrack.fwdDcaX());
      histos.fill(HIST("GlobalFwdTrack_DCAy"), fwdtrack.fwdDcaY());
    }

    // Create AmbTrackTable
    vector<uint64_t> ambTrackIds;
    ambTrackIds.clear();
    for(auto& amfttrack : amfttracks){
      ambTrackIds.push_back(amfttrack.mfttrackId());
    }

    for(auto& mfttrack : mfttracks){
      if(!mfttrack.has_mcParticle()) continue;
      auto mcParticle_mft = mfttrack.mcParticle();
      int mcCol_id;
      double Col_x, Col_y, Col_z, mcCol_x, mcCol_y, mcCol_z;
      bool MC_Col = false;
      for(auto& t_col : collisions){
        if(!t_col.has_mcCollision()) continue;
        if(mcParticle_mft.mcCollisionId()==t_col.mcCollisionId()){
          mcCol_id = t_col.mcCollisionId();
          Col_x = t_col.posX();
          Col_y = t_col.posY();
          Col_z = t_col.posZ();
          mcCol_x = t_col.mcCollision().posX();
          mcCol_y = t_col.mcCollision().posY();
          mcCol_z = t_col.mcCollision().posZ();
          MC_Col = true;
          break;
        }
      }
      if(fabs(Col_z)>10) continue;
      if(MC_Col==false) continue;
      if(mcParticle_mft.eta()>=(-2.5) || mcParticle_mft.eta()<=(-3.6)) continue;
      double decay_pos = 0;
      if(mcParticle_mft.has_daughters()){
        auto Daughters_mft = mcParticle_mft.daughters_as<aod::McParticles>();
        for(auto daughter_mft : Daughters_mft){
          decay_pos = daughter_mft.vz();
        }
      }
      if(decay_pos<(-46) && decay_pos>(-77)) continue;

      double mftchi2 = mfttrack.chi2();
      SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
      vector<double> mftv2;
      SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
      o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs2, mftchi2};
      mftpars1.propagateToZlinear(mcCol_z);
      double DCAT = sqrt(pow(mftpars1.getX()-Col_x,2)+pow(mftpars1.getY()-Col_y,2));
      double DCAX = mftpars1.getX()-mcCol_x;
      double DCAY = mftpars1.getY()-mcCol_y;
      
      double mc_sign = 0;
      if(mcParticle_mft.pdgCode()>0) mc_sign = 1;
      if(mcParticle_mft.pdgCode()<0) mc_sign = -1;
      double convert_phi = 0;
      if(mcParticle_mft.phi()>o2::math_utils::pi() && mcParticle_mft.phi()<2*o2::math_utils::pi()){
        convert_phi = mcParticle_mft.phi() - 2*(o2::math_utils::pi());
      }else{
        convert_phi = mcParticle_mft.phi();
      }
      SMatrix5 mcmftpars(mcParticle_mft.vx(), mcParticle_mft.vy(), convert_phi, tan((o2::math_utils::pi()/2)-2*atan(exp(-mcParticle_mft.eta()))), mc_sign/mcParticle_mft.pt());
      o2::track::TrackParCovFwd mcmftpars1{mcParticle_mft.vz(), mcmftpars, mftcovs2, mftchi2};
      mcmftpars1.propagateToZlinear(mcCol_z);
      double mcDCAX = mcmftpars1.getX() - mcCol_x;
      double mcDCAY = mcmftpars1.getY() - mcCol_y;
      //double mcDCAT = sqrt(pow(mcDCAX,2)+pow(mcDCAY,2));

      auto mft_pdg = fabs(mcParticle_mft.pdgCode());
      double mc_p = mcParticle_mft.p();
      double mc_e = mcParticle_mft.e();
      double mc_pt = mcParticle_mft.pt();

      mcmftpars1.propagateToZlinear(mfttrack.z());

      histos.fill(HIST("Diff_Phi_Eta"), convert_phi-mfttrack.phi(), mcParticle_mft.eta()-mfttrack.eta());
      histos.fill(HIST("Diff_PropagateMFT_MFTParameter_perP"), mc_p, sqrt(pow(mcmftpars1.getX()-mfttrack.x(),2)+pow(mcmftpars1.getY()-mfttrack.y(),2)));

      int mcMomPDG = 0;
      bool D0Daughter=false;
      if(mcParticle_mft.has_mothers()){
        auto mcMom = mcParticle_mft.mothers_first_as<aod::McParticles>();
        mcMomPDG = fabs(mcMom.pdgCode());
        if(mcMomPDG==421){
          D0Daughter=true;
          histos.fill(HIST("D0_p"), mcMom.p());
          histos.fill(HIST("D0_pT"), mcMom.pt());
        }
      }

      int PartType = 0;
      bool IsSecondary = false;
      bool HasLightParent = false;
      bool HasCharmParent = false;
      bool HasBeautyParent = false;
      int chain_number = 0;
      while(mcParticle_mft.has_mothers()){
        mcParticle_mft = *(mcParticle_mft.mothers_first_as<aod::McParticles>());
        const int pdgAbs(fabs(mcParticle_mft.pdgCode()));
        if(pdgAbs<10) break; // Quark
        chain_number++;
        if(!mcParticle_mft.producedByGenerator()){ // Produced in transport code
          IsSecondary = true;
          continue;
        }
        const int pdgRem(pdgAbs%100000);
        if(pdgRem==2212){ // Proton(Beam Particle)
          continue;
        }
        if((pdgRem<100) || (pdgRem>=10000)){ // Elementary Particles or etc.
          continue;
        }
        // Compute the flavor of constituent quark
        const int flv(pdgRem/pow(10, static_cast<int>(TMath::Log10(pdgRem))));
        if(flv>6){ // No more than 6 flavors
          LOGF(debug, "Found the bug in computing the flavor of constituent quark");
          continue;
        }
        if(flv<4){ // Light flavor
          HasLightParent = true;
          continue;
        }
        if(flv==4){
          HasCharmParent = true;
          continue;
        }
        if(flv==5){
          HasBeautyParent = true;
          continue;
        }
      }
      if(HasLightParent==false && HasBeautyParent==true && HasCharmParent==false){
        PartType = 1; // Beauty Decay Muon
      }else if(HasLightParent==false && HasBeautyParent==true && HasCharmParent==true){
        PartType = 2; // Non Prompt Charm Muon
      }else if(HasLightParent==false && HasBeautyParent==false && HasCharmParent==true){
        PartType = 3; // Prompt Charm Muon
      }else if(HasLightParent==true && IsSecondary==false){
        PartType = 4; // LF Muon(BG)
      }
      
      if(mft_pdg==13 || mft_pdg==211 || mft_pdg==321){
        histos.fill(HIST("AllMFTTrack_DCAT"), DCAT);
        histos.fill(HIST("AllMFTTrack_DCAx"), DCAX);
        histos.fill(HIST("AllMFTTrack_DCAy"), DCAY);
        if(mft_pdg==13){
          histos.fill(HIST("Muon_DCAT"), DCAT);
          histos.fill(HIST("Muon_DCAx"), DCAX);
          histos.fill(HIST("Muon_DCAy"), DCAY);
          if(mcMomPDG==421 && PartType==3 && chain_number==1){
            histos.fill(HIST("Muon_HF_DCAT"), DCAT);
            histos.fill(HIST("Muon_HF_p"), mc_p);
            histos.fill(HIST("Muon_HF_pT"), mc_pt);
            histos.fill(HIST("Muon_HF_DCAx"), mcDCAX);
            histos.fill(HIST("Muon_HF_DCAy"), mcDCAY);
            histos.fill(HIST("Muon_HF_DCAreso_p"), ((mcDCAX-DCAX)/mcDCAX), mc_p);
            if(mc_p<2){
              histos.fill(HIST("Muon_HF_DCAreso_p02"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Muon_HF_Chi2_p02"), mftchi2);
            }else if(mc_p>=2 && mc_p<5){
              histos.fill(HIST("Muon_HF_DCAreso_p25"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Muon_HF_Chi2_p25"), mftchi2);
            }else if(mc_p>=5 && mc_p<10){
              histos.fill(HIST("Muon_HF_DCAreso_p510"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Muon_HF_Chi2_p510"), mftchi2);
            }else if(mc_p>=10 && mc_p<15){
              histos.fill(HIST("Muon_HF_DCAreso_p1015"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Muon_HF_Chi2_p1015"), mftchi2);
            }else if(mc_p>=15 && mc_p<25){
              histos.fill(HIST("Muon_HF_DCAreso_p1525"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Muon_HF_Chi2_p1525"), mftchi2);
            }else if(mc_p>=25 && mc_p<50){
              histos.fill(HIST("Muon_HF_DCAreso_p2550"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Muon_HF_Chi2_p2550"), mftchi2);
            }
            histos.fill(HIST("Muon_HF_beta"), mc_p/mc_e);

            if((mc_p/mc_e)>0.99){
              histos.fill(HIST("Muon_beta099_DCAx"), DCAX);
              histos.fill(HIST("Muon_beta099_DCAy"), DCAY);
            }
            
          }else if(PartType!=3){
            histos.fill(HIST("Muon_nonHF_DCAT"), DCAT);
            histos.fill(HIST("Muon_nonHF_p"), mc_p);         
            histos.fill(HIST("Muon_nonHF_DCAx"), DCAX);
            histos.fill(HIST("Muon_nonHF_DCAy"), DCAY);
          }
        }else if(mft_pdg==211){
          histos.fill(HIST("Pion_DCAT"), DCAT);
          histos.fill(HIST("Pion_DCAx"), DCAX);
          histos.fill(HIST("Pion_DCAy"), DCAY);
          if(mcMomPDG==421 && PartType==3 && chain_number==1){
            histos.fill(HIST("Pion_HF_DCAT"), DCAT);
            histos.fill(HIST("Pion_HF_p"), mc_p);
            histos.fill(HIST("Pion_HF_pT"), mc_pt);
            histos.fill(HIST("Pion_HF_DCAx"), mcDCAX);
            histos.fill(HIST("Pion_HF_DCAy"), mcDCAY);
            histos.fill(HIST("Pion_HF_DCAreso_p"), ((mcDCAX-DCAX)/mcDCAX), mc_p);
            if(mc_p<2){
              histos.fill(HIST("Pion_HF_DCAreso_p02"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Pion_HF_Chi2_p02"), mftchi2);
            }else if(mc_p>=2 && mc_p<5){
              histos.fill(HIST("Pion_HF_DCAreso_p25"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Pion_HF_Chi2_p25"), mftchi2);
            }else if(mc_p>=5 && mc_p<10){
              histos.fill(HIST("Pion_HF_DCAreso_p510"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Pion_HF_Chi2_p510"), mftchi2);
            }else if(mc_p>=10 && mc_p<15){
              histos.fill(HIST("Pion_HF_DCAreso_p1015"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Pion_HF_Chi2_p1015"), mftchi2);
            }else if(mc_p>=15 && mc_p<25){
              histos.fill(HIST("Pion_HF_DCAreso_p1525"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Pion_HF_Chi2_p1525"), mftchi2);
            }else if(mc_p>=25 && mc_p<50){
              histos.fill(HIST("Pion_HF_DCAreso_p2550"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Pion_HF_Chi2_p2550"), mftchi2);
            }
            histos.fill(HIST("Pion_HF_beta"), mc_p/mc_e);
            if((mc_p/mc_e)>0.99){
              histos.fill(HIST("Pion_beta099_DCAx"), DCAX);
              histos.fill(HIST("Pion_beta099_DCAy"), DCAY);
            }
            
          }else if(PartType!=3){
            histos.fill(HIST("Pion_nonHF_DCAT"), DCAT);
            histos.fill(HIST("Pion_nonHF_p"), mc_p);            
            histos.fill(HIST("Pion_nonHF_DCAx"), DCAX);
            histos.fill(HIST("Pion_nonHF_DCAy"), DCAY);
          }
        }else if(mft_pdg==321){
          histos.fill(HIST("Kaon_DCAT"), DCAT);
          histos.fill(HIST("Kaon_DCAx"), DCAX);
          histos.fill(HIST("Kaon_DCAy"), DCAY);
          if(mcMomPDG==421 && PartType==3 && chain_number==1){
            histos.fill(HIST("Kaon_HF_DCAT"), DCAT);
            histos.fill(HIST("Kaon_HF_p"), mc_p);
            histos.fill(HIST("Kaon_HF_pT"), mc_pt);    
            histos.fill(HIST("Kaon_HF_DCAx"), mcDCAX);
            histos.fill(HIST("Kaon_HF_DCAy"), mcDCAY);
            histos.fill(HIST("Kaon_HF_DCAreso_p"), ((mcDCAX-DCAX)/mcDCAX), mc_p);
            if(mc_p<2){
              histos.fill(HIST("Kaon_HF_DCAreso_p02"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Kaon_HF_Chi2_p02"), mftchi2);
            }else if(mc_p>=2 && mc_p<5){
              histos.fill(HIST("Kaon_HF_DCAreso_p25"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Kaon_HF_Chi2_p25"), mftchi2);
            }else if(mc_p>=5 && mc_p<10){
              histos.fill(HIST("Kaon_HF_DCAreso_p510"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Kaon_HF_Chi2_p510"), mftchi2);
            }else if(mc_p>=10 && mc_p<15){
              histos.fill(HIST("Kaon_HF_DCAreso_p1015"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Kaon_HF_Chi2_p1015"), mftchi2);
            }else if(mc_p>=15 && mc_p<25){
              histos.fill(HIST("Kaon_HF_DCAreso_p1525"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Kaon_HF_Chi2_p1525"), mftchi2);
            }else if(mc_p>=25 && mc_p<50){
              histos.fill(HIST("Kaon_HF_DCAreso_p2550"), (mcDCAX-DCAX)/mcDCAX);
              histos.fill(HIST("Kaon_HF_Chi2_p2550"), mftchi2);
            }
            histos.fill(HIST("Kaon_HF_beta"), mc_p/mc_e);
            if((mc_p/mc_e)>0.99){
              histos.fill(HIST("Kaon_beta099_DCAx"), DCAX);
              histos.fill(HIST("Kaon_beta099_DCAy"), DCAY);
            }
            
          }else if(PartType!=3){
            histos.fill(HIST("Kaon_nonHF_DCAT"), DCAT);
            histos.fill(HIST("Kaon_nonHF_p"), mc_p);            
            histos.fill(HIST("Kaon_nonHF_DCAx"), DCAX);
            histos.fill(HIST("Kaon_nonHF_DCAy"), DCAY);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DCAAnalysis>(cfgc)
  };
}