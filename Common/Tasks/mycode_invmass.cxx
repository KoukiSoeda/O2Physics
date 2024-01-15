/// \author Koki Soeda
/// \since 13/09/2023

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
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps>;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

vector<double> PCA_Cal(vector<double> &vec_mftdca)
{
  vector<double> vec_pca;
  double mumftX = vec_mftdca[0];
  double mumftY = vec_mftdca[1];
  double mumftZ = vec_mftdca[2];
  double mudcaX = vec_mftdca[3];
  double mudcaY = vec_mftdca[4];
  double canmftX = vec_mftdca[5];
  double canmftY = vec_mftdca[6];
  double canmftZ = vec_mftdca[7];
  double candcaX = vec_mftdca[8];
  double candcaY = vec_mftdca[9];
  double Col_x = vec_mftdca[10];
  double Col_y = vec_mftdca[11];
  double Col_z = vec_mftdca[12];

  auto unit_Na = sqrt(pow(mumftX-mudcaX,2)+pow(mumftY-mudcaY,2)+pow(mumftZ-Col_z,2));
  auto unit_Nc = sqrt(pow(canmftX-candcaX,2)+pow(canmftY-candcaY,2)+pow(canmftZ-Col_z,2));
  auto Nax = (mumftX-mudcaX)/unit_Na;
  auto Nay = (mumftY-mudcaY)/unit_Na;
  auto Naz = (mumftZ-Col_z)/unit_Na;
  auto Ncx = (canmftX-candcaX)/unit_Nc;
  auto Ncy = (canmftY-candcaY)/unit_Nc;
  auto Ncz = (canmftZ-Col_z)/unit_Nc;
  auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
  auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
  auto A3 = (mudcaX-candcaX)*Nax + (mudcaY-candcaY)*Nay + (Col_z-Col_z)*Naz;
  auto B1 = A2;
  auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
  auto B3 = (candcaX-mudcaX)*Ncx + (candcaY-mudcaY)*Ncy + (Col_z-Col_z)*Ncz;
  auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
  auto s = -((A2*t+A3)/A1);

  double predict_mux = mudcaX + s*Nax;
  double predict_muy = mudcaY + s*Nay;
  double predict_muz = Col_z + s*Naz;
  double predict_canx = candcaX + t*Ncx;
  double predict_cany = candcaY + t*Ncy;
  double predict_canz = Col_z + t*Ncz;
  double r_xyz = sqrt(pow(predict_canx-predict_mux,2) + pow(predict_cany-predict_muy,2) + pow(predict_canz-predict_muz,2));
  auto vecx_mu = mumftX - mudcaX;
  auto vecy_mu = mumftY - mudcaY;
  auto vecz_mu = mumftZ - Col_z;
  auto vecx_can = canmftX - candcaX;
  auto vecy_can = canmftY - candcaY;
  auto vecz_can = canmftZ - Col_z;
  auto cosxy = (vecx_mu*vecx_can + vecy_mu*vecy_can + vecz_mu*vecz_can)/((sqrt(pow(vecx_mu,2)+pow(vecy_mu,2)+pow(vecz_mu,2)))*(sqrt(pow(vecx_can,2)+pow(vecy_can,2)+pow(vecz_can,2))));

  auto pcaX = ((predict_mux+predict_canx)/2)-Col_x;
  auto pcaY = ((predict_muy+predict_cany)/2)-Col_y;
  auto pcaZ = ((predict_muz+predict_canz)/2)-Col_z;
  auto pcaD = sqrt(pow(pcaX,2)+pow(pcaY,2)+pow(pcaZ,2));
  vec_pca = {r_xyz, cosxy, pcaX, pcaY, pcaZ, pcaD};
  return vec_pca;
}

double InvMass(vector<double> &mass_para)
{
  double inv_mass=0, inv_mass2=0;
  const double mu_mass2 = 0.105658*0.105658;
  //const double pi_mass2 = 0.139571*0.139571;
  const double k_mass2 = 0.493677*0.493677;
  double muon_p2 = mass_para[0]*mass_para[0];
  double muon_px = mass_para[1];
  double muon_py = mass_para[2];
  double muon_pz = mass_para[3];
  double pair_p2 = mass_para[4]*mass_para[4];
  double pair_px = mass_para[5];
  double pair_py = mass_para[6];
  double pair_pz = mass_para[7];
  inv_mass2 = pow(sqrt(mu_mass2+muon_p2)+sqrt(k_mass2+pair_p2),2) - (pow(muon_px+pair_px,2)+pow(muon_py+pair_py,2)+pow(muon_pz+pair_pz,2));
  inv_mass = sqrt(inv_mass2);
  return inv_mass;
}

struct ProcessInvMass{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&){
    const AxisSpec axisMass{5000, 0, 50, "(GeV/c)"};
    const AxisSpec axisP{10000, 0, 100, "(GeV/c)"};
    const AxisSpec axisPT{10000, 0, 50, "(GeV/c)"};

    histos.add("HF_RecoInvMass_MC", "HF_RecoInvMass_MC", kTH1F, {axisMass});
    histos.add("LF_RecoInvMass_MC", "LF_RecoInvMass_MC", kTH1F, {axisMass});
    histos.add("HF_RecoInvMass_Det", "HF_RecoInvMass_Det", kTH1F, {axisMass});
    histos.add("LF_RecoInvMass_Det", "LF_RecoInvMass_Det", kTH1F, {axisMass});
    //MC Info.
    histos.add("HF_Muon_pT_MC", "HF_Muon_pT_MC", kTH1F, {axisPT});
    histos.add("LF_Muon_pT_MC", "LF_Muon_pT_MC", kTH1F, {axisPT});
    histos.add("HF_Muon_p_MC", "HF_Muon_p_MC", kTH1F, {axisP});
    histos.add("LF_Muon_p_MC", "LF_Muon_p_MC", kTH1F, {axisP});
    histos.add("HF_Muon_pz_MC", "HF_Muon_pz_MC", kTH1F, {axisP});
    histos.add("LF_Muon_pz_MC", "LF_Muon_pz_MC", kTH1F, {axisP});
    histos.add("HF_Kaon_pT_MC", "HF_Kaon_pT_MC", kTH1F, {axisPT});
    histos.add("MFT_Kaon_pT_MC", "MFT_Kaon_pT_MC", kTH1F, {axisPT});
    histos.add("MFT_Pion_pT_MC", "MFT_Pion_pT_MC", kTH1F, {axisPT});
    //Including Detector effects
    histos.add("HF_Muon_pT", "HF_Muon_pT", kTH1F, {axisPT});
    histos.add("LF_Muon_pT", "LF_Muon_pT", kTH1F, {axisPT});
    histos.add("HF_Muon_p", "HF_Muon_p", kTH1F, {axisP});
    histos.add("LF_Muon_p", "LF_Muon_p", kTH1F, {axisP});
    histos.add("HF_Muon_pz", "HF_Muon_pz", kTH1F, {axisP});
    histos.add("LF_Muon_pz", "LF_Muon_pz", kTH1F, {axisP});
    histos.add("HF_Kaon_pT", "HF_Kaon_pT", kTH1F, {axisPT});
    histos.add("MFT_Kaon_pT", "MFT_Kaon_pT", kTH1F, {axisPT});
    histos.add("MFT_Pion_pT", "MFT_Pion_pT", kTH1F, {axisPT});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions,
               soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA, aod::FwdTracksCov> const& fwdtracks, 
               soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
               aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
               aod::McCollisions const&)
  {
    for(auto& fwdtrack : fwdtracks){
      if(!fwdtrack.has_mcParticle()) continue;
      if(fwdtrack.trackType()==1 || fwdtrack.trackType()==3 || fwdtrack.trackType()==4) continue;
      if(fwdtrack.eta()<=(-3.6) || fwdtrack.eta()>=(-2.5)) continue;

      auto mcParticle_fwd = fwdtrack.mcParticle();
      int mcParticle_id = fwdtrack.mcParticleId();
      int mcCol_id;
      double Col_x, Col_y, Col_z;
      bool Has_MCCol = false;
      for(auto& t_col : collisions){
        if(!t_col.has_mcCollision()) continue;
        if(mcParticle_fwd.mcCollisionId()==t_col.mcCollisionId()){
          mcCol_id = t_col.mcCollisionId();
          Col_x = t_col.posX();
          Col_y = t_col.posY();
          Col_z = t_col.posZ();
          Has_MCCol = true;
          break;
        }
      }
      if(fabs(Col_z)>10 || Has_MCCol==false) continue;
      if(fabs(mcParticle_fwd.pdgCode())!=13 || !mcParticle_fwd.has_mothers()) continue;

      int momID = mcParticle_fwd.mothers_first_as<aod::McParticles>().globalIndex();
      double mc_muon_p = mcParticle_fwd.p();
      double mc_muon_px = mcParticle_fwd.px();
      double mc_muon_py = mcParticle_fwd.py();
      double mc_muon_pz = mcParticle_fwd.pz();
      double mc_muon_pT = mcParticle_fwd.pt();

      double muon_p = fwdtrack.p();
      double muon_px = fwdtrack.px();
      double muon_py = fwdtrack.py();
      double muon_pz = fwdtrack.pz();
      double muon_pT = fwdtrack.pt();

      int PartType = 0;
      bool IsSecondary = false;
      bool HasLightParent = false;
      bool HasCharmParent = false;
      bool HasBeautyParent = false;
      while(mcParticle_fwd.has_mothers()){
        mcParticle_fwd = *(mcParticle_fwd.mothers_first_as<aod::McParticles>());
        const int pdgAbs(fabs(mcParticle_fwd.pdgCode()));
        if(pdgAbs<10) break; // Quark
        if(!mcParticle_fwd.producedByGenerator()){ // Produced in transport code
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
        histos.fill(HIST("HF_Muon_pT_MC"), mc_muon_pT);
        histos.fill(HIST("HF_Muon_p_MC"), mc_muon_p);
        histos.fill(HIST("HF_Muon_pz_MC"), mc_muon_pz);
        histos.fill(HIST("HF_Muon_pT"), muon_pT);
        histos.fill(HIST("HF_Muon_p"), muon_p);
        histos.fill(HIST("HF_Muon_pz"), muon_pz);
      }else if(HasLightParent==true && IsSecondary==false){
        PartType = 4; // LF Muon(BG)
        histos.fill(HIST("LF_Muon_pT_MC"), mc_muon_pT);
        histos.fill(HIST("LF_Muon_p_MC"), mc_muon_p);
        histos.fill(HIST("LF_Muon_pz_MC"), mc_muon_pz);
        histos.fill(HIST("LF_Muon_pT"), muon_pT);
        histos.fill(HIST("LF_Muon_p"), muon_p);
        histos.fill(HIST("LF_Muon_pz"), muon_pz);
      }
      if(PartType==1 || PartType==2) continue;

      for(auto& mfttrack : mfttracks){
        if(!mfttrack.has_mcParticle()) continue;
        if(mfttrack.mcParticle().mcCollisionId()!=mcCol_id) continue;
        if(mcParticle_id==mfttrack.mcParticleId() || momID==mfttrack.mcParticleId()) continue;
        
        vector<double> forinvmass_mc;
        vector<double> forinvmass;
        forinvmass_mc.clear();
        forinvmass.clear();

        auto mcParticle_mft = mfttrack.mcParticle();
        int mfttrack_pdg = fabs(mcParticle_mft.pdgCode());
        double mc_mft_pt = mcParticle_mft.pt();
        double mc_mft_p = mcParticle_mft.p();
        double mc_mft_px = mcParticle_mft.px();
        double mc_mft_py = mcParticle_mft.py();
        double mc_mft_pz = mcParticle_mft.pz();

        double mft_pt = mfttrack.pt();
        double mft_p = mfttrack.p();
        double mft_px = mfttrack.px();
        double mft_py = mfttrack.py();
        double mft_pz = mfttrack.pz();

        int PartType_mft = 0;
        bool IsSecondary = false;
        bool HasLightParent = false;
        bool HasCharmParent = false;
        bool HasBeautyParent = false;
        while(mcParticle_mft.has_mothers()){
          mcParticle_mft = *(mcParticle_mft.mothers_first_as<aod::McParticles>());
          const int pdgAbs(fabs(mcParticle_mft.pdgCode()));
          if(pdgAbs<10) break; // Quark
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
          PartType_mft = 1; // Beauty Decay Muon
        }else if(HasLightParent==false && HasBeautyParent==true && HasCharmParent==true){
          PartType_mft = 2; // Non Prompt Charm Muon
        }else if(HasLightParent==false && HasBeautyParent==false && HasCharmParent==true){
          PartType_mft = 3; // Prompt Charm Muon
          if(mfttrack_pdg==321){
            histos.fill(HIST("HF_Kaon_pT_MC"), mc_mft_pt);
            histos.fill(HIST("HF_Kaon_pT"), mft_pt);
          }
        }else if(HasLightParent==true && IsSecondary==false){
          PartType_mft = 4; // LF Muon(BG)
        }

        if(PartType_mft!=3){
          if(mfttrack_pdg==321){
            histos.fill(HIST("MFT_Kaon_pT_MC"), mc_mft_pt);
            histos.fill(HIST("MFT_Kaon_pT"), mft_pt);
          }else if(mfttrack_pdg==211){
            histos.fill(HIST("MFT_Pion_pT_MC"), mc_mft_pt);
            histos.fill(HIST("MFT_Pion_pT"), mft_pt);
          }
        }
        forinvmass_mc = {mc_muon_p, mc_muon_px, mc_muon_py, mc_muon_pz, mc_mft_p, mc_mft_px, mc_mft_py, mc_mft_pz};
        double RecoMass_MC = InvMass(forinvmass_mc);
        forinvmass = {muon_p, muon_px, muon_py, muon_pz, mft_p, mft_px, mft_py, mft_pz};
        double RecoMass = InvMass(forinvmass);
        if(PartType==3){
          histos.fill(HIST("HF_RecoInvMass_MC"), RecoMass_MC);
          histos.fill(HIST("HF_RecoInvMass_Det"), RecoMass);
        }else if(PartType==4){
          histos.fill(HIST("LF_RecoInvMass_MC"), RecoMass_MC);
          histos.fill(HIST("LF_RecoInvMass_Det"), RecoMass);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ProcessInvMass>(cfgc)
  };
}