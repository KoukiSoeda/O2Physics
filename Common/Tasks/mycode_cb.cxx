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
  double mcColID = vec_mftdca[13];
  double add_dist = vec_mftdca[14];
  double ans = vec_mftdca[15];
  double mudcat = vec_mftdca[16];
  double mupt = vec_mftdca[17];
  double mftchi2 = vec_mftdca[18];

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
  vec_pca = {r_xyz, cosxy, pcaX, pcaY, pcaZ, pcaD, add_dist, mcColID, mudcat, mupt, mftchi2, ans};
  return vec_pca;
}

struct MyAnalysisTask{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<double> pT_cut{"pT_cut", 0, "Cut value for transverse momentum"};
  Configurable<double> MFT_chi2{"MFTTrack_chi2", 2.5, "MFT tracks chi2"};
  Configurable<double> PCA_cut{"PCA_cut", 0.1, "Cut value for PCA"};
  Configurable<double> CosSim_cut{"CosSim_cut", 0.99, "Cut value for Cosine Simiraliyty"};
  Configurable<double> Total_dist_cut{"Total_dist_cut", 30.0, "Cut value for total distance with mu track"};
  Configurable<double> DCAT_min_c{"DCAT_min_c", 0.01, "Lower limit of DCAT as D"};
  Configurable<double> DCAT_max_c{"DCAT_max_c", 0.1, "Upper limit of DCAT as D"};
  Configurable<double> DCAT_min_b{"DCAT_min_b", 0.1, "Lower limit of DCAT as B"};
  Configurable<double> DCAT_max_b{"DCAT_max_b", 0.2, "Upper limit of DCAT as B"};
  Configurable<double> PCAD_min_c{"PCAD_min_c", 0.15, "Lower limit of PCAD as D"};
  Configurable<double> PCAD_max_c{"PCAD_max_c", 1, "Upper limit of PCAD as D"};
  Configurable<double> PCAD_min_b{"PCAD_min_b", 1, "Lower limit of PCAD as B"};
  Configurable<double> PCAD_max_b{"PCAD_max_b", 5, "Upper limit of PCAD as B"};
  int aod_counts = 1;

  void init(InitContext const&){
    const AxisSpec axispT{100, 0, 50, "p_{T} (GeV/c)"};

    histos.add("AllMuon_pT", "AllMuon_pT", kTH1F, {axispT});
    histos.add("AllBeauty_pT", "AllBeauty_pT", kTH1F, {axispT});
    histos.add("AllpCharm_pT", "AllpCharm_pT", kTH1F, {axispT});
    histos.add("AllnpCharm_pT", "AllnpCharm_pT", kTH1F, {axispT});
    histos.add("AllBG_pT", "AllBG_pT", kTH1F, {axispT});

    // Using cut value for D muon
    histos.add("AsD_DCATCut_Muon_pT", "AsD_DCATCut_Muon_pT", kTH1F, {axispT});
    histos.add("AsD_DCATCut_Beauty_pT", "AsD_DCATCut_Beauty_pT", kTH1F, {axispT});
    histos.add("AsD_DCATCut_pCharm_pT", "AsD_DCATCut_pCharm_pT", kTH1F, {axispT});
    histos.add("AsD_DCATCut_npCharm_pT", "AsD_DCATCut_npCharm_pT", kTH1F, {axispT});
    histos.add("AsD_DCATCut_BG_pT", "AsD_DCATCut_BG_pT", kTH1F, {axispT});
    histos.add("AsD_DCAT_PCACut_Muon_pT", "AsD_DCAT_PCACut_Muon_pT", kTH1F, {axispT});
    histos.add("AsD_DCAT_PCACut_Beauty_pT", "AsD_DCAT_PCACut_Beauty_pT", kTH1F, {axispT});
    histos.add("AsD_DCAT_PCACut_pCharm_pT", "AsD_DCAT_PCACut_pCharm_pT", kTH1F, {axispT});
    histos.add("AsD_DCAT_PCACut_npCharm_pT", "AsD_DCAT_PCACut_npCharm_pT", kTH1F, {axispT});
    histos.add("AsD_DCAT_PCACut_BG_pT", "AsD_DCAT_PCACut_BG_pT", kTH1F, {axispT});

    // Using cut value for B and nonPrompt D muon
    histos.add("AsB_DCATCut_Muon_pT", "AsB_DCATCut_Muon_pT", kTH1F, {axispT});
    histos.add("AsB_DCATCut_Beauty_pT", "AsB_DCATCut_Beauty_pT", kTH1F, {axispT});
    histos.add("AsB_DCATCut_pCharm_pT", "AsB_DCATCut_pCharm_pT", kTH1F, {axispT});
    histos.add("AsB_DCATCut_npCharm_pT", "AsB_DCATCut_npCharm_pT", kTH1F, {axispT});
    histos.add("AsB_DCATCut_BG_pT", "AsB_DCATCut_BG_pT", kTH1F, {axispT});
    histos.add("AsB_DCAT_PCACut_Muon_pT", "AsB_DCAT_PCACut_Muon_pT", kTH1F, {axispT});
    histos.add("AsB_DCAT_PCACut_Beauty_pT", "AsB_DCAT_PCACut_Beauty_pT", kTH1F, {axispT});
    histos.add("AsB_DCAT_PCACut_pCharm_pT", "AsB_DCAT_PCACut_pCharm_pT", kTH1F, {axispT});
    histos.add("AsB_DCAT_PCACut_npCharm_pT", "AsB_DCAT_PCACut_npCharm_pT", kTH1F, {axispT});
    histos.add("AsB_DCAT_PCACut_BG_pT", "AsB_DCAT_PCACut_BG_pT", kTH1F, {axispT});

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
      if(fwdtrack.pt()<pT_cut) continue;

      auto mcParticle_fwd = fwdtrack.mcParticle();
      double mcCol_id;
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

      int muid = mcParticle_fwd.globalIndex();
      auto mumom = mcParticle_fwd.mothers_first_as<aod::McParticles>();
      auto Daughters = mumom.daughters_as<aod::McParticles>();

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
      }else if(HasLightParent==true && IsSecondary==false){
        PartType = 4; // LF Muon(BG)
      }

      vector<int> vec_daughter;
      vec_daughter.clear();
      bool HasPair = false;
      for(auto& Daughter : Daughters){
        vec_daughter.push_back(Daughter.globalIndex());
      }
      if(vec_daughter.size()>1) HasPair = true;

      bool HasMuonTrack=false, HasPair1Track=false, HasPair2Track=false, HasPair3Track=false;
      double mudcaX, mudcaY, mumftX, mumftY, mumftZ;
      double p1dcaX, p1dcaY, p1mftX, p1mftY, p1mftZ, p1chi2;
      double p2dcaX, p2dcaY, p2mftX, p2mftY, p2mftZ, p2chi2;
      double p3dcaX, p3dcaY, p3mftX, p3mftY, p3mftZ, p3chi2;
      double bgdcaX, bgdcaY, bgmftX, bgmftY, bgmftZ, bgchi2;
      int mcColId_mu;

      vector<double> muon_parameter;
      muon_parameter.clear();
      for(auto& mutrack : mfttracks){
        if(!mutrack.has_mcParticle()) continue;
        if(mutrack.mcParticle().mcCollisionId()==mcCol_id){
          auto mcParticle_mft = mutrack.mcParticle();
          if(mcParticle_mft.globalIndex()==muid){
            HasMuonTrack = true;
            double muchi2 = mutrack.chi2();
            vector<double> mucov;
            SMatrix55 mucovs2(mucov.begin(), mucov.end());
            SMatrix5 mupars(mutrack.x(), mutrack.y(), mutrack.phi(), mutrack.tgl(), mutrack.signed1Pt());
            o2::track::TrackParCovFwd mupars1{mutrack.z(), mupars, mucovs2, muchi2};
            mupars1.propagateToZlinear(Col_z);

            mudcaX = mupars1.getX();
            mudcaY = mupars1.getY();
            mumftX = mutrack.x();
            mumftY = mutrack.y();
            mumftZ = mutrack.z();
            mcColId_mu = mcParticle_mft.mcCollisionId();
            muon_parameter = {mutrack.x(), mutrack.y(), mutrack.z(), mutrack.phi(), mutrack.tgl(), mutrack.signed1Pt(), mutrack.chi2()};
          }
        }
      }
      if(HasMuonTrack==false && PartType==4){
        HasMuonTrack = true;
        double mftchi2 = fwdtrack.chi2();
        vector<double> mucov;
        SMatrix55 mucovs2(mucov.begin(), mucov.end());
        SMatrix5 mftpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
        o2::track::TrackParCovFwd mftpars1{fwdtrack.z(), mftpars, mucovs2, mftchi2};
        mftpars1.propagateToZlinear(Col_z);
        mudcaX = mftpars1.getX();
        mudcaY = mftpars1.getY();
        mumftX = fwdtrack.x();
        mumftY = fwdtrack.y();
        mumftZ = fwdtrack.z();
        muon_parameter = {fwdtrack.x(), fwdtrack.y(), fwdtrack.z(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt(), fwdtrack.chi2()};
      }

      if(HasMuonTrack==false) continue;

      vector<double> mftv1;
      SMatrix55 mftcovs1(mftv1.begin(), mftv1.end());
      SMatrix5 muonpars(muon_parameter[0], muon_parameter[1], muon_parameter[3], muon_parameter[4], muon_parameter[5]);
      o2::track::TrackParCovFwd mftpars1{muon_parameter[2], muonpars, mftcovs1, muon_parameter[6]};

      /////////////////////////////////////////////////////////////////
      //  12.14にやること　→ (完了)DCAカットだけのコードの追加。(321行目あたりに)    //
      //                (完了)PCAカットのコード追加。                       //
      //                解析コードを書く。(for_analysis)                  //
      ////////////////////////////////////////////////////////////////
      double mu_pt = fwdtrack.pt();
      double muon_dcat = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));

      vector<vector<double>> vec_pca;
      vec_pca.clear();
      for(auto& mfttrack : mfttracks){
        if(!mfttrack.has_mcParticle()) continue;
        if(mfttrack.mcParticle().mcCollisionId()==mcCol_id && mfttrack.mcParticle().mcCollisionId()==mcColId_mu){
          auto mcParticle_mft = mfttrack.mcParticle();
          if(mcParticle_mft.globalIndex()==muid) continue;
          vector<double> vec_mftdca;
          vec_mftdca.clear();

          if(find(vec_daughter.begin(), vec_daughter.end(), mcParticle_mft.globalIndex())!=vec_daughter.end() && HasPair1Track==false){
            HasPair1Track = true;
            double ans1 = 1;
            double mftchi2 = mfttrack.chi2();
            SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
            vector<double> mftv2;
            SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
            o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars, mftcovs2, mftchi2};
            mftpars2.propagateToZlinear(Col_z);
            p1dcaX = mftpars2.getX();
            p1dcaY = mftpars2.getY();
            p1mftX = mfttrack.x();
            p1mftY = mfttrack.y();
            p1mftZ = mfttrack.z();
            p1chi2 = mfttrack.chi2();
            double add_dist = 0;
            for(int i=0; i<5000; i++){
              auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
              if((cut_z-Col_z)>0) break;
              mftpars1.propagateToZlinear(cut_z);
              mftpars2.propagateToZlinear(cut_z);
              auto diff = sqrt(pow(mftpars2.getX()-mftpars1.getX(),2)+pow(mftpars2.getY()-mftpars1.getY(),2));
              add_dist += diff;
            }
            vec_mftdca = {mumftX, mumftY, mumftZ, mudcaX, mudcaY, p1mftX, p1mftY, p1mftZ, p1dcaX, p1dcaY, Col_x, Col_y, Col_z, mcCol_id, add_dist, ans1, muon_dcat, mu_pt, p1chi2};
            vector<double> p1_pca = PCA_Cal(vec_mftdca);
            p1_pca.push_back(PartType);
            p1_pca.push_back(muid);
            vec_pca.push_back(p1_pca);

          }else if(find(vec_daughter.begin(), vec_daughter.end(), mcParticle_mft.globalIndex())!=vec_daughter.end() && HasPair1Track==true && HasPair2Track==false){
            HasPair2Track = true;
            double ans2 = 2;
            double mftchi2 = mfttrack.chi2();
            SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
            vector<double> mftv2;
            SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
            o2::track::TrackParCovFwd mftpars3{mfttrack.z(), mftpars, mftcovs2, mftchi2};
            mftpars3.propagateToZlinear(Col_z);
            p2dcaX = mftpars3.getX();
            p2dcaY = mftpars3.getY();
            p2mftX = mfttrack.x();
            p2mftY = mfttrack.y();
            p2mftZ = mfttrack.z();
            p2chi2 = mfttrack.chi2();
            double add_dist = 0;
            for(int i=0; i<5000; i++){
              auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
              if((cut_z-Col_z)>0) break;
              mftpars1.propagateToZlinear(cut_z);
              mftpars3.propagateToZlinear(cut_z);
              auto diff = sqrt(pow(mftpars3.getX()-mftpars1.getX(),2)+pow(mftpars3.getY()-mftpars1.getY(),2));
              add_dist += diff;
            }
            vec_mftdca = {mumftX, mumftY, mumftZ, mudcaX, mudcaY, p2mftX, p2mftY, p2mftZ, p2dcaX, p2dcaY, Col_x, Col_y, Col_z, mcCol_id, add_dist, ans2, muon_dcat, mu_pt, p2chi2};
            vector<double> p2_pca = PCA_Cal(vec_mftdca);
            p2_pca.push_back(PartType);
            p2_pca.push_back(muid);
            vec_pca.push_back(p2_pca);

          }else if(find(vec_daughter.begin(), vec_daughter.end(), mcParticle_mft.globalIndex())!=vec_daughter.end() && HasPair1Track==true && HasPair2Track==true && HasPair3Track==false){
            HasPair3Track = true;
            double ans3 = 3;
            double mftchi2 = mfttrack.chi2();
            SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
            vector<double> mftv2;
            SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
            o2::track::TrackParCovFwd mftpars4{mfttrack.z(), mftpars, mftcovs2, mftchi2};
            mftpars4.propagateToZlinear(Col_z);
            p3dcaX = mftpars4.getX();
            p3dcaY = mftpars4.getY();
            p3mftX = mfttrack.x();
            p3mftY = mfttrack.y();
            p3mftZ = mfttrack.z();
            p3chi2 = mfttrack.chi2();
            double add_dist = 0;
            for(int i=0; i<5000; i++){
              auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
              if((cut_z-Col_z)>0) break;
              mftpars1.propagateToZlinear(cut_z);
              mftpars4.propagateToZlinear(cut_z);
              auto diff = sqrt(pow(mftpars4.getX()-mftpars1.getX(),2)+pow(mftpars4.getY()-mftpars1.getY(),2));
              add_dist += diff;
            }
            vec_mftdca = {mumftX, mumftY, mumftZ, mudcaX, mudcaY, p3mftX, p3mftY, p3mftZ, p3dcaX, p3dcaY, Col_x, Col_y, Col_z, mcCol_id, add_dist, ans3, muon_dcat, mu_pt, p3chi2};
            vector<double> p3_pca = PCA_Cal(vec_mftdca);
            p3_pca.push_back(PartType);
            p3_pca.push_back(muid);
            vec_pca.push_back(p3_pca);

          }else if(find(vec_daughter.begin(), vec_daughter.end(), mcParticle_mft.globalIndex())==vec_daughter.end() || HasPair==false){
            double mftchi2 = mfttrack.chi2();
            double ans0 = 0;
            SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
            vector<double> mftv2;
            SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
            o2::track::TrackParCovFwd mftpars0{mfttrack.z(), mftpars, mftcovs2, mftchi2};
            mftpars0.propagateToZlinear(Col_z);
            bgdcaX = mftpars0.getX();
            bgdcaY = mftpars0.getY();
            bgmftX = mfttrack.x();
            bgmftY = mfttrack.y();
            bgmftZ = mfttrack.z();
            bgchi2 = mfttrack.chi2();
            double add_dist = 0;
            for(int i=0; i<5000; i++){
              auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
              if((cut_z-Col_z)>0) break;
              mftpars1.propagateToZlinear(cut_z);
              mftpars0.propagateToZlinear(cut_z);
              auto diff = sqrt(pow(mftpars0.getX()-mftpars1.getX(),2)+pow(mftpars0.getY()-mftpars1.getY(),2));
              add_dist += diff;
            }
            vec_mftdca = {mumftX, mumftY, mumftZ, mudcaX, mudcaY, bgmftX, bgmftY, bgmftZ, bgdcaX, bgdcaY, Col_x, Col_y, Col_z, mcCol_id, add_dist, ans0, muon_dcat, mu_pt, bgchi2};
            vector<double> bg_pca = PCA_Cal(vec_mftdca);
            bg_pca.push_back(PartType);
            bg_pca.push_back(muid);
            vec_pca.push_back(bg_pca);
          }
        }
      }
      if(vec_pca.empty()) continue;
      sort(vec_pca.begin(), vec_pca.end());

      histos.fill(HIST("AllMuon_pT"), mu_pt);
      if(PartType==1){
        histos.fill(HIST("AllBeauty_pT"), mu_pt);
      }else if(PartType==2){
        histos.fill(HIST("AllnpCharm_pT"), mu_pt);
      }else if(PartType==3){
        histos.fill(HIST("AllpCharm_pT"), mu_pt);
      }else if(PartType==4){
        histos.fill(HIST("AllBG_pT"), mu_pt);
      }
      
      bool c_DCATCut = false;
      bool b_DCATCut = false;
      if(PartType==1){
        if(muon_dcat>DCAT_min_c && muon_dcat<DCAT_max_c){ // Cut value for muons as charm 
          histos.fill(HIST("AsD_DCATCut_Beauty_pT"), mu_pt);
          histos.fill(HIST("AsD_DCATCut_Muon_pT"), mu_pt);
          c_DCATCut = true;
        }else if(muon_dcat>DCAT_min_b && muon_dcat<DCAT_max_b){ // Cut value for muons as beauty or non prompt charm
          histos.fill(HIST("AsB_DCATCut_Beauty_pT"), mu_pt);
          histos.fill(HIST("AsB_DCATCut_Muon_pT"), mu_pt);
          b_DCATCut = true;
        }
      }else if(PartType==2){
        if(muon_dcat>DCAT_min_c && muon_dcat<DCAT_max_c){ // Cut value for muons as charm 
          histos.fill(HIST("AsD_DCATCut_npCharm_pT"), mu_pt);
          histos.fill(HIST("AsD_DCATCut_Muon_pT"), mu_pt);
          c_DCATCut = true;
        }else if(muon_dcat>DCAT_min_b && muon_dcat<DCAT_max_b){ // Cut value for muons as beauty or non prompt charm
          histos.fill(HIST("AsB_DCATCut_npCharm_pT"), mu_pt);
          histos.fill(HIST("AsB_DCATCut_Muon_pT"), mu_pt);
          b_DCATCut = true;
        }
      }else if(PartType==3){
        if(muon_dcat>DCAT_min_c && muon_dcat<DCAT_max_c){ // Cut value for muons as charm 
          histos.fill(HIST("AsD_DCATCut_pCharm_pT"), mu_pt);
          histos.fill(HIST("AsD_DCATCut_Muon_pT"), mu_pt);
          c_DCATCut = true;
        }else if(muon_dcat>DCAT_min_b && muon_dcat<DCAT_max_b){ // Cut value for muons as beauty or non prompt charm
          histos.fill(HIST("AsB_DCATCut_pCharm_pT"), mu_pt);
          histos.fill(HIST("AsB_DCATCut_Muon_pT"), mu_pt);
          b_DCATCut = true;
        }
      }else if(PartType==4){
        if(muon_dcat>DCAT_min_c && muon_dcat<DCAT_max_c){ // Cut value for muons as charm 
          histos.fill(HIST("AsD_DCATCut_BG_pT"), mu_pt);
          histos.fill(HIST("AsD_DCATCut_Muon_pT"), mu_pt);
          c_DCATCut = true;
        }else if(muon_dcat>DCAT_min_b && muon_dcat<DCAT_max_b){ // Cut value for muons as beauty or non prompt charm
          histos.fill(HIST("AsB_DCATCut_BG_pT"), mu_pt);
          histos.fill(HIST("AsB_DCATCut_Muon_pT"), mu_pt);
          b_DCATCut = true;
        }
      }

      if(c_DCATCut==true){
        bool HasGoodMuon_asCharm = false;
        for(int i=0; i<vec_pca.size(); i++){
          if(vec_pca.at(i).at(10)<MFT_chi2){
            if(vec_pca.at(i).at(0)<PCA_cut){
              if(vec_pca.at(i).at(1)>CosSim_cut){
                if(vec_pca.at(i).at(6)<Total_dist_cut){
                  if(vec_pca.at(i).at(5)>PCAD_min_c && vec_pca.at(i).at(5)<PCAD_max_c){
                    HasGoodMuon_asCharm = true;
                  }
                }
              }
            }
          }
        }
        if(HasGoodMuon_asCharm==true){
          histos.fill(HIST("AsD_DCAT_PCACut_Muon_pT"), mu_pt);
          if(PartType==1){
            histos.fill(HIST("AsD_DCAT_PCACut_Beauty_pT"), mu_pt);
          }else if(PartType==2){
            histos.fill(HIST("AsD_DCAT_PCACut_npCharm_pT"), mu_pt);
          }else if(PartType==3){
            histos.fill(HIST("AsD_DCAT_PCACut_pCharm_pT"), mu_pt);
          }else if(PartType==4){
            histos.fill(HIST("AsD_DCAT_PCACut_BG_pT"), mu_pt);
          }
        }
      }else if(b_DCATCut==true){
        bool HasGoodMuon_asBeauty = false;
        for(int i=0; i<vec_pca.size(); i++){   
          if(vec_pca.at(i).at(10)<MFT_chi2){
            if(vec_pca.at(i).at(0)<PCA_cut){
              if(vec_pca.at(i).at(1)>CosSim_cut){
                if(vec_pca.at(i).at(6)<Total_dist_cut){
                  if(vec_pca.at(i).at(5)>PCAD_min_b && vec_pca.at(i).at(5)<PCAD_max_b){
                    HasGoodMuon_asBeauty = true;
                  }
                }
              }
            }
          }
        }
        if(HasGoodMuon_asBeauty==true){
          histos.fill(HIST("AsB_DCAT_PCACut_Muon_pT"), mu_pt);
          if(PartType==1){
            histos.fill(HIST("AsB_DCAT_PCACut_Beauty_pT"), mu_pt);
          }else if(PartType==2){
            histos.fill(HIST("AsB_DCAT_PCACut_npCharm_pT"), mu_pt);
          }else if(PartType==3){
            histos.fill(HIST("AsB_DCAT_PCACut_pCharm_pT"), mu_pt);
          }else if(PartType==4){
            histos.fill(HIST("AsB_DCAT_PCACut_BG_pT"), mu_pt);
          }
        }
      }
    }
  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyAnalysisTask>(cfgc)
  };
}