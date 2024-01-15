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

struct CutValueAnalysis{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<double> pT_cut{"pT_cut", 0.5, "Cut value for transverse momentum"};

  void init(InitContext const&){
    const AxisSpec axisDCA{10000, -5, 5, "cm"};
    const AxisSpec axisDCAT{6000, 0, 30, "cm"};
    const AxisSpec axispT{100, 0, 50, "p_{T} (GeV/c)"};
    const AxisSpec axisPCAR{5001, -0.0005, 5.0005, "PCAR(cm)"};
    const AxisSpec axisPCAD{1001, -0.005, 10.005, "PCAD(cm)"};
    const AxisSpec axisPCAZ{5001, -25.005, 25.005, "z(cm)"};
    const AxisSpec axisChi2{5010, -1, 500, "chi2"};
    const AxisSpec axisDist{3000, 0, 300, "cm"};
    const AxisSpec axisCounter{2, -0.5, 1.5, ""};
    const AxisSpec axisCosSim{500, 0.5, 1, "cosine_similarity"};

    histos.add("Beauty_DCAT", "Beauty_DCAT", kTH2F, {axisDCAT, axispT});
    histos.add("Beauty_Pair_DCAT", "Beauty_Pair_DCAT", kTH2F, {axisDCAT, axispT}); //pT is muon's pT
    histos.add("Beauty_Pair_Chi2", "Beauty_Pair_Chi2", kTH2F, {axisChi2, axispT}); //pT is muon's pT
    histos.add("Beauty_Pair_CosSim", "Beauty_Pair_CosSim", kTH1F, {axisCosSim});
    histos.add("Beauty_PCAR", "Beauty_PCAR", kTH1F, {axisPCAR});
    histos.add("Beauty_PCAZ", "Beauty_PCAZ", kTH1F, {axisPCAZ});
    histos.add("Beauty_PCAD", "Beauty_PCAD", kTH1F, {axisPCAD});
    histos.add("Beauty_TotalDist", "Beauty_TotalDist", kTH1F, {axisDist});
    histos.add("Beauty_Fake_Pair_Chi2", "Beauty_Fake_Pair_Chi2", kTH1F, {axisChi2});
    histos.add("Beauty_Fake_Pair_CosSim", "Beauty_Fake_Pair_CosSim", kTH1F, {axisCosSim});
    histos.add("Beauty_Fake_PCAR", "Beauty_Fake_PCAR", kTH1F, {axisPCAR});
    histos.add("Beauty_Fake_PCAZ", "Beauty_Fake_PCAZ", kTH1F, {axisPCAZ});
    histos.add("Beauty_Fake_PCAD", "Beauty_Fake_PCAD", kTH1F, {axisPCAD});
    histos.add("Beauty_Fake_TotalDist", "Beauty_Fake_TotalDist", kTH1F, {axisDist});
    histos.add("Beauty_PCA_Answer", "Beauty_PCA_Answer", kTH1F, {axisCounter});

    histos.add("pCharm_DCAT", "pCharm_DCAT", kTH2F, {axisDCAT, axispT});
    histos.add("pCharm_Pair_DCAT", "pCharm_Pair_DCAT", kTH2F, {axisDCAT, axispT}); //pT is muon's pT
    histos.add("pCharm_Pair_Chi2", "pCharm_Pair_Chi2", kTH2F, {axisChi2, axispT}); //pT is muon's pT
    histos.add("pCharm_Pair_CosSim", "pCharm_Pair_CosSim", kTH1F, {axisCosSim});
    histos.add("pCharm_PCAR", "pCharm_PCAR", kTH1F, {axisPCAR});
    histos.add("pCharm_PCAZ", "pCharm_PCAZ", kTH1F, {axisPCAZ});
    histos.add("pCharm_PCAD", "pCharm_PCAD", kTH1F, {axisPCAD});
    histos.add("pCharm_TotalDist", "pCharm_TotalDist", kTH1F, {axisDist});
    histos.add("pCharm_Fake_Pair_Chi2", "pCharm_Fake_Pair_Chi2", kTH1F, {axisChi2});
    histos.add("pCharm_Fake_Pair_CosSim", "pCharm_Fake_Pair_CosSim", kTH1F, {axisCosSim});
    histos.add("pCharm_Fake_PCAR", "pCharm_Fake_PCAR", kTH1F, {axisPCAR});
    histos.add("pCharm_Fake_PCAZ", "pCharm_Fake_PCAZ", kTH1F, {axisPCAZ});
    histos.add("pCharm_Fake_PCAD", "pCharm_Fake_PCAD", kTH1F, {axisPCAD});
    histos.add("pCharm_Fake_TotalDist", "pCharm_Fake_TotalDist", kTH1F, {axisDist});
    histos.add("pCharm_PCA_Answer", "pCharm_PCA_Answer", kTH1F, {axisCounter});

    histos.add("npCharm_DCAT", "npCharm_DCAT", kTH2F, {axisDCAT, axispT});
    histos.add("npCharm_Pair_DCAT", "npCharm_Pair_DCAT", kTH2F, {axisDCAT, axispT}); //pT is muon's pT
    histos.add("npCharm_Pair_Chi2", "npCharm_Pair_Chi2", kTH2F, {axisChi2, axispT}); //pT is muon's pT
    histos.add("npCharm_Pair_CosSim", "npCharm_Pair_CosSim", kTH1F, {axisCosSim});
    histos.add("npCharm_PCAR", "npCharm_PCAR", kTH1F, {axisPCAR});
    histos.add("npCharm_PCAZ", "npCharm_PCAZ", kTH1F, {axisPCAZ});
    histos.add("npCharm_PCAD", "npCharm_PCAD", kTH1F, {axisPCAD});
    histos.add("npCharm_TotalDist", "npCharm_TotalDist", kTH1F, {axisDist});
    histos.add("npCharm_Fake_Pair_Chi2", "npCharm_Fake_Pair_Chi2", kTH1F, {axisChi2});
    histos.add("npCharm_Fake_Pair_CosSim", "npCharm_Fake_Pair_CosSim", kTH1F, {axisCosSim});
    histos.add("npCharm_Fake_PCAR", "npCharm_Fake_PCAR", kTH1F, {axisPCAR});
    histos.add("npCharm_Fake_PCAZ", "npCharm_Fake_PCAZ", kTH1F, {axisPCAZ});
    histos.add("npCharm_Fake_PCAD", "npCharm_Fake_PCAD", kTH1F, {axisPCAD});
    histos.add("npCharm_Fake_TotalDist", "npCharm_Fake_TotalDist", kTH1F, {axisDist});
    histos.add("npCharm_PCA_Answer", "npCharm_PCA_Answer", kTH1F, {axisCounter});

    histos.add("LF_DCAT", "LF_DCAT", kTH2F, {axisDCAT, axispT});
    histos.add("LF_Pair_DCAT", "LF_Pair_DCAT", kTH2F, {axisDCAT, axispT});
    histos.add("LF_Pair_Chi2", "LF_Pair_Chi2", kTH2F, {axisChi2, axispT});
    histos.add("LF_Pair_CosSim", "LF_Pair_CosSim", kTH1F, {axisCosSim});
    histos.add("LF_PCAR", "LF_PCAR", kTH1F, {axisPCAR});
    histos.add("LF_PCAZ", "LF_PCAZ", kTH1F, {axisPCAZ});
    histos.add("LF_PCAD", "LF_PCAD", kTH1F, {axisPCAD});
    histos.add("LF_TotalDist", "LF_TotalDist", kTH1F, {axisDist});
    histos.add("LF_Fake_Pair_Chi2", "LF_Fake_Pair_Chi2", kTH1F, {axisChi2});
    histos.add("LF_Fake_Pair_CosSim", "LF_Fake_Pair_CosSim", kTH1F, {axisCosSim});
    histos.add("LF_Fake_PCAR", "LF_Fake_PCAR", kTH1F, {axisPCAR});
    histos.add("LF_Fake_PCAZ", "LF_Fake_PCAZ", kTH1F, {axisPCAZ});
    histos.add("LF_Fake_PCAD", "LF_Fake_PCAD", kTH1F, {axisPCAD});
    histos.add("LF_Fake_TotalDist", "LF_Fake_TotalDist", kTH1F, {axisDist});

    histos.add("BG_MFT_DCAT", "BG_MFT_DCAT", kTH1F, {axisDCAT});

    auto pca_ans_pcharm = histos.get<TH1>(HIST("pCharm_PCA_Answer"));
    auto* x = pca_ans_pcharm->GetXaxis();
    x->SetBinLabel(1, "Incorrect");
    x->SetBinLabel(2, "Correct");
    auto pca_ans_npcharm = histos.get<TH1>(HIST("npCharm_PCA_Answer"));
    auto* x1 = pca_ans_npcharm->GetXaxis();
    x1->SetBinLabel(1, "Incorrect");
    x1->SetBinLabel(2, "Correct");
    auto pca_ans_beauty = histos.get<TH1>(HIST("Beauty_PCA_Answer"));
    auto* x2 = pca_ans_beauty->GetXaxis();
    x2->SetBinLabel(1, "Incorrect");
    x2->SetBinLabel(2, "Correct");
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
      
      int muid = mcParticle_fwd.globalIndex();
      auto mumom = mcParticle_fwd.mothers_first_as<aod::McParticles>();
      auto Daughters = mumom.daughters_as<aod::McParticles>();
      double mu_pt = fwdtrack.pt();

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

      double secverx, secvery, secverz;
      vector<int> vec_daughter;
      vec_daughter.clear();
      bool HasPair = false;
      int lf_counter = 0;
      for(auto& Daughter : Daughters){
        int d_pdg = fabs(Daughter.pdgCode());
        if(d_pdg==13 && Daughter.globalIndex()==muid){
          secverx = Daughter.vx();
          secvery = Daughter.vy();
          secverz = Daughter.vz();
        }
        vec_daughter.push_back(Daughter.globalIndex());
      }
      if(lf_counter!=0) cout << endl;
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

      double mu_dcat = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
      if(PartType==1){
        histos.fill(HIST("Beauty_DCAT"), mu_dcat, mu_pt);
      }else if(PartType==2){
        histos.fill(HIST("npCharm_DCAT"), mu_dcat, mu_pt);
      }else if(PartType==3){
        histos.fill(HIST("pCharm_DCAT"), mu_dcat, mu_pt);
      }else if(PartType==4){
        histos.fill(HIST("LF_DCAT"), mu_dcat, mu_pt);
      }

      vector<double> mftv1;
      SMatrix55 mftcovs1(mftv1.begin(), mftv1.end());
      SMatrix5 muonpars(muon_parameter[0], muon_parameter[1], muon_parameter[3], muon_parameter[4], muon_parameter[5]);
      o2::track::TrackParCovFwd mftpars1{muon_parameter[2], muonpars, mftcovs1, muon_parameter[6]};

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
            double p1dcat = sqrt(pow(p1dcaX-Col_x,2)+pow(p1dcaY-Col_y,2));
            vec_mftdca = {mumftX, mumftY, mumftZ, mudcaX, mudcaY, p1mftX, p1mftY, p1mftZ, p1dcaX, p1dcaY, Col_x, Col_y, Col_z};
            vector<double> p1_pca = PCA_Cal(vec_mftdca);
            double add_dist = 0;
            for(int i=0; i<5000; i++){
              auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
              if((cut_z-Col_z)>0) break;
              mftpars1.propagateToZlinear(cut_z);
              mftpars2.propagateToZlinear(cut_z);
              auto diff = sqrt(pow(mftpars2.getX()-mftpars1.getX(),2)+pow(mftpars2.getY()-mftpars1.getY(),2));
              add_dist += diff;
            }
            p1_pca.push_back(add_dist);
            p1_pca.push_back(ans1);
            vec_pca.push_back(p1_pca);

            if(PartType==1){
              histos.fill(HIST("Beauty_Pair_Chi2"), p1chi2, mu_pt);
              histos.fill(HIST("Beauty_Pair_DCAT"), p1dcat, mu_pt);
              histos.fill(HIST("Beauty_Pair_CosSim"), p1_pca[1]);
              histos.fill(HIST("Beauty_PCAR"), p1_pca[0]);
              histos.fill(HIST("Beauty_PCAZ"), p1_pca[4]);
              histos.fill(HIST("Beauty_PCAD"), p1_pca[5]);
              histos.fill(HIST("Beauty_TotalDist"), p1_pca[6]);
            }else if(PartType==2){
              histos.fill(HIST("npCharm_Pair_Chi2"), p1chi2, mu_pt);
              histos.fill(HIST("npCharm_Pair_DCAT"), p1dcat, mu_pt);
              histos.fill(HIST("npCharm_Pair_CosSim"), p1_pca[1]);
              histos.fill(HIST("npCharm_PCAR"), p1_pca[0]);
              histos.fill(HIST("npCharm_PCAZ"), p1_pca[4]);
              histos.fill(HIST("npCharm_PCAD"), p1_pca[5]);
              histos.fill(HIST("npCharm_TotalDist"), p1_pca[6]);
            }else if(PartType==3){
              histos.fill(HIST("pCharm_Pair_Chi2"), p1chi2, mu_pt);
              histos.fill(HIST("pCharm_Pair_DCAT"), p1dcat, mu_pt);
              histos.fill(HIST("pCharm_Pair_CosSim"), p1_pca[1]);
              histos.fill(HIST("pCharm_PCAR"), p1_pca[0]);
              histos.fill(HIST("pCharm_PCAZ"), p1_pca[4]);
              histos.fill(HIST("pCharm_PCAD"), p1_pca[5]);
              histos.fill(HIST("pCharm_TotalDist"), p1_pca[6]);
            }else if(PartType==4){
              histos.fill(HIST("LF_Pair_Chi2"), p1chi2, mu_pt);
              histos.fill(HIST("LF_Pair_DCAT"), p1dcat, mu_pt);
              histos.fill(HIST("LF_Pair_CosSim"), p1_pca[1]);
              histos.fill(HIST("LF_PCAR"), p1_pca[0]);
              histos.fill(HIST("LF_PCAZ"), p1_pca[4]);
              histos.fill(HIST("LF_PCAD"), p1_pca[5]);
              histos.fill(HIST("LF_TotalDist"), p1_pca[6]);
            }

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
            double p2dcat = sqrt(pow(p2dcaX-Col_x,2)+pow(p2dcaY-Col_y,2));
            vec_mftdca = {mumftX, mumftY, mumftZ, mudcaX, mudcaY, p2mftX, p2mftY, p2mftZ, p2dcaX, p2dcaY, Col_x, Col_y, Col_z};
            vector<double> p2_pca = PCA_Cal(vec_mftdca);
            double add_dist = 0;
            for(int i=0; i<5000; i++){
              auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
              if((cut_z-Col_z)>0) break;
              mftpars1.propagateToZlinear(cut_z);
              mftpars3.propagateToZlinear(cut_z);
              auto diff = sqrt(pow(mftpars3.getX()-mftpars1.getX(),2)+pow(mftpars3.getY()-mftpars1.getY(),2));
              add_dist += diff;
            }
            p2_pca.push_back(add_dist);
            p2_pca.push_back(ans2);
            vec_pca.push_back(p2_pca);

            if(PartType==1){
              histos.fill(HIST("Beauty_Pair_Chi2"), p2chi2, mu_pt);
              histos.fill(HIST("Beauty_Pair_DCAT"), p2dcat, mu_pt);
              histos.fill(HIST("Beauty_Pair_CosSim"), p2_pca[1]);
              histos.fill(HIST("Beauty_PCAR"), p2_pca[0]);
              histos.fill(HIST("Beauty_PCAZ"), p2_pca[4]);
              histos.fill(HIST("Beauty_PCAD"), p2_pca[5]);
              histos.fill(HIST("Beauty_TotalDist"), p2_pca[6]);
            }else if(PartType==2){
              histos.fill(HIST("npCharm_Pair_Chi2"), p2chi2, mu_pt);
              histos.fill(HIST("npCharm_Pair_DCAT"), p2dcat, mu_pt);
              histos.fill(HIST("npCharm_Pair_CosSim"), p2_pca[1]);
              histos.fill(HIST("npCharm_PCAR"), p2_pca[0]);
              histos.fill(HIST("npCharm_PCAZ"), p2_pca[4]);
              histos.fill(HIST("npCharm_PCAD"), p2_pca[5]);
              histos.fill(HIST("npCharm_TotalDist"), p2_pca[6]);
            }else if(PartType==3){
              histos.fill(HIST("pCharm_Pair_Chi2"), p2chi2, mu_pt);
              histos.fill(HIST("pCharm_Pair_DCAT"), p2dcat, mu_pt);
              histos.fill(HIST("pCharm_Pair_CosSim"), p2_pca[1]);
              histos.fill(HIST("pCharm_PCAR"), p2_pca[0]);
              histos.fill(HIST("pCharm_PCAZ"), p2_pca[4]);
              histos.fill(HIST("pCharm_PCAD"), p2_pca[5]);
              histos.fill(HIST("pCharm_TotalDist"), p2_pca[6]);
            }else if(PartType==4){
              histos.fill(HIST("LF_Pair_Chi2"), p2chi2, mu_pt);
              histos.fill(HIST("LF_Pair_DCAT"), p2dcat, mu_pt);
              histos.fill(HIST("LF_Pair_CosSim"), p2_pca[1]);
              histos.fill(HIST("LF_PCAR"), p2_pca[0]);
              histos.fill(HIST("LF_PCAZ"), p2_pca[4]);
              histos.fill(HIST("LF_PCAD"), p2_pca[5]);
              histos.fill(HIST("LF_TotalDist"), p2_pca[6]);
            }

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
            double p3dcat = sqrt(pow(p3dcaX-Col_x,2)+pow(p3dcaY-Col_y,2));
            vec_mftdca = {mumftX, mumftY, mumftZ, mudcaX, mudcaY, p3mftX, p3mftY, p3mftZ, p3dcaX, p3dcaY, Col_x, Col_y, Col_z};
            vector<double> p3_pca = PCA_Cal(vec_mftdca);
            double add_dist = 0;
            for(int i=0; i<5000; i++){
              auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
              if((cut_z-Col_z)>0) break;
              mftpars1.propagateToZlinear(cut_z);
              mftpars4.propagateToZlinear(cut_z);
              auto diff = sqrt(pow(mftpars4.getX()-mftpars1.getX(),2)+pow(mftpars4.getY()-mftpars1.getY(),2));
              add_dist += diff;
            }
            p3_pca.push_back(add_dist);
            p3_pca.push_back(ans3);
            vec_pca.push_back(p3_pca);

            if(PartType==1){
              histos.fill(HIST("Beauty_Pair_Chi2"), p3chi2, mu_pt);
              histos.fill(HIST("Beauty_Pair_DCAT"), p3dcat, mu_pt);
              histos.fill(HIST("Beauty_Pair_CosSim"), p3_pca[1]);
              histos.fill(HIST("Beauty_PCAR"), p3_pca[0]);
              histos.fill(HIST("Beauty_PCAZ"), p3_pca[4]);
              histos.fill(HIST("Beauty_PCAD"), p3_pca[5]);
              histos.fill(HIST("Beauty_TotalDist"), p3_pca[6]);
            }else if(PartType==2){
              histos.fill(HIST("npCharm_Pair_Chi2"), p3chi2, mu_pt);
              histos.fill(HIST("npCharm_Pair_DCAT"), p3dcat, mu_pt);
              histos.fill(HIST("npCharm_Pair_CosSim"), p3_pca[1]);
              histos.fill(HIST("npCharm_PCAR"), p3_pca[0]);
              histos.fill(HIST("npCharm_PCAZ"), p3_pca[4]);
              histos.fill(HIST("npCharm_PCAD"), p3_pca[5]);
              histos.fill(HIST("npCharm_TotalDist"), p3_pca[6]);
            }else if(PartType==3){
              histos.fill(HIST("pCharm_Pair_Chi2"), p3chi2, mu_pt);
              histos.fill(HIST("pCharm_Pair_DCAT"), p3dcat, mu_pt);
              histos.fill(HIST("pCharm_Pair_CosSim"), p3_pca[1]);
              histos.fill(HIST("pCharm_PCAR"), p3_pca[0]);
              histos.fill(HIST("pCharm_PCAZ"), p3_pca[4]);
              histos.fill(HIST("pCharm_PCAD"), p3_pca[5]);
              histos.fill(HIST("pCharm_TotalDist"), p3_pca[6]);
            }else if(PartType==4){
              histos.fill(HIST("LF_Pair_Chi2"), p3chi2, mu_pt);
              histos.fill(HIST("LF_Pair_DCAT"), p3dcat, mu_pt);
              histos.fill(HIST("LF_Pair_CosSim"), p3_pca[1]);
              histos.fill(HIST("LF_PCAR"), p3_pca[0]);
              histos.fill(HIST("LF_PCAZ"), p3_pca[4]);
              histos.fill(HIST("LF_PCAD"), p3_pca[5]);
              histos.fill(HIST("LF_TotalDist"), p3_pca[6]);
            }

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
            double bgdcat = sqrt(pow(bgdcaX-Col_x,2)+pow(bgdcaY-Col_y,2));
            histos.fill(HIST("BG_MFT_DCAT"), bgdcat);
            vec_mftdca = {mumftX, mumftY, mumftZ, mudcaX, mudcaY, bgmftX, bgmftY, bgmftZ, bgdcaX, bgdcaY, Col_x, Col_y, Col_z};
            vector<double> bg_pca = PCA_Cal(vec_mftdca);
            double add_dist = 0;
            for(int i=0; i<5000; i++){
              auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
              if((cut_z-Col_z)>0) break;
              mftpars1.propagateToZlinear(cut_z);
              mftpars0.propagateToZlinear(cut_z);
              auto diff = sqrt(pow(mftpars0.getX()-mftpars1.getX(),2)+pow(mftpars0.getY()-mftpars1.getY(),2));
              add_dist += diff;
            }
            bg_pca.push_back(add_dist);
            bg_pca.push_back(ans0);
            vec_pca.push_back(bg_pca);

            if(PartType==1){
              histos.fill(HIST("Beauty_Fake_Pair_Chi2"), bgchi2);
              histos.fill(HIST("Beauty_Fake_Pair_CosSim"), bg_pca[1]);
              histos.fill(HIST("Beauty_Fake_PCAR"), bg_pca[0]);
              histos.fill(HIST("Beauty_Fake_PCAZ"), bg_pca[4]);
              histos.fill(HIST("Beauty_Fake_PCAD"), bg_pca[5]);
              histos.fill(HIST("Beauty_Fake_TotalDist"), bg_pca[6]);
            }else if(PartType==2){
              histos.fill(HIST("npCharm_Fake_Pair_Chi2"), bgchi2);
              histos.fill(HIST("npCharm_Fake_Pair_CosSim"), bg_pca[1]);
              histos.fill(HIST("npCharm_Fake_PCAR"), bg_pca[0]);
              histos.fill(HIST("npCharm_Fake_PCAZ"), bg_pca[4]);
              histos.fill(HIST("npCharm_Fake_PCAD"), bg_pca[5]);
              histos.fill(HIST("npCharm_Fake_TotalDist"), bg_pca[6]);
            }else if(PartType==3){
              histos.fill(HIST("pCharm_Fake_Pair_Chi2"), bgchi2);
              histos.fill(HIST("pCharm_Fake_Pair_CosSim"), bg_pca[1]);
              histos.fill(HIST("pCharm_Fake_PCAR"), bg_pca[0]);
              histos.fill(HIST("pCharm_Fake_PCAZ"), bg_pca[4]);
              histos.fill(HIST("pCharm_Fake_PCAD"), bg_pca[5]);
              histos.fill(HIST("pCharm_Fake_TotalDist"), bg_pca[6]);
            }else if(PartType==4){
              histos.fill(HIST("LF_Fake_Pair_Chi2"), bgchi2);
              histos.fill(HIST("LF_Fake_Pair_CosSim"), bg_pca[1]);
              histos.fill(HIST("LF_Fake_PCAR"), bg_pca[0]);
              histos.fill(HIST("LF_Fake_PCAZ"), bg_pca[4]);
              histos.fill(HIST("LF_Fake_PCAD"), bg_pca[5]);
              histos.fill(HIST("LF_Fake_TotalDist"), bg_pca[6]);
            }

          }
        }
      }
      if(vec_pca.empty()) continue;
      if(HasPair1Track==true || HasPair2Track==true || HasPair3Track==true){
        sort(vec_pca.begin(), vec_pca.end());
        if(vec_pca.at(0).at(7)==0){
          if(PartType==1){
            histos.fill(HIST("Beauty_PCA_Answer"), 0);
          }else if(PartType==2){
            histos.fill(HIST("npCharm_PCA_Answer"), 0);
          }else if(PartType==3){
            histos.fill(HIST("pCharm_PCA_Answer"), 0);
          }
        }else{
          if(PartType==1){
            histos.fill(HIST("Beauty_PCA_Answer"), 1);
          }else if(PartType==2){
            histos.fill(HIST("npCharm_PCA_Answer"), 1);
          }else if(PartType==3){
            histos.fill(HIST("pCharm_PCA_Answer"), 1);
          }
        }
      }
    }
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CutValueAnalysis>(cfgc)
  };
}