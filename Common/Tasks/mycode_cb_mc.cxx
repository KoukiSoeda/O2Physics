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

  if(pcaZ>0){
    r_xyz = 10000;
    for(int z_mu=0; z_mu<100; z_mu++){
      predict_mux = mudcaX + (z_mu*0.001)*Nax;
      predict_muy = mudcaY + (z_mu*0.001)*Nay;
      predict_muz = Col_z + (z_mu*0.001)*Naz;
      for(int z_pair=0; z_pair<100; z_pair++){
        predict_canx = candcaX + (z_pair*0.001)*Ncx;
        predict_cany = candcaY + (z_pair*0.001)*Ncy;
        predict_canz = Col_z + (z_pair*0.001)*Ncz;
        double r_xyz_tmp = sqrt(pow(predict_canx-predict_mux,2) + pow(predict_cany-predict_muy,2) + pow(predict_canz-predict_muz,2));
        if(r_xyz>r_xyz_tmp){
          r_xyz=r_xyz_tmp;
          pcaX = ((predict_mux+predict_canx)/2)-Col_x;
          pcaY = ((predict_muy+predict_cany)/2)-Col_y;
          pcaZ = ((predict_muz+predict_canz)/2)-Col_z;
        }
      }
    }
  }
  auto pcaD = sqrt(pow(pcaX,2)+pow(pcaY,2)+pow(pcaZ,2));
  double dcat = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
  vec_pca = {r_xyz, cosxy, pcaX, pcaY, pcaZ, pcaD, dcat};
  return vec_pca;
}

struct McCharmBeauty{
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

  void init(InitContext const&){
    const AxisSpec axisPCAR{10010, -0.0005, 1.0005, "PCAR(cm)"};
    const AxisSpec axisPCAD{1001, -0.005, 10.005, "PCAD(cm)"};
    const AxisSpec axisPCAZ{5001, -25.005, 25.005, "z(cm)"};
    const AxisSpec axisDCA{10000, -5, 5, "cm"};
    const AxisSpec axisDCAT{30000, 0, 30, "cm"};
    const AxisSpec axispT{100, 0, 50, "p_{T}(GeV/c)"};
    const AxisSpec axisp{1000, 0, 100, "p(GeV/c)"};
    const AxisSpec axisCounter{2, -0.5, 1.5, ""};
    const AxisSpec axisDeltaX{4000, -1, 1, "#Delta(x_{mc}-x_{det})"};
    const AxisSpec axisDeltaY{4000, -1, 1, "#Delta(y_{mc}-y_{det})"};
    const AxisSpec axisDeltaZ{4000, -1, 1, "#Delta(z_{mc}-z_{det})"};
    const AxisSpec axisMultiplicity{250, -0.5, 250.5, "Multiplicity"};

    histos.add("Beauty_PCAR", "Beauty_PCAR", kTH1F, {axisPCAR});
    histos.add("Beauty_PCAD", "Beauty_PCAD", kTH1F, {axisPCAD});
    histos.add("Beauty_PCAZ", "Beauty_PCAZ", kTH1F, {axisPCAZ});
    histos.add("Beauty_Fake_PCAR", "Beauty_Fake_PCAR", kTH1F, {axisPCAR});
    histos.add("Beauty_Fake_PCAD", "Beauty_Fake_PCAD", kTH1F, {axisPCAD});
    histos.add("Beauty_Fake_PCAZ", "Beauty_Fake_PCAZ", kTH1F, {axisPCAZ});
    histos.add("Beauty_DeltaPCAX", "Beauty_DeltaPCAX", kTH1F, {axisDeltaX});
    histos.add("Beauty_DeltaPCAY", "Beauty_DeltaPCAY", kTH1F, {axisDeltaY});
    histos.add("Beauty_DeltaPCAZ", "Beauty_DeltaPCAZ", kTH1F, {axisDeltaZ});
    histos.add("Beauty_PCA_Answer", "Beauty_PCA_Answer", kTH1F, {axisCounter});
    histos.add("Beauty_DCAx", "Beauty_DCAx", kTH1F, {axisDCA});
    histos.add("Beauty_DCAy", "Beauty_DCAy", kTH1F, {axisDCA});
    histos.add("Beauty_DCAT", "Beauty_DCAT", kTH1F, {axisDCAT});
    histos.add("Beauty_Muon_pT", "Beauty_Muon_pT", kTH1F, {axispT});
    histos.add("Beauty_Muon_p", "Beauty_Muon_p", kTH1F, {axisp});

    histos.add("pCharm_PCAR", "pCharm_PCAR", kTH1F, {axisPCAR});
    histos.add("pCharm_PCAD", "pCharm_PCAD", kTH1F, {axisPCAD});
    histos.add("pCharm_PCAZ", "pCharm_PCAZ", kTH1F, {axisPCAZ});
    histos.add("pCharm_Fake_PCAR", "pCharm_Fake_PCAR", kTH1F, {axisPCAR});
    histos.add("pCharm_Fake_PCAD", "pCharm_Fake_PCAD", kTH1F, {axisPCAD});
    histos.add("pCharm_Fake_PCAZ", "pCharm_Fake_PCAZ", kTH1F, {axisPCAZ});
    histos.add("pCharm_DeltaPCAX", "pCharm_DeltaPCAX", kTH1F, {axisDeltaX});
    histos.add("pCharm_DeltaPCAY", "pCharm_DeltaPCAY", kTH1F, {axisDeltaY});
    histos.add("pCharm_DeltaPCAZ", "pCharm_DeltaPCAZ", kTH1F, {axisDeltaZ});
    histos.add("pCharm_PCA_Answer", "pCharm_PCA_Answer", kTH1F, {axisCounter});
    histos.add("pCharm_DCAx", "pCharm_DCAx", kTH1F, {axisDCA});
    histos.add("pCharm_DCAy", "pCharm_DCAy", kTH1F, {axisDCA});
    histos.add("pCharm_DCAT", "pCharm_DCAT", kTH1F, {axisDCAT});
    histos.add("pCharm_Muon_pT", "pCharm_Muon_pT", kTH1F, {axispT});
    histos.add("pCharm_Muon_p", "pCharm_Muon_p", kTH1F, {axisp});

    histos.add("npCharm_PCAR", "npCharm_PCAR", kTH1F, {axisPCAR});
    histos.add("npCharm_PCAD", "npCharm_PCAD", kTH1F, {axisPCAD});
    histos.add("npCharm_PCAZ", "npCharm_PCAZ", kTH1F, {axisPCAZ});
    histos.add("npCharm_Fake_PCAR", "npCharm_Fake_PCAR", kTH1F, {axisPCAR});
    histos.add("npCharm_Fake_PCAD", "npCharm_Fake_PCAD", kTH1F, {axisPCAD});
    histos.add("npCharm_Fake_PCAZ", "npCharm_Fake_PCAZ", kTH1F, {axisPCAZ});
    histos.add("npCharm_DeltaPCAX", "npCharm_DeltaPCAX", kTH1F, {axisDeltaX});
    histos.add("npCharm_DeltaPCAY", "npCharm_DeltaPCAY", kTH1F, {axisDeltaY});
    histos.add("npCharm_DeltaPCAZ", "npCharm_DeltaPCAZ", kTH1F, {axisDeltaZ});
    histos.add("npCharm_PCA_Answer", "npCharm_PCA_Answer", kTH1F, {axisCounter});
    histos.add("npCharm_DCAx", "npCharm_DCAx", kTH1F, {axisDCA});
    histos.add("npCharm_DCAy", "npCharm_DCAy", kTH1F, {axisDCA});
    histos.add("npCharm_DCAT", "npCharm_DCAT", kTH1F, {axisDCAT});
    histos.add("npCharm_Muon_pT", "npCharm_Muon_pT", kTH1F, {axispT});
    histos.add("npCharm_Muon_p", "npCharm_Muon_p", kTH1F, {axisp});

    histos.add("BG_PCAR", "BG_PCAR", kTH1F, {axisPCAR});
    histos.add("BG_PCAD", "BG_PCAD", kTH1F, {axisPCAD});
    histos.add("BG_PCAZ", "BG_PCAZ", kTH1F, {axisPCAZ});

    histos.add("nPair_Multiplicity", "nPair_Multiplicity", kTH1F, {axisMultiplicity});

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
    const double mftz = -46.0; // cm
    for(auto& fwdtrack : fwdtracks){
      if(!fwdtrack.has_mcParticle()) continue;
      if(fwdtrack.trackType()==1 || fwdtrack.trackType()==3 || fwdtrack.trackType()==4) continue;
      if(fwdtrack.eta()<=(-3.6) || fwdtrack.eta()>=(-2.5)) continue;
      if(fwdtrack.pt()<pT_cut) continue;

      auto mcParticle_fwd = fwdtrack.mcParticle();
      int mcCol_id = -1;
      double Col_x, Col_y, Col_z;
      double mcCol_x, mcCol_y, mcCol_z;
      bool Has_MCCol = false;
      for(auto& t_col : collisions){
        if(!t_col.has_mcCollision()) continue;
        if(mcParticle_fwd.mcCollisionId()==t_col.mcCollisionId()){
          mcCol_id = t_col.mcCollisionId();
          Col_x = t_col.posX();
          Col_y = t_col.posY();
          Col_z = t_col.posZ();
          mcCol_x = t_col.mcCollision().posX();
          mcCol_y = t_col.mcCollision().posY();
          mcCol_z = t_col.mcCollision().posZ();
          Has_MCCol = true;
          break;
        }
      }
      if(fabs(mcCol_z)>10 || Has_MCCol==false) continue;
      if(fabs(mcParticle_fwd.pdgCode())!=13 || !mcParticle_fwd.has_mothers()) continue;

      int muid = mcParticle_fwd.globalIndex();
      double mu_pt = mcParticle_fwd.pt();
      double mu_p = mcParticle_fwd.p();
      auto mumom = mcParticle_fwd.mothers_first_as<aod::McParticles>();
      auto Daughters = mumom.daughters_as<aod::McParticles>();
      double convert_phi = 0;
      if(mcParticle_fwd.phi()>o2::math_utils::pi() && mcParticle_fwd.phi()<2*o2::math_utils::pi()){
        convert_phi = mcParticle_fwd.phi() - 2*(o2::math_utils::pi());
      }else{
        convert_phi = mcParticle_fwd.phi();
      }
      double mu_signedpt = mcParticle_fwd.pdgCode()/(fabs(mcParticle_fwd.pdgCode())*mcParticle_fwd.pt());

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

      double secverx, secvery, secverz, mueta;
      vector<int> vec_daughter;
      vec_daughter.clear();
      bool HasPair = false;
      for(auto& Daughter : Daughters){
        int d_pdg = fabs(Daughter.pdgCode());
        if(d_pdg==13 && Daughter.globalIndex()==muid){
          secverx = Daughter.vx();
          secvery = Daughter.vy();
          secverz = Daughter.vz();
          mueta = Daughter.eta();
        }
        vec_daughter.push_back(Daughter.globalIndex());
      }
      if(vec_daughter.size()>1) HasPair = true;
      //if(!HasPair) continue;

      vector<double> muon_vec;
      SMatrix55 muoncovs(muon_vec.begin(), muon_vec.end());
      SMatrix5 muonpars(secverx, secvery, convert_phi, tan((o2::math_utils::pi()/2)-2*atan(exp(-mueta))), mu_signedpt);
      o2::track::TrackParCovFwd muonpropa(secverz, muonpars, muoncovs, fwdtrack.chi2());
      muonpropa.propagateToZlinear(mftz);
      double mumftX = muonpropa.getX();
      double mumftY = muonpropa.getY();
      muonpropa.propagateToZlinear(mcCol_z);
      double mudcaX = muonpropa.getX();
      double mudcaY = muonpropa.getY();

      double dcaX = mudcaX - mcCol_x;
      double dcaY = mudcaY - mcCol_y;
      double dcat = sqrt(pow(dcaX,2)+pow(dcaY,2));

      bool Is_LFmuon = false;
      double bg_x, bg_y, bg_z, bg_phi, bg_tgl, bg_signedpt, bg_chi2;
      if(PartType==1){
        histos.fill(HIST("Beauty_DCAx"), dcaX);
        histos.fill(HIST("Beauty_DCAy"), dcaY);
        histos.fill(HIST("Beauty_DCAT"), dcat);
        histos.fill(HIST("Beauty_Muon_pT"), mu_pt);
        histos.fill(HIST("Beauty_Muon_p"), mu_p);
      }else if(PartType==2){
        histos.fill(HIST("npCharm_DCAx"), dcaX);
        histos.fill(HIST("npCharm_DCAy"), dcaY);
        histos.fill(HIST("npCharm_DCAT"), dcat);
        histos.fill(HIST("npCharm_Muon_pT"), mu_pt);
        histos.fill(HIST("npCharm_Muon_p"), mu_p);
      }else if(PartType==3){
        histos.fill(HIST("pCharm_DCAx"), dcaX);
        histos.fill(HIST("pCharm_DCAy"), dcaY);
        histos.fill(HIST("pCharm_DCAT"), dcat);
        histos.fill(HIST("pCharm_Muon_pT"), mu_pt);
        histos.fill(HIST("pCharm_Muon_p"), mu_p);
      }else if(PartType==4){
        Is_LFmuon = true;
        bg_x = fwdtrack.x();
        bg_y = fwdtrack.y();
        bg_z = fwdtrack.z();
        bg_phi = fwdtrack.phi();
        bg_tgl = fwdtrack.tgl();
        bg_signedpt = fwdtrack.signed1Pt();
        bg_chi2 = fwdtrack.chi2();
      }
      double bgdcaX=0, bgdcaY=0;
      if(Is_LFmuon==true){
        vector<double> bg_vec;
        SMatrix55 bgcovs(bg_vec.begin(), bg_vec.end());
        SMatrix5 bgpars(bg_x, bg_y, bg_phi, bg_tgl, bg_signedpt);
        o2::track::TrackParCovFwd bgpropa(bg_z, bgpars, bgcovs, bg_chi2);
        bgpropa.propagateToZlinear(mcCol_z);
        bgdcaX = bgpropa.getX();
        bgdcaY = bgpropa.getY();
      }

      bool HasMuonTrack=false, HasPair1Track=false, HasPair2Track=false, HasPair3Track=false;
      double p1mftX, p1mftY, p1dcaX, p1dcaY, p1chi2, p1phi, p1tgl, p1signedpt;
      double p2mftX, p2mftY, p2dcaX, p2dcaY, p2chi2, p2phi, p2tgl, p2signedpt;
      double p3mftX, p3mftY, p3dcaX, p3dcaY, p3chi2, p3phi, p3tgl, p3signedpt;
      double canmftX, canmftY, candcaX, candcaY, canchi2, canphi, cantgl, cansignedpt;
      vector<double> vec_pca_p1;
      vector<double> vec_pca_p2;
      vector<double> vec_pca_p3;
      vector<vector<double>> vec_pca;
      vec_pca_p1.clear();
      vec_pca_p2.clear();
      vec_pca_p3.clear();
      vec_pca.clear();

      int npair = 0;
      for(auto& mfttrack : mfttracks){
        if(!mfttrack.has_mcParticle()) continue;
        if(mfttrack.mcParticle().mcCollisionId()!=mcCol_id) continue;
        auto mcParticle_mft = mfttrack.mcParticle();
        int mftid = mcParticle_mft.globalIndex();
        vector<double> vec_mftdca;
        vec_mftdca.clear();

        if(mftid==muid){
          HasMuonTrack = true;
          continue;
        }

        double mc_sign = 0;
        npair++;
        if(find(vec_daughter.begin(), vec_daughter.end(), mftid)!=vec_daughter.end() && HasPair1Track==false){
          HasPair1Track = true;
          double ans1 = 1;
          p1chi2 = mfttrack.chi2();
          double tmp_phi = mcParticle_mft.phi();
          if(tmp_phi>o2::math_utils::pi() && tmp_phi<2*o2::math_utils::pi()){
            p1phi = tmp_phi - 2*(o2::math_utils::pi());
          }else{
            p1phi = tmp_phi;
          }
          p1tgl = tan((o2::math_utils::pi()/2)-2*atan(exp(-mcParticle_mft.eta())));
          if(mcParticle_mft.pdgCode()>0){
            mc_sign = 1;
          }else if(mcParticle_mft.pdgCode()<0){
            mc_sign = -1;
          }
          p1signedpt = mc_sign/mcParticle_mft.pt();

          vector<double> p1_vec;
          SMatrix55 p1covs(p1_vec.begin(), p1_vec.end());
          SMatrix5 p1pars(mcParticle_mft.vx(), mcParticle_mft.vy(), p1phi, p1tgl, p1signedpt);
          o2::track::TrackParCovFwd p1propa(mcParticle_mft.vz(), p1pars, p1covs, p1chi2);
          p1propa.propagateToZlinear(mftz);
          p1mftX = p1propa.getX();
          p1mftY = p1propa.getY();
          p1propa.propagateToZlinear(mcCol_z);
          p1dcaX = p1propa.getX();
          p1dcaY = p1propa.getY();
          vec_mftdca = {mumftX, mumftY, mftz, mudcaX, mudcaY, p1mftX, p1mftY, mftz, p1dcaX, p1dcaY, mcCol_x, mcCol_y, mcCol_z};
          vec_pca_p1 = PCA_Cal(vec_mftdca);
          vec_pca_p1.push_back(ans1);
          vec_pca.push_back(vec_pca_p1);

        }else if(find(vec_daughter.begin(), vec_daughter.end(), mftid)!=vec_daughter.end() && HasPair1Track==true && HasPair2Track==false){
          HasPair2Track = true;
          double ans2 = 2;
          p2chi2 = mfttrack.chi2();
          double tmp_phi = mcParticle_mft.phi();
          if(tmp_phi>o2::math_utils::pi() && tmp_phi<2*o2::math_utils::pi()){
            p2phi = tmp_phi - 2*(o2::math_utils::pi());
          }else{
            p2phi = tmp_phi;
          }
          p2tgl = tan((o2::math_utils::pi()/2)-2*atan(exp(-mcParticle_mft.eta())));
          if(mcParticle_mft.pdgCode()>0){
            mc_sign = 1;
          }else if(mcParticle_mft.pdgCode()<0){
            mc_sign = -1;
          }
          p2signedpt = mc_sign/mcParticle_mft.pt();

          vector<double> p2_vec;
          SMatrix55 p2covs(p2_vec.begin(), p2_vec.end());
          SMatrix5 p2pars(mcParticle_mft.vx(), mcParticle_mft.vy(), p2phi, p2tgl, p2signedpt);
          o2::track::TrackParCovFwd p2propa(mcParticle_mft.vz(), p2pars, p2covs, p2chi2);
          p2propa.propagateToZlinear(mftz);
          p2mftX = p2propa.getX();
          p2mftY = p2propa.getY();
          p2propa.propagateToZlinear(mcCol_z);
          p2dcaX = p2propa.getX();
          p2dcaY = p2propa.getY();
          vec_mftdca = {mumftX, mumftY, mftz, mudcaX, mudcaY, p2mftX, p2mftY, mftz, p2dcaX, p2dcaY, mcCol_x, mcCol_y, mcCol_z};
          vec_pca_p2 = PCA_Cal(vec_mftdca);
          vec_pca_p2.push_back(ans2);
          vec_pca.push_back(vec_pca_p2);

        }else if(find(vec_daughter.begin(), vec_daughter.end(), mftid)!=vec_daughter.end() && HasPair1Track==true && HasPair2Track==true && HasPair3Track==false){
          HasPair3Track = true;
          double ans3 = 3;
          p3chi2 = mfttrack.chi2();
          double tmp_phi = mcParticle_mft.phi();
          if(tmp_phi>o2::math_utils::pi() && tmp_phi<2*o2::math_utils::pi()){
            p3phi = tmp_phi - 2*(o2::math_utils::pi());
          }else{
            p3phi = tmp_phi;
          }
          p3tgl = tan((o2::math_utils::pi()/2)-2*atan(exp(-mcParticle_mft.eta())));
          if(mcParticle_mft.pdgCode()>0){
            mc_sign = 1;
          }else if(mcParticle_mft.pdgCode()<0){
            mc_sign = -1;
          }
          p3signedpt = mc_sign/mcParticle_mft.pt();

          vector<double> p3_vec;
          SMatrix55 p3covs(p3_vec.begin(), p3_vec.end());
          SMatrix5 p3pars(mcParticle_mft.vx(), mcParticle_mft.vy(), p3phi, p3tgl, p3signedpt);
          o2::track::TrackParCovFwd p3propa(mcParticle_mft.vz(), p3pars, p3covs, p3chi2);
          p3propa.propagateToZlinear(mftz);
          p3mftX = p3propa.getX();
          p3mftY = p3propa.getY();
          p3propa.propagateToZlinear(mcCol_z);
          p3dcaX = p3propa.getX();
          p3dcaY = p3propa.getY();
          vec_mftdca = {mumftX, mumftY, mftz, mudcaX, mudcaY, p3mftX, p3mftY, mftz, p3dcaX, p3dcaY, mcCol_x, mcCol_y, mcCol_z};
          vec_pca_p3 = PCA_Cal(vec_mftdca);
          vec_pca_p3.push_back(ans3);
          vec_pca.push_back(vec_pca_p3);

        }else if(find(vec_daughter.begin(), vec_daughter.end(), mftid)==vec_daughter.end() && HasPair==true && Is_LFmuon==false){
          canchi2 = mfttrack.chi2();
          double ans0 = 0;
          double tmp_phi = mcParticle_mft.phi();
          if(tmp_phi>o2::math_utils::pi() && tmp_phi<2*o2::math_utils::pi()){
            canphi = tmp_phi - 2*(o2::math_utils::pi());
          }else{
            canphi = tmp_phi;
          }
          cantgl = tan((o2::math_utils::pi()/2)-2*atan(exp(-mcParticle_mft.eta())));
          if(mcParticle_mft.pdgCode()>0){
            mc_sign = 1;
          }else if(mcParticle_mft.pdgCode()<0){
            mc_sign = -1;
          }
          cansignedpt = mc_sign/mcParticle_mft.pt();

          vector<double> can_vec;
          SMatrix55 cancovs(can_vec.begin(), can_vec.end());
          SMatrix5 canpars(mcParticle_mft.vx(), mcParticle_mft.vy(), canphi, cantgl, cansignedpt);
          o2::track::TrackParCovFwd canpropa(mcParticle_mft.vz(), canpars, cancovs, canchi2);
          canpropa.propagateToZlinear(mftz);
          canmftX = canpropa.getX();
          canmftY = canpropa.getY();
          canpropa.propagateToZlinear(mcCol_z);
          candcaX = canpropa.getX();
          candcaY = canpropa.getY();
          vec_mftdca = {mumftX, mumftY, mftz, mudcaX, mudcaY, canmftX, canmftY, mftz, candcaX, candcaY, mcCol_x, mcCol_y, mcCol_z};
          vector<double> vec_pca_can;
          vec_pca_can.clear();
          vec_pca_can = PCA_Cal(vec_mftdca);
          vec_pca_can.push_back(ans0);
          vec_pca.push_back(vec_pca_can);
        }else if(Is_LFmuon==true){
          canchi2 = mfttrack.chi2();
          double ans0 = 0;
          double tmp_phi = mcParticle_mft.phi();
          if(tmp_phi>o2::math_utils::pi() && tmp_phi<2*o2::math_utils::pi()){
            canphi = tmp_phi - 2*(o2::math_utils::pi());
          }else{
            canphi = tmp_phi;
          }
          cantgl = tan((o2::math_utils::pi()/2)-2*atan(exp(-mcParticle_mft.eta())));
          if(mcParticle_mft.pdgCode()>0){
            mc_sign = 1;
          }else if(mcParticle_mft.pdgCode()<0){
            mc_sign = -1;
          }
          cansignedpt = mc_sign/mcParticle_mft.pt();

          vector<double> can_vec;
          SMatrix55 cancovs(can_vec.begin(), can_vec.end());
          SMatrix5 canpars(mcParticle_mft.vx(), mcParticle_mft.vy(), canphi, cantgl, cansignedpt);
          o2::track::TrackParCovFwd canpropa(mcParticle_mft.vz(), canpars, cancovs, canchi2);
          canpropa.propagateToZlinear(mftz);
          canmftX = canpropa.getX();
          canmftY = canpropa.getY();
          canpropa.propagateToZlinear(mcCol_z);
          candcaX = canpropa.getX();
          candcaY = canpropa.getY();
          vec_mftdca = {bg_x, bg_y, bg_z, bgdcaX, bgdcaY, canmftX, canmftY, mftz, candcaX, candcaY, mcCol_x, mcCol_y, mcCol_z};
          vector<double> vec_pca_can;
          vec_pca_can.clear();
          vec_pca_can = PCA_Cal(vec_mftdca);
          vec_pca_can.push_back(ans0);
          vec_pca.push_back(vec_pca_can);
        }
      }
      if(vec_pca.empty()) continue;
      histos.fill(HIST("nPair_Multiplicity"), npair);
      sort(vec_pca.begin(), vec_pca.end());

      if(PartType==4){
        histos.fill(HIST("BG_PCAR"), vec_pca.at(0).at(0));
        histos.fill(HIST("BG_PCAD"), vec_pca.at(0).at(5));
        histos.fill(HIST("BG_PCAZ"), vec_pca.at(0).at(4));
      }

      if(HasMuonTrack==true){
        if(HasPair1Track==true && HasPair2Track==false){
          double delta_x = (vec_pca_p1[2]+mcCol_x)-secverx;
          double delta_y = (vec_pca_p1[3]+mcCol_y)-secvery;
          double delta_z = (vec_pca_p1[4]+mcCol_z)-secverz;
          if(PartType==1){
            histos.fill(HIST("Beauty_PCAR"), vec_pca_p1[0]);
            histos.fill(HIST("Beauty_PCAD"), vec_pca_p1[5]);
            histos.fill(HIST("Beauty_PCAZ"), vec_pca_p1[4]);
            histos.fill(HIST("Beauty_DeltaPCAX"), delta_x);
            histos.fill(HIST("Beauty_DeltaPCAY"), delta_y);
            histos.fill(HIST("Beauty_DeltaPCAZ"), delta_z);
            if(vec_pca.at(0).at(6)!=0){
              histos.fill(HIST("Beauty_PCA_Answer"), 1);
              if(vec_pca.size()>=2){
                histos.fill(HIST("Beauty_Fake_PCAR"), vec_pca.at(1).at(0));
                histos.fill(HIST("Beauty_Fake_PCAD"), vec_pca.at(1).at(5));
                histos.fill(HIST("Beauty_Fake_PCAZ"), vec_pca.at(1).at(4));
              }
            }else{
              histos.fill(HIST("Beauty_PCA_Answer"), 0);
            }
          }else if(PartType==2){
            histos.fill(HIST("npCharm_PCAR"), vec_pca_p1[0]);
            histos.fill(HIST("npCharm_PCAD"), vec_pca_p1[5]);
            histos.fill(HIST("npCharm_PCAZ"), vec_pca_p1[4]);
            histos.fill(HIST("npCharm_DeltaPCAX"), delta_x);
            histos.fill(HIST("npCharm_DeltaPCAY"), delta_y);
            histos.fill(HIST("npCharm_DeltaPCAZ"), delta_z);
            if(vec_pca.at(0).at(6)!=0){
              histos.fill(HIST("npCharm_PCA_Answer"), 1);
              if(vec_pca.size()>=2){
                histos.fill(HIST("npCharm_Fake_PCAR"), vec_pca.at(1).at(0));
                histos.fill(HIST("npCharm_Fake_PCAD"), vec_pca.at(1).at(5));
                histos.fill(HIST("npCharm_Fake_PCAZ"), vec_pca.at(1).at(4));
              }
            }else{
              histos.fill(HIST("npCharm_PCA_Answer"), 0);
            }
          }else if(PartType==3){
            histos.fill(HIST("pCharm_PCAR"), vec_pca_p1[0]);
            histos.fill(HIST("pCharm_PCAD"), vec_pca_p1[5]);
            histos.fill(HIST("pCharm_PCAZ"), vec_pca_p1[4]);
            histos.fill(HIST("pCharm_DeltaPCAX"), delta_x);
            histos.fill(HIST("pCharm_DeltaPCAY"), delta_y);
            histos.fill(HIST("pCharm_DeltaPCAZ"), delta_z);
            if(vec_pca.at(0).at(6)!=0){
              histos.fill(HIST("pCharm_PCA_Answer"), 1);
              if(vec_pca.size()>=2){
                histos.fill(HIST("pCharm_Fake_PCAR"), vec_pca.at(1).at(0));
                histos.fill(HIST("pCharm_Fake_PCAD"), vec_pca.at(1).at(5));
                histos.fill(HIST("pCharm_Fake_PCAZ"), vec_pca.at(1).at(4));
              }
            }else{
              histos.fill(HIST("pCharm_PCA_Answer"), 0);
            }
          }
        }else if(HasPair1Track==true && HasPair2Track==true && HasPair3Track==false){
          if(PartType==1){
            histos.fill(HIST("Beauty_PCAR"), vec_pca_p1[0]);
            histos.fill(HIST("Beauty_PCAD"), vec_pca_p1[5]);
            histos.fill(HIST("Beauty_PCAZ"), vec_pca_p1[4]);
            if(vec_pca.at(0).at(6)!=0){
              histos.fill(HIST("Beauty_PCA_Answer"), 1);
            }else{
              histos.fill(HIST("Beauty_PCA_Answer"), 0);
            }
          }else if(PartType==2){
            histos.fill(HIST("npCharm_PCAR"), vec_pca_p1[0]);
            histos.fill(HIST("npCharm_PCAD"), vec_pca_p1[5]);
            histos.fill(HIST("npCharm_PCAZ"), vec_pca_p1[4]);
            if(vec_pca.at(0).at(6)!=0){
              histos.fill(HIST("npCharm_PCA_Answer"), 1);
            }else{
              histos.fill(HIST("npCharm_PCA_Answer"), 0);
            }
          }else if(PartType==3){
            histos.fill(HIST("pCharm_PCAR"), vec_pca_p1[0]);
            histos.fill(HIST("pCharm_PCAD"), vec_pca_p1[5]);
            histos.fill(HIST("pCharm_PCAZ"), vec_pca_p1[4]);
            if(vec_pca.at(0).at(6)!=0){
              histos.fill(HIST("pCharm_PCA_Answer"), 1);
            }else{
              histos.fill(HIST("pCharm_PCA_Answer"), 0);
            }
          }
        }else if(HasPair1Track==true && HasPair2Track==true && HasPair3Track==true){
          if(PartType==1){
            histos.fill(HIST("Beauty_PCAR"), vec_pca_p1[0]);
            histos.fill(HIST("Beauty_PCAD"), vec_pca_p1[5]);
            histos.fill(HIST("Beauty_PCAZ"), vec_pca_p1[4]);
            if(vec_pca.at(0).at(6)!=0){
              histos.fill(HIST("Beauty_PCA_Answer"), 1);
            }else{
              histos.fill(HIST("Beauty_PCA_Answer"), 0);
            }
          }else if(PartType==2){
            histos.fill(HIST("npCharm_PCAR"), vec_pca_p1[0]);
            histos.fill(HIST("npCharm_PCAD"), vec_pca_p1[5]);
            histos.fill(HIST("npCharm_PCAZ"), vec_pca_p1[4]);
            if(vec_pca.at(0).at(6)!=0){
              histos.fill(HIST("npCharm_PCA_Answer"), 1);
            }else{
              histos.fill(HIST("npCharm_PCA_Answer"), 0);
            }
          }else if(PartType==3){
            histos.fill(HIST("pCharm_PCAR"), vec_pca_p1[0]);
            histos.fill(HIST("pCharm_PCAD"), vec_pca_p1[5]);
            histos.fill(HIST("pCharm_PCAZ"), vec_pca_p1[4]);
            if(vec_pca.at(0).at(6)!=0){
              histos.fill(HIST("pCharm_PCA_Answer"), 1);
            }else{
              histos.fill(HIST("pCharm_PCA_Answer"), 0);
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<McCharmBeauty>(cfgc)
  };
}