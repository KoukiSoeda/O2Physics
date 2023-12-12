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

namespace o2::aod
{
namespace truepair
{
DECLARE_SOA_COLUMN(ColId, colid, int);
DECLARE_SOA_COLUMN(ColX, colx, double);
DECLARE_SOA_COLUMN(ColY, coly, double);
DECLARE_SOA_COLUMN(ColZ, colz, double);
DECLARE_SOA_COLUMN(MomPDG, mompdg, int);
DECLARE_SOA_COLUMN(MuId, muid, int);
DECLARE_SOA_COLUMN(PairId_1, pairid_1, int);
DECLARE_SOA_COLUMN(PairPDG_1, pairpdg_1, int);
DECLARE_SOA_COLUMN(PairDcaX_1, pairdcax_1, double);
DECLARE_SOA_COLUMN(PairDcaY_1, pairdcay_1, double);
DECLARE_SOA_COLUMN(PairId_2, pairid_2, int);
DECLARE_SOA_COLUMN(PairPDG_2, pairpdg_2, int);
DECLARE_SOA_COLUMN(PairDcaX_2, pairdcax_2, double);
DECLARE_SOA_COLUMN(PairDcaY_2, pairdcay_2, double);
DECLARE_SOA_COLUMN(SecVerX, secverx, double);
DECLARE_SOA_COLUMN(SecVerY, secvery, double);
DECLARE_SOA_COLUMN(SecVerZ, secverz, double);
DECLARE_SOA_COLUMN(PcaX, pcax, double);
DECLARE_SOA_COLUMN(PcaY, pcay, double);
DECLARE_SOA_COLUMN(PcaZ, pcaz, double);
DECLARE_SOA_COLUMN(PcaD, pcad, double);
DECLARE_SOA_COLUMN(PcaR, pcar, double);
DECLARE_SOA_COLUMN(MftX, x, double);
DECLARE_SOA_COLUMN(MftY, y, double);
DECLARE_SOA_COLUMN(MftZ, z, double);
DECLARE_SOA_COLUMN(Pt, pt, double);
DECLARE_SOA_COLUMN(Px, px, double);
DECLARE_SOA_COLUMN(Py, py, double);
DECLARE_SOA_COLUMN(Pz, pz, double);
DECLARE_SOA_COLUMN(P, p, double);
DECLARE_SOA_COLUMN(Chi2, chi2, double);
DECLARE_SOA_COLUMN(Phi, phi, double);
DECLARE_SOA_COLUMN(Tgl, tgl, double);
DECLARE_SOA_COLUMN(SignedPt, signed1Pt, double);
DECLARE_SOA_COLUMN(Eta, eta, double);
DECLARE_SOA_COLUMN(PairMftX, pairx, double);
DECLARE_SOA_COLUMN(PairMftY, pairy, double);
DECLARE_SOA_COLUMN(PairMftZ, pairz, double);
DECLARE_SOA_COLUMN(PairChi2, pairchi2, double);
DECLARE_SOA_COLUMN(PairPhi, pairphi, double);
DECLARE_SOA_COLUMN(PairTgl, pairtgl, double);
DECLARE_SOA_COLUMN(PairSignedPt, pairsigned1Pt, double);
DECLARE_SOA_COLUMN(ParticleType, particletype, int); // 0 - No quarks(bug)
                                                     // 1 - Beauty Decay
                                                     // 2 - Non Prompt Charm
                                                     // 3 - Prompt Charm
                                                     // 4 - LF Decay
}
DECLARE_SOA_TABLE(MCPair, "AOD", "MCPAIR",
                  truepair::ColId,
                  truepair::ColX,
                  truepair::ColY,
                  truepair::ColZ,
                  truepair::MomPDG,
                  truepair::MuId,
                  truepair::PairId_1,
                  truepair::PairPDG_1,
                  truepair::PairDcaX_1,
                  truepair::PairDcaY_1,
                  truepair::PairId_2,
                  truepair::PairPDG_2,
                  truepair::PairDcaX_2,
                  truepair::PairDcaY_2,
                  truepair::SecVerX,
                  truepair::SecVerY,
                  truepair::SecVerZ,
                  truepair::PcaX,
                  truepair::PcaY,
                  truepair::PcaZ,
                  truepair::PcaD,
                  truepair::PcaR,
                  truepair::MftX,
                  truepair::MftY,
                  truepair::MftZ,
                  truepair::Pt,
                  truepair::Px,
                  truepair::Py,
                  truepair::Pz,
                  truepair::P,
                  truepair::Chi2,
                  truepair::Phi,
                  truepair::Tgl,
                  truepair::SignedPt,
                  truepair::Eta,
                  truepair::PairMftX,
                  truepair::PairMftY,
                  truepair::PairMftZ,
                  truepair::PairChi2,
                  truepair::PairPhi,
                  truepair::PairTgl,
                  truepair::PairSignedPt,
                  truepair::ParticleType);
}

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

struct McInformation{
  Produces<aod::MCPair> mcpairtable;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<double> PCA_cut{"PCA_cut", 0.1, "Cut value for PCA"};
  Configurable<double> DCAT_min{"DCAT_min", 0.01, "Lower limit of DCAT"};
  Configurable<double> DCAT_max{"DCAT_max", 0.2, "Upper limit of DCAT"};
  Configurable<double> CosSim_cut{"CosSim_cut", 0.99, "Cut value for Cosine Simiraliyty"};
  Configurable<double> Total_dist_cut{"Total_dist_cut", 30.0, "Cut value for total distance with mu track"};
  Configurable<double> pT_cut{"pT_cut", 0.5, "Cut value for transverse momentum"};
  Configurable<double> MUON_chi2{"MUON_chi2", 1.5, "Muons chi2 in the MUON arm"};

  void init(InitContext const&){
    const AxisSpec axisDCAT{3001, -0.005, 30.005, "DCAT(cm)"};
    const AxisSpec axisPCAR{5001, -0.0005, 5.0005, "PCA(cm)"};
    const AxisSpec axisPCAD{10001, -0.0005, 10.0005, "PCAD(cm)"};
    const AxisSpec axispT{100, 0, 50, "p_{T}(GeV/c)"};
    const AxisSpec axisp{1000, 0, 100, "p(GeV/c)"};

    histos.add("Beauty_DCAT", "Beauty_DCAT", kTH1F, {axisDCAT});
    histos.add("pCharm_DCAT", "pCharm_DCAT", kTH1F, {axisDCAT});
    histos.add("npCharm_DCAT", "npCharm_DCAT", kTH1F, {axisDCAT});
    histos.add("LF_DCAT", "LF_DCAT", kTH1F, {axisDCAT});
    histos.add("Beauty_p", "Beauty_p", kTH1F, {axisp});
    histos.add("pCharm_p", "pCharm_p", kTH1F, {axisp});
    histos.add("npCharm_p", "npCharm_p", kTH1F, {axisp});
    histos.add("LF_p", "LF_p", kTH1F, {axisp});

    histos.add("Beauty_pT", "Beauty_pT", kTH1F, {axispT});
    histos.add("pCharm_pT", "pCharm_pT", kTH1F, {axispT});
    histos.add("npCharm_pT", "npCharm_pT", kTH1F, {axispT});
    histos.add("LF_pT", "LF_pT", kTH1F, {axispT});
    histos.add("Kaon_pair_pT", "Kaon_pair_pT", kTH1F, {axispT});
    histos.add("Pion_bg_pT", "Pion_bg_pT", kTH1F, {axispT});

    histos.add("Beauty_PCAR", "Beauty_PCAR", kTH1F, {axisPCAR});
    histos.add("npCharm_PCAR", "npCharm_PCAR", kTH1F, {axisPCAR});
    histos.add("pCharm_PCAR", "pCharm_PCAR", kTH1F, {axisPCAR});
    histos.add("LF_PCAR", "LF_PCAR", kTH1F, {axisPCAR});
    histos.add("Beauty_PCAD", "Beauty_PCAD", kTH1F, {axisPCAD});
    histos.add("npCharm_PCAD", "npCharm_PCAD", kTH1F, {axisPCAD});
    histos.add("pCharm_PCAD", "pCharm_PCAD", kTH1F, {axisPCAD});
    histos.add("LF_PCAD", "LF_PCAD", kTH1F, {axisPCAD});
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
      int mom_pdg = mumom.pdgCode();
      auto Daughters = mumom.daughters_as<aod::McParticles>();
      double muon_p = fwdtrack.p();

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
        histos.fill(HIST("Beauty_pT"), fwdtrack.pt());
        histos.fill(HIST("Beauty_p"), fwdtrack.p());
      }else if(HasLightParent==false && HasBeautyParent==true && HasCharmParent==true){
        PartType = 2; // Non Prompt Charm Muon
        histos.fill(HIST("npCharm_pT"), fwdtrack.pt());
        histos.fill(HIST("npCharm_p"), fwdtrack.p());
      }else if(HasLightParent==false && HasBeautyParent==false && HasCharmParent==true){
        PartType = 3; // Prompt Charm Muon
        histos.fill(HIST("pCharm_pT"), fwdtrack.pt());
        histos.fill(HIST("pCharm_p"), fwdtrack.p());
      }else if(HasLightParent==true && IsSecondary==false){
        PartType = 4; // LF Muon(BG)
        histos.fill(HIST("LF_pT"), fwdtrack.pt());
        histos.fill(HIST("LF_p"), fwdtrack.p());
      }

      double secverx, secvery, secverz;
      vector<int> vec_daughter;
      vec_daughter.clear();
      bool HasPair = false;
      for(auto& Daughter : Daughters){
        int d_pdg = fabs(Daughter.pdgCode());
        if(d_pdg==13 && Daughter.globalIndex()==muid){
          secverx = Daughter.vx();
          secvery = Daughter.vy();
          secverz = Daughter.vz();
        }
        vec_daughter.push_back(Daughter.globalIndex());
      }
      if(vec_daughter.size()>1) HasPair = true;
      if(!HasPair) continue;

      bool HasMuonTrack=false, HasPair1Track=false, HasPair2Track=false, HasPair3Track=false;
      double mudcaX, mudcaY, mumftX, mumftY, mumftZ, muchi2, muphi, mutgl, musignedpt, mupt, mup, mupx, mupy, mupz,mueta;
      double p1dcaX, p1dcaY, p1mftX, p1mftY, p1mftZ, p1chi2, p1phi, p1tgl, p1signedpt;
      double p2dcaX, p2dcaY, p2mftX, p2mftY, p2mftZ, p2chi2, p2phi, p2tgl, p2signedpt;
      double p3dcaX, p3dcaY, p3mftX, p3mftY, p3mftZ, p3chi2, p3phi, p3tgl, p3signedpt;
      int mcColId_mu, mcColId_p1, mcColId_p2;
      for(auto& mfttrack : mfttracks){
        if(!mfttrack.has_mcParticle()) continue;
        if(mfttrack.mcParticle().mcCollisionId()==mcCol_id){
          auto mcParticle_mft = mfttrack.mcParticle();
          if(mcParticle_mft.globalIndex()==muid){
            double mftchi2 = mfttrack.chi2();
            SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
            vector<double> mftv2;
            SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
            o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs2, mftchi2};
            mftpars1.propagateToZlinear(Col_z);

            mudcaX = mftpars1.getX();
            mudcaY = mftpars1.getY();
            mumftX = mfttrack.x();
            mumftY = mfttrack.y();
            mumftZ = mfttrack.z();
            muchi2 = mfttrack.chi2();
            muphi = mfttrack.phi();
            mutgl = mfttrack.tgl();
            mueta = mfttrack.eta();
            musignedpt = mfttrack.signed1Pt();
            mupt = mfttrack.pt();
            mup = mfttrack.p();
            mupx = mfttrack.px();
            mupy = mfttrack.py();
            mupz = mfttrack.pz();
            mcColId_mu = mcParticle_mft.mcCollisionId();
            double mu_dcat = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
            
            HasMuonTrack = true;
            if(PartType==1){
              histos.fill(HIST("Beauty_DCAT"), mu_dcat);
            }else if(PartType==2){
              histos.fill(HIST("npCharm_DCAT"), mu_dcat);
            }else if(PartType==3){
              histos.fill(HIST("pCharm_DCAT"), mu_dcat);
            }else if(PartType==4){
              histos.fill(HIST("LF_DCAT"), mu_dcat);
            }
            
          }
          if(find(vec_daughter.begin(), vec_daughter.end(), mcParticle_mft.globalIndex())!=vec_daughter.end() && HasPair1Track==false){
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
            p1phi = mfttrack.phi();
            p1tgl = mfttrack.tgl();
            p1signedpt = mfttrack.signed1Pt();
            mcColId_p1 = mcParticle_mft.mcCollisionId();
            histos.fill(HIST("Kaon_pair_pT"), mfttrack.pt());
            HasPair1Track = true;
          }else if(find(vec_daughter.begin(), vec_daughter.end(), mcParticle_mft.globalIndex())!=vec_daughter.end() && HasPair1Track==true && HasPair2Track==false){
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
            p2phi = mfttrack.phi();
            p2tgl = mfttrack.tgl();
            p2signedpt = mfttrack.signed1Pt();
            mcColId_p2 = mcParticle_mft.mcCollisionId();
            HasPair2Track = true;
          }else if(find(vec_daughter.begin(), vec_daughter.end(), mcParticle_mft.globalIndex())!=vec_daughter.end() && HasPair1Track==true && HasPair2Track==true && HasPair3Track==false){
            double mftchi2 = mfttrack.chi2();
            SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
            vector<double> mftv2;
            SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
            o2::track::TrackParCovFwd mftpars3{mfttrack.z(), mftpars, mftcovs2, mftchi2};
            mftpars3.propagateToZlinear(Col_z);

            p3dcaX = mftpars3.getX();
            p3dcaY = mftpars3.getY();
            p3mftX = mfttrack.x();
            p3mftY = mfttrack.y();
            p3mftZ = mfttrack.z();
            p3chi2 = mfttrack.chi2();
            p3phi = mfttrack.phi();
            p3tgl = mfttrack.tgl();
            p3signedpt = mfttrack.signed1Pt();
            mcColId_p3 = mcParticle_mft.mcCollisionId();
            HasPair3Track = true;
          }else{
            histos.fill(HIST("Pion_bg_pT"), mfttrack.pt());
          }
        }
      }
      if(HasMuonTrack && HasPair1Track && mcColId_mu==mcColId_p1){
        vector<double> vec_mftdca;
        vec_mftdca.clear();
        vec_mftdca = {mumftX,mumftY,mumftZ,mudcaX,mudcaY,p1mftX,p1mftY,p1mftZ,p1dcaX,p1dcaY,Col_x,Col_y,Col_z};
        auto vec_pca = PCA_Cal(vec_mftdca);

        if(PartType==1){
          histos.fill(HIST("Beauty_PCAR"), vec_pca[0]);
          histos.fill(HIST("Beauty_PCAD"), vec_pca[5]);
        }else if(PartType==2){
          histos.fill(HIST("npCharm_PCAR"), vec_pca[0]);
          histos.fill(HIST("npCharm_PCAD"), vec_pca[5]);
        }else if(PartType==3){
          histos.fill(HIST("pCharm_PCAR"), vec_pca[0]);
          histos.fill(HIST("pCharm_PCAD"), vec_pca[5]);
        }else if(PartType==4){
          histos.fill(HIST("LF_PCAR"), vec_pca[0]);
          histos.fill(HIST("LF_PCAD"), vec_pca[5]);
        }
        mcpairtable(mcColId_mu,Col_x,Col_y,Col_z,mom_pdg,muid,pid_1,pdg_1,p1dcaX,p1dcaY,pid_2,pdg_2,p2dcaX,p2dcaY,secverx,secvery,secverz,vec_pca[2],vec_pca[3],vec_pca[4],vec_pca[5],vec_pca[0],mumftX,mumftY,mumftZ,mupt,mupx,mupy,mupz,mup,muchi2,muphi,mutgl,musignedpt,mueta,p2mftX,p2mftY,p2mftZ,p2chi2,p2phi,p2tgl,p2signedpt,PartType);
      }
    }
  }
};

struct MyAnalysisTask{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<double> PCA_cut{"PCA_cut", 0.1, "Cut value for PCA"};
  Configurable<double> DCAT_min{"DCAT_min", 0.01, "Lower limit of DCAT"};
  Configurable<double> DCAT_max{"DCAT_max", 0.2, "Upper limit of DCAT"};
  Configurable<double> CosSim_cut{"CosSim_cut", 0.99, "Cut value for Cosine Simiraliyty"};
  Configurable<double> Total_dist_cut{"Total_dist_cut", 30.0, "Cut value for total distance with mu track"};
  Configurable<double> PCAD_cut_max{"PCAD_cut_max", 1, "Max cut value for PCAD"};
  Configurable<double> PCAD_cut_min{"PCAD_cut_min", 0.15, "Min cut value for PCAD"};
  Configurable<double> pT_cut{"pT_cut", 0.5, "Cut value for transverse momentum"};
  Configurable<double> pDCAT_cut{"pDCAT_cut", 12000, "The unit of pDCAT is cm*GeV/c"};
  Configurable<double> MUON_chi2{"MUON_chi2", 1.5, "Muons chi2 in the MUON arm"};

  void init(InitContext const&){
    const AxisSpec axisDCAT{3001, -0.005, 30.005, "DCAT(cm)"};
    const AxisSpec axispT{100, 0, 50, "p_{T}(GeV/c)"};
    const AxisSpec axisnEvent{100, 0, 1, "nGoodEvents/nEvents"};
    const AxisSpec axisMass{1000, 0, 10, "InvMass(GeV/c^2)"};

    histos.add("AllMuon_pT", "AllMuon_pT", kTH1F, {axispT});
    histos.add("Beauty_pT", "Beauty_pT", kTH1F, {axispT});
    histos.add("pCharm_pT", "pCharm_pT", kTH1F, {axispT});
    histos.add("npCharm_pT", "npCharm_pT", kTH1F, {axispT});
    histos.add("BG_pT", "BG_pT", kTH1F, {axispT});

    histos.add("pCharm_DCATcut", "pCharm_DCATcut", kTH1F, {axispT});
    histos.add("BG_DCATcut", "BG_DCATcut", kTH1F, {axispT});
    histos.add("pCharm_DCAT_PCAcut", "pCharm_DCAT_PCAcut", kTH1F, {axispT});
    histos.add("BG_DCAT_PCAcut", "BG_DCAT_PCAcut", kTH1F, {axispT});

    histos.add("AllEvent_GoodEvent_forPCAcut_Signal", "AllEvent_GoodEvent_forPCAcut_Signal", kTH1F, {axisnEvent});
    histos.add("AllEvent_GoodEvent_forPCAcut_BG", "AllEvent_GoodEvent_forPCAcut_BG", kTH1F, {axisnEvent});
    histos.add("DecayMuon_InvMass", "DecayMuon_InvMass", kTH1F, {axisMass});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions,
							 soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
							 soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
							 aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
							 aod::McCollisions const&,
               aod::MCPair const& truepairs)
  {
    const double k_mass2 = pow(0.493677,2);
    const double mu_mass2 = pow(0.105658,2);
    if(truepairs.size()!=0){
      for(auto& truepair : truepairs){
        int colid = -1;
        double Col_x, Col_y, Col_z;
        bool Has_MCCol = false;
        for(auto& t_col : collisions){
          if(!t_col.has_mcCollision()) continue;
          if(truepair.colid()==t_col.mcCollisionId()){
            Col_x = t_col.posX();
            Col_y = t_col.posY();
            Col_z = t_col.posZ();
            colid = t_col.mcCollisionId();
            Has_MCCol = true;
            break;
          }
        }
        if(fabs(Col_z)>10 || Has_MCCol==false) continue;
        for(auto& fwdtrack : fwdtracks){
          if(!fwdtrack.has_mcParticle()) continue;
          if(fwdtrack.mcParticle().mcCollisionId()!=colid) continue;
          if(fwdtrack.trackType()==1 || fwdtrack.trackType()==3 || fwdtrack.trackType()==4) continue;
          if(fwdtrack.eta()<=(-3.6) || fwdtrack.eta()>=(-2.5)) continue;
          if(fwdtrack.pt()<pT_cut) continue;
          auto mcParticle_fwd = fwdtrack.mcParticle();
          if(fabs(mcParticle_fwd.pdgCode())!=13) continue;

          double mudcaX,mudcaY,mumftX,mumftY,mumftZ,mueta,muphi;
          double fwdp, fwdpt, fwdpx, fwdpy, fwdpz;
          int PartType = 0;
          bool Signal_Muon = false;
          bool BG_NonPair_Muon = false;
          int trackId = -1;
          SMatrix5 mupars(truepair.x(), truepair.y(), truepair.phi(), truepair.tgl(), truepair.signed1Pt());
          vector<double> muv2;
          SMatrix55 mucovs2(muv2.begin(), muv2.end());
          o2::track::TrackParCovFwd mupars2{truepair.z(), mupars, mucovs2, truepair.chi2()};

          if(truepair.muid()==mcParticle_fwd.globalIndex()){
            Signal_Muon = true;
            trackId = truepair.muid();
            mupars2.propagateToZlinear(Col_z);
            mudcaX = mupars2.getX();
            mudcaY = mupars2.getY();
            mumftX = truepair.x();
            mumftY = truepair.y();
            mumftZ = truepair.z();
            mueta = truepair.eta();
            muphi = truepair.phi();
            fwdp = truepair.p();
            fwdpt = truepair.pt();
            fwdpx = truepair.px();
            fwdpy = truepair.py();
            fwdpz = truepair.pz();
            PartType = truepair.particletype();

          }else{
            trackId = fwdtrack.matchMFTTrackId();
            SMatrix5 mupars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
            vector<double> muv2;
            SMatrix55 mucovs2(muv2.begin(), muv2.end());
            o2::track::TrackParCovFwd mupars2{fwdtrack.z(), mupars, mucovs2, fwdtrack.chi2()};
            mupars2.propagateToZlinear(Col_z);
            mudcaX = mupars2.getX();
            mudcaY = mupars2.getY();
            mumftX = fwdtrack.x();
            mumftY = fwdtrack.y();
            mumftZ = fwdtrack.z();
            mueta = fwdtrack.eta();
            muphi = fwdtrack.phi();
            fwdp = fwdtrack.p();
            fwdpt = fwdtrack.pt();
            fwdpx = fwdtrack.px();
            fwdpy = fwdtrack.py();
            fwdpz = fwdtrack.pz();
            double fwd_pDcat = fwdp*sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
            if(fwd_pDcat>pDCAT_cut) continue;

            BG_NonPair_Muon = true;
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
          }

          histos.fill(HIST("AllMuon_pT"), fwdpt);
          double dcat_mu = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
          bool Good_muon = false;
          if(PartType==1){
            histos.fill(HIST("Beauty_pT"), fwdpt);
          }else if(PartType==2){
            histos.fill(HIST("npCharm_pT"), fwdpt);
          }else if(PartType==3){
            histos.fill(HIST("pCharm_pT"), fwdpt);
            if(dcat_mu>DCAT_min && dcat_mu<DCAT_max){
              histos.fill(HIST("pCharm_DCATcut"), fwdpt);
              Good_muon = true;
            }
          }else if(PartType==4){
            histos.fill(HIST("BG_pT"), fwdpt);
            if(dcat_mu>DCAT_min && dcat_mu<DCAT_max){
              histos.fill(HIST("BG_DCATcut"), fwdpt);
              Good_muon = true;
            }
          }

          if(Good_muon==false) continue;
          //if(PartType==4) cout << "BG_";
          vector<vector<double>> vec_pca;
          for(auto& mfttrack : mfttracks){
            if(!mfttrack.has_mcParticle()) continue;
            if(mfttrack.mcParticle().mcCollisionId()!=colid) continue;
            if(Signal_Muon==true && mfttrack.mcParticleId()==trackId) continue;
            if(BG_NonPair_Muon==true && mfttrack.globalIndex()==trackId) continue;
            auto mcParticle_mft = mfttrack.mcParticle();
            if(fabs(mcParticle_mft.pdgCode())==2212) continue;
            double mftchi2 = mfttrack.chi2();

            if(mftchi2>2.5) continue;
            SMatrix5 mftpars_2(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
            vector<double> mftv2;
            SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
            o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars_2, mftcovs2, mftchi2};
            mftpars2.propagateToZlinear(Col_z);
            double candcaX = mftpars2.getX();
            double candcaY = mftpars2.getY();
            double canmftX = mfttrack.x();
            double canmftY = mfttrack.y();
            double canmftZ = mfttrack.z();
            double pair_id = mfttrack.mcParticleId();

            vector<double> vec_mftdca;
            vec_mftdca.clear();
            vec_mftdca = {mumftX,mumftY,mumftZ,mudcaX,mudcaY,canmftX,canmftY,canmftZ,candcaX,candcaY,Col_x,Col_y,Col_z};
            auto vec_pca_par = PCA_Cal(vec_mftdca);

            double add_dist = 0;
            for(int i=0; i<5000; i++){
              auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
              if((cut_z-Col_z)>0) break;
              mupars2.propagateToZlinear(cut_z);
              mftpars2.propagateToZlinear(cut_z);
              auto diff = sqrt(pow(mftpars2.getX()-mupars2.getX(),2)+pow(mftpars2.getY()-mupars2.getY(),2));
              add_dist += diff;
            }

            double ans = 0;
            if(pair_id==truepair.pairid_1()){
              ans = 1;
            }else if(pair_id==truepair.pairid_2()){
              ans = 2;
            }
            vec_pca_par.push_back(add_dist);
            vec_pca_par.push_back(ans);
            //vec_pca = {r_xyz, cosxy, pcaX, pcaY, pcaZ, pcaD, add_dist, ans}
            vec_pca.push_back(vec_pca_par);
          }
          if(vec_pca.empty()) continue;
          sort(vec_pca.begin(), vec_pca.end());
          int n_goodevents = 0;
          bool HasGoodMuon_Charm = false;
          bool HasGoodMuon_BG = false;
          //cout << "Data = ";
          for(int i=0; i<vec_pca.size(); i++){
            //cout << "{ " << vec_pca.at(i).at(0) << ", " <<  vec_pca.at(i).at(1) << ", " << vec_pca.at(i).at(2) << ", " << vec_pca.at(i).at(3) << ", " << vec_pca.at(i).at(4) << ", " << vec_pca.at(i).at(5) << ", " << vec_pca.at(i).at(6) << ", " << vec_pca.at(i).at(7);
            if(vec_pca.at(i).at(0)<PCA_cut){
              if(vec_pca.at(i).at(1)>CosSim_cut){
                if(vec_pca.at(i).at(6)<Total_dist_cut){
                  if(vec_pca.at(i).at(5)>PCAD_cut_min &&vec_pca.at(i).at(5)<PCAD_cut_max){
                    if(PartType==3){
                      HasGoodMuon_Charm = true;
                    }else if(PartType==4){
                      HasGoodMuon_BG = true;
                    }
                    //cout << " Good Event!";
                    n_goodevents++;
                  }
                }
              }
            }
            //cout << " }" << endl;
          }
          double d_n_goodevents = n_goodevents;
          double d_n_vecsize = vec_pca.size();
          if(PartType==3){            
            if(HasGoodMuon_Charm==true){
              histos.fill(HIST("pCharm_DCAT_PCAcut"), fwdpt);
              histos.fill(HIST("AllEvent_GoodEvent_forPCAcut_Signal"), d_n_goodevents/d_n_vecsize);
            }
          }else if(PartType==4){            
            if(HasGoodMuon_BG==true){
              histos.fill(HIST("BG_DCAT_PCAcut"), fwdpt);
              histos.fill(HIST("AllEvent_GoodEvent_forPCAcut_BG"), d_n_goodevents/d_n_vecsize);
            }
          }
        }
      }
    }else if(truepairs.size()==0){
      for(auto& fwdtrack : fwdtracks){
        if(!fwdtrack.has_mcParticle()) continue;
        if(fwdtrack.eta()<=(-3.6) || fwdtrack.eta()>=(-2.5)) continue;
        if(fwdtrack.trackType()==1 || fwdtrack.trackType()==3 || fwdtrack.trackType()==4) continue;
        if(fwdtrack.pt()<0.5) continue;
        auto mcParticle_fwd = fwdtrack.mcParticle();
        if(fabs(mcParticle_fwd.pdgCode())!=13) continue;
        double Col_x, Col_y, Col_z;
        int colid;
        bool HasCol = false;
        for(auto& t_col : collisions){
          if(!t_col.has_mcCollision()) continue;
          if(t_col.mcCollisionId()==mcParticle_fwd.mcCollisionId()){
            Col_x = t_col.posX();
            Col_y = t_col.posY();
            Col_z = t_col.posZ();
            colid = t_col.mcCollisionId();
            HasCol = true;
            break;
          }
        }
        if(fabs(Col_z)>10 || HasCol==false) continue;
        SMatrix5 mupars_bg(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
        vector<double> muv2;
        SMatrix55 mucovs2(muv2.begin(), muv2.end());
        o2::track::TrackParCovFwd mupars2{fwdtrack.z(), mupars_bg, mucovs2, fwdtrack.chi2()};
        mupars2.propagateToZlinear(Col_z);
        int trackId = fwdtrack.matchMFTTrackId();
        double bg_mudcaX = mupars2.getX();
        double bg_mudcaY = mupars2.getY();
        double bg_mumftX = fwdtrack.x();
        double bg_mumftY = fwdtrack.y();
        double bg_mumftZ = fwdtrack.z();
        double bg_pt = fwdtrack.pt();
        double fwdp = fwdtrack.p();
        double fwdpx = fwdtrack.px();
        double fwdpy = fwdtrack.py();
        double fwdpz = fwdtrack.pz();
        double fwd_pDcat = fwdtrack.p()*sqrt(pow(bg_mudcaX-Col_x,2)+pow(bg_mudcaY-Col_y,2));
        if(fwd_pDcat>pDCAT_cut) continue;

        bool IsSecondary = false;
        bool HasLightParent = false;
        bool HasCharmParent = false;
        bool HasBeautyParent = false;
        int PartType = 0;
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
        histos.fill(HIST("AllMuon_pT"), bg_pt);
        double dcat_mu = sqrt(pow(bg_mudcaX-Col_x,2)+pow(bg_mudcaY-Col_y,2));
        bool Good_muon = false;
        if(PartType==1){
          histos.fill(HIST("Beauty_pT"), bg_pt);
        }else if(PartType==2){
          histos.fill(HIST("npCharm_pT"), bg_pt);
        }else if(PartType==3){
          histos.fill(HIST("pCharm_pT"), bg_pt);
          if(dcat_mu>DCAT_min && dcat_mu<DCAT_max){
            histos.fill(HIST("pCharm_DCATcut"), bg_pt);
            Good_muon = true;
          }
        }else if(PartType==4){
          histos.fill(HIST("BG_pT"), bg_pt);
          if(dcat_mu>DCAT_min && dcat_mu<DCAT_max){
            histos.fill(HIST("BG_DCATcut"), bg_pt);
            Good_muon = true;
          }
        }

        if(Good_muon==false) continue;
        vector<vector<double>> vec_pca;
        for(auto& mfttrack : mfttracks){
          if(!mfttrack.has_mcParticle()) continue;
          if(mfttrack.mcParticle().mcCollisionId()!=colid) continue;
          if(mfttrack.globalIndex()==trackId) continue;
          auto mcParticle_mft = mfttrack.mcParticle();
          if(fabs(mcParticle_mft.pdgCode())==2212) continue;
          double mftchi2 = mfttrack.chi2();
          if(PartType==4){
            double bg_invmass2 = pow(sqrt(mu_mass2+pow(fwdp,2))+sqrt(k_mass2+pow(mfttrack.p(),2)),2) - (pow(fwdpx+mfttrack.px(),2)+pow(fwdpy+mfttrack.py(),2)+pow(fwdpz+mfttrack.pz(),2));
            double bg_invmass = sqrt(bg_invmass2);
            histos.fill(HIST("DecayMuon_InvMass"), bg_invmass);
          }
          if(mftchi2>2.5) continue;
          SMatrix5 mftpars_2(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
          vector<double> mftv2;
          SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
          o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars_2, mftcovs2, mftchi2};
          mftpars2.propagateToZlinear(Col_z);

          double candcaX = mftpars2.getX();
          double candcaY = mftpars2.getY();
          double canmftX = mfttrack.x();
          double canmftY = mfttrack.y();
          double canmftZ = mfttrack.z();

          vector<double> vec_mftdca;
          vec_mftdca.clear();
          vec_mftdca = {bg_mumftX,bg_mumftY,bg_mumftZ,bg_mudcaX,bg_mudcaY,canmftX,canmftY,canmftZ,candcaX,candcaY,Col_x,Col_y,Col_z};
          auto vec_pca_par = PCA_Cal(vec_mftdca);

          double add_dist = 0;
          for(int i=0; i<5000; i++){
            auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
            if((cut_z-Col_z)>0) break;
            mupars2.propagateToZlinear(cut_z);
            mftpars2.propagateToZlinear(cut_z);
            auto diff = sqrt(pow(mftpars2.getX()-mupars2.getX(),2)+pow(mftpars2.getY()-mupars2.getY(),2));
            add_dist += diff;
          }

          double ans = 0;
          
          vec_pca_par.push_back(add_dist);
          vec_pca_par.push_back(ans);
          //vec_pca_par = {r_xyz, cosxy, pcaX, pcaY, pcaZ, pcaD, add_dist, ans}
          vec_pca.push_back(vec_pca_par);
        }
        if(vec_pca.empty()) continue;
        sort(vec_pca.begin(), vec_pca.end());
        int n_goodevents = 0;
        bool HasGoodMuon_Charm = false;
        bool HasGoodMuon_BG = false;
        //if(PartType==3) cout << "NonPairCharm_";
        //if(PartType==4) cout << "BG_";
        //cout << "Data = ";
        for(int i=0; i<vec_pca.size(); i++){
          //cout << "{ " << vec_pca.at(i).at(0) << ", " <<  vec_pca.at(i).at(1) << ", " << vec_pca.at(i).at(2) << ", " << vec_pca.at(i).at(3) << ", " << vec_pca.at(i).at(4) << ", " << vec_pca.at(i).at(5) << ", " << vec_pca.at(i).at(6) << ", " << vec_pca.at(i).at(7);
          if(vec_pca.at(i).at(0)<PCA_cut){
            if(vec_pca.at(i).at(1)>CosSim_cut){
              if(vec_pca.at(i).at(6)<Total_dist_cut){
                if(vec_pca.at(i).at(5)>PCAD_cut_min &&vec_pca.at(i).at(5)<PCAD_cut_max){
                  if(PartType==3){
                    HasGoodMuon_Charm =true;
                  }else if(PartType==4){
                    HasGoodMuon_BG = true;
                  }
                  //cout << " Good Event!";
                  n_goodevents++;
                }
              }
            }
          }
          //cout << " }" << endl;
        }
        double d_n_goodevents = n_goodevents;
        double d_n_vecsize = vec_pca.size();
        if(PartType==3){
          if(HasGoodMuon_Charm==true){
            histos.fill(HIST("pCharm_DCAT_PCAcut"), bg_pt);
            histos.fill(HIST("AllEvent_GoodEvent_forPCAcut_Signal"), d_n_goodevents/d_n_vecsize);
          }
        }else if(PartType==4){
          if(HasGoodMuon_BG==true){
            histos.fill(HIST("BG_DCAT_PCAcut"), bg_pt);
            histos.fill(HIST("AllEvent_GoodEvent_forPCAcut_BG"), d_n_goodevents/d_n_vecsize);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<McInformation>(cfgc),
    adaptAnalysisTask<MyAnalysisTask>(cfgc)
  };
}