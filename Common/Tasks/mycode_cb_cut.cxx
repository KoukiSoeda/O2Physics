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
struct McInformation{
  Produces<aod::MCPair> mcpairtable;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&){
    const AxisSpec axisDCAT{3001, -0.005, 30.005, "cm"};
    const AxisSpec axisPCAR{5001, -0.0005, 5.0005, "cm"};
    const AxisSpec axisPCAD{10001, -0.0005, 10.0005, "cm"};

    histos.add("Beauty_DCAT", "Beauty_DCAT", kTH1F, {axisDCAT});
    histos.add("pCharm_DCAT", "pCharm_DCAT", kTH1F, {axisDCAT});
    histos.add("npCharm_DCAT", "npCharm_DCAT", kTH1F, {axisDCAT});
    histos.add("BG_inMFT_DCAT", "BG_inMFT_DCAT", kTH1F, {axisDCAT});

    histos.add("D0_Distance", "D0_Distance", kTH1F, {axisPCAD});
    histos.add("Beauty_PCAR_Dist", "Beauty_PCAR_Dist", kTH1F, {axisPCAR});
    histos.add("npCharm_PCAR_Dist", "npCharm_PCAR_Dist", kTH1F, {axisPCAR});
    histos.add("pCharm_PCAR_Dist", "pCharm_PCAR_Dist", kTH1F, {axisPCAR});
    histos.add("BG_PCAR_Dist", "BG_PCAR_Dist", kTH1F, {axisPCAR});
    histos.add("Beauty_PCAD_Dist", "Beauty_PCAD_Dist", kTH1F, {axisPCAD});
    histos.add("npCharm_PCAD_Dist", "npCharm_PCAD_Dist", kTH1F, {axisPCAD});
    histos.add("pCharm_PCAD_Dist", "pCharm_PCAD_Dist", kTH1F, {axisPCAD});
    histos.add("BG_PCAD_Dist", "BG_PCAD_Dist", kTH1F, {axisPCAD});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions,
               soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA, aod::FwdTracksCov> const& fwdtracks, 
               soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
               aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
               aod::McCollisions const&,
               ExtBCs const& ebcs)
  {
    // Create AmbTrackTable
    vector<uint64_t> ambTrackIds;
    ambTrackIds.clear();
    for(auto& amfttrack : amfttracks){
      ambTrackIds.push_back(amfttrack.mfttrackId());
    }

    // Muon Source and N of daughters
    for(auto& fwdtrack : fwdtracks){
      if(!fwdtrack.has_collision() || !fwdtrack.has_mcParticle()) continue;
      auto mcParticle_fwd = fwdtrack.mcParticle();
      if(fwdtrack.trackType()==1 || fwdtrack.trackType()==3 || fwdtrack.trackType()==4 ) continue;
      if(fwdtrack.eta()<=(-3.6) || fwdtrack.eta()>=(-2.5)) continue;

      int mcCol_id;
      double Col_x, Col_y, Col_z;
      for(auto& t_col : collisions){
        if(!t_col.has_mcCollision()) continue;
        if(mcParticle_fwd.mcCollisionId()==t_col.mcCollisionId()){
          mcCol_id = t_col.mcCollisionId();
          Col_x = t_col.posX();
          Col_y = t_col.posY();
          Col_z = t_col.posZ();
        }
      }
      if(fabs(Col_z)>10) continue;
      
      if(fabs(mcParticle_fwd.pdgCode())!=13) continue;
      int muid = mcParticle_fwd.globalIndex();
      auto mumom = mcParticle_fwd.mothers_first_as<aod::McParticles>();
      int mom_pdg = mumom.pdgCode();
      auto Daughters = mumom.daughters_as<aod::McParticles>();
      
      // Muon Source
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

      // N of daughters
      double secverx, secvery, secverz;
      int pid_1=0, pid_2=0, pdg_1=0, pdg_2=0;
      bool HasPair = false;
      for(auto& Daughter : Daughters){
        int d_pdg = fabs(Daughter.pdgCode());
        if(d_pdg==13 && Daughter.globalIndex()==muid){
          secverx = Daughter.vx();
          secvery = Daughter.vy();
          secverz = Daughter.vz();
        }else if(d_pdg!=13 && d_pdg>20 && pid_1==0){
          pid_1 = Daughter.globalIndex();
          pdg_1 = Daughter.pdgCode();
          HasPair = true;
        }else if(d_pdg!=13 && d_pdg>20 && pid_1!=0){
          pid_2 = Daughter.globalIndex();
          pdg_2 = Daughter.pdgCode();
        }
      }

      if(!HasPair) continue;

      bool HasMuonTrack=false, HasPair1Track=false, HasPair2Track=false;
      double mudcaX, mudcaY, mumftX, mumftY, mumftZ, muchi2, muphi, mutgl, musignedpt, mupt, mup, mupx, mupy, mupz,mueta;
      double p1dcaX, p1dcaY, p1mftX, p1mftY, p1mftZ, p1chi2, p1phi, p1tgl, p1signedpt;
      double p2dcaX, p2dcaY, p2mftX, p2mftY, p2mftZ, p2chi2, p2phi, p2tgl, p2signedpt;
      int mcColId_mu, mcColId_p1, mcColId_p2;
      for(auto& mfttrack : mfttracks){
        if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
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
            if(PartType==1){
              histos.fill(HIST("Beauty_DCAT"), mu_dcat);
            }else if(PartType==2){
              histos.fill(HIST("npCharm_DCAT"), mu_dcat);
            }else if(PartType==3){
              histos.fill(HIST("pCharm_DCAT"), mu_dcat);
            }else if(PartType==4){
              histos.fill(HIST("BG_inMFT_DCAT"), mu_dcat);
            }
            HasMuonTrack = true;
          }
          if(mcParticle_mft.globalIndex()==pid_1){
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
            HasPair1Track = true;
          }
          if(mcParticle_mft.globalIndex()==pid_2){
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
          }
        }
      }

      if(HasMuonTrack && HasPair1Track && mcColId_mu==mcColId_p1){
        auto unit_Na = sqrt(pow(mumftX-mudcaX,2)+pow(mumftY-mudcaY,2)+pow(mumftZ-Col_z,2));
        auto unit_Nc = sqrt(pow(p1mftX-p1dcaX,2)+pow(p1mftY-p1dcaY,2)+pow(p1mftZ-Col_z,2));
        auto Nax = (mumftX-mudcaX)/unit_Na;
        auto Nay = (mumftY-mudcaY)/unit_Na;
        auto Naz = (mumftZ-Col_z)/unit_Na;
        auto Ncx = (p1mftX-p1dcaX)/unit_Nc;
        auto Ncy = (p1mftY-p1dcaY)/unit_Nc;
        auto Ncz = (p1mftZ-Col_z)/unit_Nc;
        auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
        auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
        auto A3 = (mudcaX-p1dcaX)*Nax + (mudcaY-p1dcaY)*Nay + (Col_z-Col_z)*Naz;
        auto B1 = A2;
        auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
        auto B3 = (p1dcaX-mudcaX)*Ncx + (p1dcaY-mudcaY)*Ncy + (Col_z-Col_z)*Ncz;
        auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
        auto s = -((A2*t+A3)/A1);
        double predict_mux = mudcaX + s*Nax;
        double predict_muy = mudcaY + s*Nay;
        double predict_muz = Col_z + s*Naz;
        double predict_p1x = p1dcaX + t*Ncx;
        double predict_p1y = p1dcaY + t*Ncy;
        double predict_p1z = Col_z + t*Ncz;

        double predict_r = sqrt(pow(predict_mux-predict_p1x,2)+pow(predict_muy-predict_p1y,2)+pow(predict_muz-predict_p1z,2));
        double pcaX = ((predict_mux+predict_p1x)/2)-Col_x;
        double pcaY = ((predict_muy+predict_p1y)/2)-Col_y;
        double pcaZ = ((predict_muz+predict_p1z)/2)-Col_z;
        double pcaD = sqrt(pow(pcaX,2)+pow(pcaY,2)+pow(pcaZ,2));

        if(PartType==1){
          histos.fill(HIST("Beauty_PCAR_Dist"), predict_r);
          histos.fill(HIST("Beauty_PCAD_Dist"), pcaD);
        }else if(PartType==2){
          histos.fill(HIST("npCharm_PCAR_Dist"), predict_r);
          histos.fill(HIST("npCharm_PCAD_Dist"), pcaD);
        }else if(PartType==3){
          histos.fill(HIST("pCharm_PCAR_Dist"), predict_r);
          histos.fill(HIST("pCharm_PCAD_Dist"), pcaD);
        }else if(PartType==4){
          histos.fill(HIST("BG_PCAR_Dist"), predict_r);
          histos.fill(HIST("BG_PCAD_Dist"), pcaD);
        }

        if(fabs(mom_pdg)==421) histos.fill(HIST("D0_Distance"), sqrt(pow(secverx-Col_x,2)+pow(secvery-Col_y,2)+pow(secverz-Col_z,2)));

        mcpairtable(mcColId_mu,Col_x,Col_y,Col_z,mom_pdg,muid,pid_1,pdg_1,p1dcaX,p1dcaY,pid_2,pdg_2,p2dcaX,p2dcaY,secverx,secvery,secverz,pcaX,pcaY,pcaZ,pcaD,predict_r,mumftX,mumftY,mumftZ,mupt,mupx,mupy,mupz,mup,muchi2,muphi,mutgl,musignedpt,mueta,p2mftX,p2mftY,p2mftZ,p2chi2,p2phi,p2tgl,p2signedpt,PartType);
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

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  int runNumber = -1;
  float Bz = 0; 
  float Bz_global = 0;                                        // Magnetic field for MFT
  static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
  static constexpr double centerGlobal[3] = {0, 0, -25};
  o2::parameters::GRPMagField* grpmag = nullptr;

  void init(InitContext const&){
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    const AxisSpec axisDCAT{3001, -0.005, 30.005, "cm"};
    const AxisSpec axispT{501, -0.05, 50.05, "pT(GeV/c)"};
    const AxisSpec axisP{501, -0.05, 50.05, "p(GeV/c)"};
    const AxisSpec axisZ{6201, -600.05, 20.05, "produced_z(cm)"};

    histos.add("AllMuon_pT", "AllMuon_pT", kTH1F, {axispT});
    histos.add("Beauty_pT", "Beauty_pT", kTH1F, {axispT});
    histos.add("npCharm_pT", "npCharm_pT", kTH1F, {axispT});
    histos.add("pCharm_pT", "pCharm_pT", kTH1F, {axispT});

    histos.add("True_pCharm_DCAT_pT", "True_pCharm_DCAT_pT", kTH1F, {axispT});
    histos.add("Fake_pCharm_DCAT_pT", "Fake_pCharm_DCAT_pT", kTH1F, {axispT});
    histos.add("pCharm_DCAT_PCA_pT", "pCharm_DCAT_PCA_pT", kTH1F, {axispT});

    histos.add("BG_pT", "BG_pT", kTH1F, {axispT});
    histos.add("BG_DCAT_pT", "BG_DCAT_pT", kTH1F, {axispT});
    histos.add("BG_DCAT_PCA_pT", "BG_DCAT_PCA_pT", kTH1F, {axispT});
    histos.add("BG_p_vz", "BG_p_vz", kTH2F, {axisP, axisZ});
    histos.add("BG_vz", "BG_vz", kTH1F, {axisZ});
    histos.add("BG_DCAT", "BG_DCAT", kTH1F, {axisDCAT});
  }

  void initCCDB(ExtBCs::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current()
              << " A for run " << bc.runNumber()
              << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    runNumber = bc.runNumber();

    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    Bz = field->getBz(centerMFT);
    Bz_global = field->getBz(centerGlobal);
    LOG(info) << "The field at the center of the MFT is Bz = " << Bz << " and Global is Bz_Global = " << Bz_global;
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions,
							 soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
							 soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
							 aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
							 aod::McCollisions const&,
               ExtBCs const& ebcs,
               aod::MCPair const& truepairs)
  {
    //const double k_mass = 0.493677;
    //const double mu_mass = 0.105658;
    initCCDB(ebcs.begin());

    vector<uint64_t> ambTrackIds;
    ambTrackIds.clear();
    for(auto& amfttrack : amfttracks){
      ambTrackIds.push_back(amfttrack.mfttrackId());
    }

    // For LHC22f5a & LHC22f5b
    if(truepairs.size()!=0){
      for(auto& truepair : truepairs){
        double Col_x, Col_y, Col_z;
        for(auto& t_col : collisions){
          if(!t_col.has_mcCollision()) continue;
          if(t_col.mcCollisionId()==truepair.colid()){
            Col_x = t_col.posX();
            Col_y = t_col.posY();
            Col_z = t_col.posZ();
          }
        }
        for(auto& fwdtrack : fwdtracks){
          if(!fwdtrack.has_mcParticle()) continue;
          auto mcParticle_fwd = fwdtrack.mcParticle();
          if(mcParticle_fwd.mcCollisionId()!=truepair.colid()) continue;

          double mudcaX,mudcaY,mumftX,mumftY,mumftZ,mueta,muphi;
          double fwdp, fwdpt, fwdpx, fwdpy, fwdpz;
          int PartType = 0;
          bool Signal_Muon = false;
          bool BG_NonPair_Muon = false;
          int trackId = 0;

          SMatrix5 mupars(truepair.x(), truepair.y(), truepair.phi(), truepair.tgl(), truepair.signed1Pt());
          vector<double> muv2;
          SMatrix55 mucovs2(muv2.begin(), muv2.end());
          o2::track::TrackParCovFwd mupars2{truepair.z(), mupars, mucovs2, truepair.chi2()};

          if(mcParticle_fwd.globalIndex()==truepair.muid()){
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

            histos.fill(HIST("AllMuon_pT"), fwdpt);
            if(truepair.particletype()==1){
              histos.fill(HIST("Beauty_pT"), fwdpt);
            }else if(truepair.particletype()==2){
              histos.fill(HIST("npCharm_pT"), fwdpt);
            }else if(truepair.particletype()==3){
              histos.fill(HIST("pCharm_pT"), fwdpt);
            }else if(truepair.particletype()==4){
              histos.fill(HIST("BG_pT"), fwdpt);
            }

            //Only DCAT Cut
            double fwd_dcat = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
            if(fwd_dcat<DCAT_max && fwd_dcat>DCAT_min){
              if(truepair.particletype()==3){
                histos.fill(HIST("True_pCharm_DCAT_pT"), fwdpt);
              }else{
                histos.fill(HIST("Fake_pCharm_DCAT_pT"), fwdpt);
              }
            }
          }

          if(mcParticle_fwd.globalIndex()!=truepair.muid()){
            BG_NonPair_Muon = true;
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
            double mc_p = mcParticle_fwd.p();
            double mc_z = mcParticle_fwd.vz();

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
              histos.fill(HIST("Beauty_pT"), fwdpt);
            }else if(HasLightParent==false && HasBeautyParent==true && HasCharmParent==true){
              PartType = 2; // Non Prompt Charm Muon
              histos.fill(HIST("npCharm_pT"), fwdpt);
            }else if(HasLightParent==false && HasBeautyParent==false && HasCharmParent==true){
              PartType = 3; // Prompt Charm Muon
              histos.fill(HIST("pCharm_pT"), fwdpt);
            }else if(HasLightParent==true && IsSecondary==false){
              PartType = 4; // LF Muon(BG)
            }

            histos.fill(HIST("AllMuon_pT"), fwdpt);

            //Only DCAT Cut
            double fwd_dcat = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
            if(PartType==4){
              if(fwd_dcat<DCAT_max && fwd_dcat>DCAT_min){
                histos.fill(HIST("Fake_pCharm_DCAT_pT"), fwdpt);
                histos.fill(HIST("BG_DCAT_pT"), fwdpt);
              }
              histos.fill(HIST("BG_p_vz"), mc_p, mc_z);
              histos.fill(HIST("BG_vz"), mc_z);
              histos.fill(HIST("BG_pT"), fwdpt);
              histos.fill(HIST("BG_DCAT"), fwd_dcat);
            }if(PartType==3){
              if(fwd_dcat<DCAT_max && fwd_dcat>DCAT_min){
                histos.fill(HIST("True_pCharm_DCAT_pT"), fwdpt);
              }
            }
          }

          vector<vector<double>> data_r;
          if(Signal_Muon==true || BG_NonPair_Muon==true){
            for(auto& mfttrack : mfttracks){
              if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
              if(mfttrack.mcParticle().mcCollisionId()==truepair.colid() && mfttrack.mcParticleId()!=truepair.muid()){
                if(mfttrack.globalIndex()==trackId) continue;
                auto mcParticle_mft = mfttrack.mcParticle();
                if(fabs(mcParticle_mft.pdgCode())==2212) continue;
                double mftchi2 = mfttrack.chi2();
                if(mftchi2>1) continue;
                double caneta = mfttrack.eta();
                double canphi = mfttrack.phi();

                SMatrix5 mftpars_2(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                vector<double> mftv2;
                SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
                o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars_2, mftcovs2, mftchi2};
                
                double candcaX = mftpars2.getX();
                double candcaY = mftpars2.getY();
                double canmftX = mfttrack.x();
                double canmftY = mfttrack.y();
                double canmftZ = mfttrack.z();

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
                double deltaR = sqrt(pow(mueta-caneta,2)+pow(muphi-canphi,2));
                double dcat_mu = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
                double dcat_pair = sqrt(pow(candcaX-Col_x,2)+pow(candcaY-Col_y,2));

                double pair_id = mfttrack.mcParticleId();

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

                vector<double> pardata_r = {r_xyz, fabs(mcParticle_mft.pdgCode()), cosxy, add_dist, deltaR, pcaD, dcat_mu, dcat_pair, mftchi2, ans};
                data_r.push_back(pardata_r);
              }
            }
            if(data_r.empty()) continue;
            sort(data_r.begin(), data_r.end());
          }

        }
      }
    }else{
      for(auto& fwdtrack : fwdtracks){
        if(fwdtrack.has_collision() && fwdtrack.has_mcParticle()){
          if(fwdtrack.eta()<=(-3.6) || fwdtrack.eta()>=(-2.5)) continue;
          auto mcParticle_fwd = fwdtrack.mcParticle();
          if(fabs(mcParticle_fwd.pdgCode())!=13) continue;
          LOGF(info, "Process is working...");
          double Col_x, Col_y, Col_z;
          for(auto& t_col : collisions){
            if(!t_col.has_mcCollision()) continue;
            if(t_col.mcCollisionId()==mcParticle_fwd.mcCollisionId()){
              Col_x = t_col.posX();
              Col_y = t_col.posY();
              Col_z = t_col.posZ();
            }
          }
          if(fabs(Col_z)>10) continue;

          auto trackId = fwdtrack.matchMFTTrackId();
          SMatrix5 mupars_bg(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
          vector<double> muv2;
          SMatrix55 mucovs2(muv2.begin(), muv2.end());
          o2::track::TrackParCovFwd mupars2{fwdtrack.z(), mupars_bg, mucovs2, fwdtrack.chi2()};
          mupars2.propagateToZlinear(Col_z);
          double bg_mudcaX = mupars2.getX();
          double bg_mudcaY = mupars2.getY();
          double bg_mumftX = fwdtrack.x();
          double bg_mumftY = fwdtrack.y();
          double bg_mumftZ = fwdtrack.z();
          double bg_eta = fwdtrack.eta();
          double bg_phi = fwdtrack.phi();
          double bg_pt = fwdtrack.pt();
          double mc_bg_p = mcParticle_fwd.p();
          double mc_bg_z = mcParticle_fwd.vz();

          vector<double> bg_muv2;
          SMatrix55 bg_mucovs2(bg_muv2.begin(), bg_muv2.end());
          SMatrix5 bg_mupars(bg_mumftX, bg_mumftY, fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
          o2::track::TrackParCovFwd bg_mupars2{bg_mumftZ, bg_mupars, bg_mucovs2, fwdtrack.chi2()};

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
            histos.fill(HIST("Beauty_pT"), bg_pt);
          }else if(HasLightParent==false && HasBeautyParent==true && HasCharmParent==true){
            PartType = 2; // Non Prompt Charm Muon
            histos.fill(HIST("npCharm_pT"), bg_pt);
          }else if(HasLightParent==false && HasBeautyParent==false && HasCharmParent==true){
            PartType = 3; // Prompt Charm Muon
            histos.fill(HIST("pCharm_pT"), bg_pt);
          }else if(HasLightParent==true && IsSecondary==false){
            PartType = 4; // LF Muon(BG)
          }

          histos.fill(HIST("AllMuon_pT"), bg_pt);

          if(PartType==4){
            histos.fill(HIST("BG_pT"), bg_pt);
            double fwd_dcat = sqrt(pow(bg_mudcaX-Col_x,2)+pow(bg_mudcaY-Col_y,2));
            if(fwd_dcat<DCAT_max && fwd_dcat>DCAT_min){
              histos.fill(HIST("Fake_pCharm_DCAT_pT"), bg_pt);
              histos.fill(HIST("BG_DCAT_pT"), bg_pt);
            }
            histos.fill(HIST("BG_vz"), mc_bg_z);
            histos.fill(HIST("BG_p_vz"), mc_bg_p, mc_bg_z);
            histos.fill(HIST("BG_DCAT"), fwd_dcat);
            
            vector<vector<double>> data_bg;
            for(auto& mfttrack : mfttracks){
              if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
              if(mfttrack.mcParticle().mcCollisionId()!=fwdtrack.mcParticle().mcCollisionId()) continue;
              if(mfttrack.globalIndex()==trackId) continue;
              auto mcParticle_mft = mfttrack.mcParticle();
              if(fabs(mcParticle_mft.pdgCode())==2212) continue;
              double mftchi2 = mfttrack.chi2();
              if(mftchi2>1) continue;
              SMatrix5 mftpars_2(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
              vector<double> mftv2;
              SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
              o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars_2, mftcovs2, mftchi2};
              
              double candcaX = mftpars2.getX();
              double candcaY = mftpars2.getY();
              double canmftX = mfttrack.x();
              double canmftY = mfttrack.y();
              double canmftZ = mfttrack.z();

              double mfteta = mfttrack.eta();
              double mftphi = mfttrack.phi();

              auto unit_Na = sqrt(pow(bg_mumftX-bg_mudcaX,2)+pow(bg_mumftY-bg_mudcaY,2)+pow(bg_mumftZ-Col_z,2));
              auto unit_Nc = sqrt(pow(canmftX-candcaX,2)+pow(canmftY-candcaY,2)+pow(canmftZ-Col_z,2));
              auto Nax = (bg_mumftX-bg_mudcaX)/unit_Na;
              auto Nay = (bg_mumftY-bg_mudcaY)/unit_Na;
              auto Naz = (bg_mumftZ-Col_z)/unit_Na;
              auto Ncx = (canmftX-candcaX)/unit_Nc;
              auto Ncy = (canmftY-candcaY)/unit_Nc;
              auto Ncz = (canmftZ-Col_z)/unit_Nc;
              auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
              auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
              auto A3 = (bg_mudcaX-candcaX)*Nax + (bg_mudcaY-candcaY)*Nay + (Col_z-Col_z)*Naz;
              auto B1 = A2;
              auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
              auto B3 = (candcaX-bg_mudcaX)*Ncx + (candcaY-bg_mudcaY)*Ncy + (Col_z-Col_z)*Ncz;
              auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
              auto s = -((A2*t+A3)/A1);

              double predict_bg_mux = bg_mudcaX + s*Nax;
              double predict_bg_muy = bg_mudcaY + s*Nay;
              double predict_bg_muz = Col_z + s*Naz;
              double predict_canx = candcaX + t*Ncx;
              double predict_cany = candcaY + t*Ncy;
              double predict_canz = Col_z + t*Ncz;
              double r_xyz = sqrt(pow(predict_canx-predict_bg_mux,2) + pow(predict_cany-predict_bg_muy,2) + pow(predict_canz-predict_bg_muz,2));

              auto vecx_bg_mu = bg_mumftX - bg_mudcaX;
              auto vecy_bg_mu = bg_mumftY - bg_mudcaY;
              auto vecz_bg_mu = bg_mumftZ - Col_z;
              auto vecx_can = canmftX - candcaX;
              auto vecy_can = canmftY - candcaY;
              auto vecz_can = canmftZ - Col_z;
              auto cosxy = (vecx_bg_mu*vecx_can + vecy_bg_mu*vecy_can + vecz_bg_mu*vecz_can)/((sqrt(pow(vecx_bg_mu,2)+pow(vecy_bg_mu,2)+pow(vecz_bg_mu,2)))*(sqrt(pow(vecx_can,2)+pow(vecy_can,2)+pow(vecz_can,2))));

              auto pcaX = ((predict_bg_mux+predict_canx)/2)-Col_x;
              auto pcaY = ((predict_bg_muy+predict_cany)/2)-Col_y;
              auto pcaZ = ((predict_bg_muz+predict_canz)/2)-Col_z;
              auto pcaD = sqrt(pow(pcaX,2)+pow(pcaY,2)+pow(pcaZ,2));
              double deltaR = sqrt(pow(bg_eta-mfteta,2)+pow(bg_phi-mftphi,2));
              double dcat_bg_mu = sqrt(pow(bg_mudcaX-Col_x,2)+pow(bg_mudcaY-Col_y,2));
              double dcat_pair = sqrt(pow(candcaX-Col_x,2)+pow(candcaY-Col_y,2));

              double add_dist = 0;
              for(int i=0; i<5000; i++){
                auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                if((cut_z-Col_z)>0) break;
                bg_mupars2.propagateToZlinear(cut_z);
                mftpars2.propagateToZlinear(cut_z);
                auto diff = sqrt(pow(mftpars2.getX()-bg_mupars2.getX(),2)+pow(mftpars2.getY()-bg_mupars2.getY(),2));
                add_dist += diff;
              }
              vector<double> pardata_r = {r_xyz, fabs(mcParticle_mft.pdgCode()), cosxy, add_dist, deltaR, pcaD, dcat_bg_mu, dcat_pair, mftchi2};
              data_bg.push_back(pardata_r);
            }
            sort(data_bg.begin(), data_bg.end());
          }
        }else if(!fwdtrack.has_collision() || !fwdtrack.has_mcParticle()){
          if(!fwdtrack.has_collision() && fwdtrack.has_mcParticle()){
            if(fwdtrack.eta()<=(-3.6) || fwdtrack.eta()>=(-2.5)) continue;
            auto mcParticle_fwd = fwdtrack.mcParticle();
            if(fabs(mcParticle_fwd.pdgCode())!=13) continue;
            LOGF(info, "Second Process is working...");
            double Col_x, Col_y, Col_z;
            for(auto& t_col : collisions){
              if(!t_col.has_mcCollision()) continue;
              if(t_col.mcCollisionId()==mcParticle_fwd.mcCollisionId()){
                Col_x = t_col.posX();
                Col_y = t_col.posY();
                Col_z = t_col.posZ();
              }
            }
            if(fabs(Col_z)>10) continue;

            auto trackId = fwdtrack.matchMFTTrackId();
            SMatrix5 mupars_bg(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
            vector<double> muv2;
            SMatrix55 mucovs2(muv2.begin(), muv2.end());
            o2::track::TrackParCovFwd mupars2{fwdtrack.z(), mupars_bg, mucovs2, fwdtrack.chi2()};
            mupars2.propagateToZlinear(Col_z);
            double bg_mudcaX = mupars2.getX();
            double bg_mudcaY = mupars2.getY();
            double bg_mumftX = fwdtrack.x();
            double bg_mumftY = fwdtrack.y();
            double bg_mumftZ = fwdtrack.z();
            double bg_eta = fwdtrack.eta();
            double bg_phi = fwdtrack.phi();
            double bg_pt = fwdtrack.pt();
            double mc_bg_p = mcParticle_fwd.p();
            double mc_bg_z = mcParticle_fwd.vz();

            vector<double> bg_muv2;
            SMatrix55 bg_mucovs2(bg_muv2.begin(), bg_muv2.end());
            SMatrix5 bg_mupars(bg_mumftX, bg_mumftY, fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
            o2::track::TrackParCovFwd bg_mupars2{bg_mumftZ, bg_mupars, bg_mucovs2, fwdtrack.chi2()};

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
              histos.fill(HIST("Beauty_pT"), bg_pt);
            }else if(HasLightParent==false && HasBeautyParent==true && HasCharmParent==true){
              PartType = 2; // Non Prompt Charm Muon
              histos.fill(HIST("npCharm_pT"), bg_pt);
            }else if(HasLightParent==false && HasBeautyParent==false && HasCharmParent==true){
              PartType = 3; // Prompt Charm Muon
              histos.fill(HIST("pCharm_pT"), bg_pt);
            }else if(HasLightParent==true && IsSecondary==false){
              PartType = 4; // LF Muon(BG)
            }
            histos.fill(HIST("AllMuon_pT"), bg_pt);
            if(PartType==4){
              histos.fill(HIST("BG_pT"), bg_pt);
              double fwd_dcat = sqrt(pow(bg_mudcaX-Col_x,2)+pow(bg_mudcaY-Col_y,2));
              if(fwd_dcat<DCAT_max && fwd_dcat>DCAT_min){
                histos.fill(HIST("Fake_pCharm_DCAT_pT"), bg_pt);
                histos.fill(HIST("BG_DCAT_pT"), bg_pt);
              }
              histos.fill(HIST("BG_vz"), mc_bg_z);
              histos.fill(HIST("BG_p_vz"), mc_bg_p, mc_bg_z);
              histos.fill(HIST("BG_DCAT"), fwd_dcat);
              
              vector<vector<double>> data_bg;
              for(auto& mfttrack : mfttracks){
                if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
                if(mfttrack.mcParticle().mcCollisionId()!=fwdtrack.mcParticle().mcCollisionId()) continue;
                if(mfttrack.globalIndex()==trackId) continue;
                auto mcParticle_mft = mfttrack.mcParticle();
                if(fabs(mcParticle_mft.pdgCode())==2212) continue;
                double mftchi2 = mfttrack.chi2();
                if(mftchi2>1) continue;
                SMatrix5 mftpars_2(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                vector<double> mftv2;
                SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
                o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars_2, mftcovs2, mftchi2};
                
                double candcaX = mftpars2.getX();
                double candcaY = mftpars2.getY();
                double canmftX = mfttrack.x();
                double canmftY = mfttrack.y();
                double canmftZ = mfttrack.z();

                double mfteta = mfttrack.eta();
                double mftphi = mfttrack.phi();

                auto unit_Na = sqrt(pow(bg_mumftX-bg_mudcaX,2)+pow(bg_mumftY-bg_mudcaY,2)+pow(bg_mumftZ-Col_z,2));
                auto unit_Nc = sqrt(pow(canmftX-candcaX,2)+pow(canmftY-candcaY,2)+pow(canmftZ-Col_z,2));
                auto Nax = (bg_mumftX-bg_mudcaX)/unit_Na;
                auto Nay = (bg_mumftY-bg_mudcaY)/unit_Na;
                auto Naz = (bg_mumftZ-Col_z)/unit_Na;
                auto Ncx = (canmftX-candcaX)/unit_Nc;
                auto Ncy = (canmftY-candcaY)/unit_Nc;
                auto Ncz = (canmftZ-Col_z)/unit_Nc;
                auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
                auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
                auto A3 = (bg_mudcaX-candcaX)*Nax + (bg_mudcaY-candcaY)*Nay + (Col_z-Col_z)*Naz;
                auto B1 = A2;
                auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
                auto B3 = (candcaX-bg_mudcaX)*Ncx + (candcaY-bg_mudcaY)*Ncy + (Col_z-Col_z)*Ncz;
                auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
                auto s = -((A2*t+A3)/A1);

                double predict_bg_mux = bg_mudcaX + s*Nax;
                double predict_bg_muy = bg_mudcaY + s*Nay;
                double predict_bg_muz = Col_z + s*Naz;
                double predict_canx = candcaX + t*Ncx;
                double predict_cany = candcaY + t*Ncy;
                double predict_canz = Col_z + t*Ncz;
                double r_xyz = sqrt(pow(predict_canx-predict_bg_mux,2) + pow(predict_cany-predict_bg_muy,2) + pow(predict_canz-predict_bg_muz,2));

                auto vecx_bg_mu = bg_mumftX - bg_mudcaX;
                auto vecy_bg_mu = bg_mumftY - bg_mudcaY;
                auto vecz_bg_mu = bg_mumftZ - Col_z;
                auto vecx_can = canmftX - candcaX;
                auto vecy_can = canmftY - candcaY;
                auto vecz_can = canmftZ - Col_z;
                auto cosxy = (vecx_bg_mu*vecx_can + vecy_bg_mu*vecy_can + vecz_bg_mu*vecz_can)/((sqrt(pow(vecx_bg_mu,2)+pow(vecy_bg_mu,2)+pow(vecz_bg_mu,2)))*(sqrt(pow(vecx_can,2)+pow(vecy_can,2)+pow(vecz_can,2))));

                auto pcaX = ((predict_bg_mux+predict_canx)/2)-Col_x;
                auto pcaY = ((predict_bg_muy+predict_cany)/2)-Col_y;
                auto pcaZ = ((predict_bg_muz+predict_canz)/2)-Col_z;
                auto pcaD = sqrt(pow(pcaX,2)+pow(pcaY,2)+pow(pcaZ,2));
                double deltaR = sqrt(pow(bg_eta-mfteta,2)+pow(bg_phi-mftphi,2));
                double dcat_bg_mu = sqrt(pow(bg_mudcaX-Col_x,2)+pow(bg_mudcaY-Col_y,2));
                double dcat_pair = sqrt(pow(candcaX-Col_x,2)+pow(candcaY-Col_y,2));

                double add_dist = 0;
                for(int i=0; i<5000; i++){
                  auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                  if((cut_z-Col_z)>0) break;
                  bg_mupars2.propagateToZlinear(cut_z);
                  mftpars2.propagateToZlinear(cut_z);
                  auto diff = sqrt(pow(mftpars2.getX()-bg_mupars2.getX(),2)+pow(mftpars2.getY()-bg_mupars2.getY(),2));
                  add_dist += diff;
                }
                vector<double> pardata_r = {r_xyz, fabs(mcParticle_mft.pdgCode()), cosxy, add_dist, deltaR, pcaD, dcat_bg_mu, dcat_pair, mftchi2};
                data_bg.push_back(pardata_r);
              }
              sort(data_bg.begin(), data_bg.end());
            }
          }else if(!fwdtrack.has_mcParticle()){
            LOGF(info, "No MC Particle");
            continue;
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