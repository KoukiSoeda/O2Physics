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

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "TMath.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include "Framework/ASoAHelpers.h"
#include <math.h>
#include <TLorentzVector.h>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace evsel;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
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
DECLARE_SOA_COLUMN(MuVec, muvec, int);
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
DECLARE_SOA_COLUMN(CosSim, cossim, double);
DECLARE_SOA_COLUMN(TotalD, totald, double);
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
DECLARE_SOA_COLUMN(PairMftX, pairx, double);
DECLARE_SOA_COLUMN(PairMftY, pairy, double);
DECLARE_SOA_COLUMN(PairMftZ, pairz, double);
DECLARE_SOA_COLUMN(PairChi2, pairchi2, double);
DECLARE_SOA_COLUMN(PairPhi, pairphi, double);
DECLARE_SOA_COLUMN(PairTgl, pairtgl, double);
DECLARE_SOA_COLUMN(PairSignedPt, pairsigned1Pt, double);
DECLARE_SOA_COLUMN(ParticleType, particletype, int);
}
namespace muinfo
{
DECLARE_SOA_COLUMN(ColId, colid, int);
DECLARE_SOA_COLUMN(MuId, muid, int);
DECLARE_SOA_COLUMN(MftX, x, double);
DECLARE_SOA_COLUMN(MftY, y, double);
DECLARE_SOA_COLUMN(MftZ, z, double);
DECLARE_SOA_COLUMN(Chi2, chi2, double);
DECLARE_SOA_COLUMN(Phi, phi, double);
DECLARE_SOA_COLUMN(Tgl, tgl, double);
DECLARE_SOA_COLUMN(SignedPt, signed1Pt, double);
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
                  truepair::MuVec,
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
                  truepair::CosSim,
                  truepair::TotalD,
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
                  truepair::PairMftX,
                  truepair::PairMftY,
                  truepair::PairMftZ,
                  truepair::PairChi2,
                  truepair::PairPhi,
                  truepair::PairTgl,
                  truepair::PairSignedPt,
                  truepair::ParticleType);
DECLARE_SOA_TABLE(MuInfo, "ADO", "MUINFO",
                  muinfo::ColId,
                  muinfo::MuId,
                  muinfo::MftX,
                  muinfo::MftY,
                  muinfo::MftZ,
                  muinfo::Chi2,
                  muinfo::Phi,
                  muinfo::Tgl,
                  muinfo::SignedPt,
                  muinfo::ParticleType);
}

struct McInformation{
  Produces<aod::MCPair> mcpairtable;
  Produces<aod::MuInfo> muinfotable;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsZ{"nBinsZ", 110010, "N bins in z-axis"};

  void init(InitContext const&){
    const AxisSpec axisParType{5, -0.5, 4.5, "ParticleType"};
    const AxisSpec axisDCAT{100010, -0.0005, 10.0005, "DCAT"};
    const AxisSpec axispT{501, -0.05, 50.05, "pT(GeV/c)"};
    const AxisSpec axisZ{nBinsZ, -60.0005, 50.0005, "z(cm)"};
    const AxisSpec axisAdd{5000, 0, 500, "sum(cm)"};
    const AxisSpec axisPCAD{30001, -0.0005, 30.0005, "cm"};
    const AxisSpec axisVecCor{2, -0.5, 1.5, ""};
    const AxisSpec axisDelta{5001, -0.005, 50.005};
    const AxisSpec axisEta{111, -10.05, -1.05, "#eta"};
    const AxisSpec axisPhi{700, -0.5, 7, "#phi"};
    const AxisSpec axisMass{10001, -0.0005, 10.0005, "GeV/c^2"};
    const AxisSpec axisCos{1000, 0, 1, ""};

// To do -> create histograms for pcad by each pT regions.

    histos.add("ParticleType", "ParticleType", kTH1F, {axisParType});
    histos.add("DecayL", "DecayL", kTH1F, {axisZ});
    histos.add("VectorCorrelation", "VectorCorrelation", kTH2F, {axisVecCor, axisVecCor});
    histos.add("Delta_XY_mu", "Delta_XY_mu", kTH1F, {axisDelta});
    histos.add("Delta_XY_pair", "Delta_XY_pair", kTH1F, {axisDelta});
    histos.add("Delta_XY_BG", "Delta_XY_BG", kTH1F, {axisDelta});

    histos.add("Beauty_mu_DCAT", "Beauty_mu_DCAT", kTH1F, {axisDCAT});
    histos.add("Beauty_mu_DCAT_1_pT_15", "Beauty_mu_DCAT_1_pT_15", kTH1F, {axisDCAT});
    histos.add("Beauty_mu_DCAT_15_pT_2", "Beauty_mu_DCAT_15_pT_2", kTH1F, {axisDCAT});
    histos.add("Beauty_mu_DCAT_2_pT_3", "Beauty_mu_DCAT_2_pT_3", kTH1F, {axisDCAT});
    histos.add("Beauty_mu_DCAT_3_pT_4", "Beauty_mu_DCAT_3_pT_4", kTH1F, {axisDCAT});
    histos.add("Beauty_mu_DCAT_4_pT_6", "Beauty_mu_DCAT_4_pT_6", kTH1F, {axisDCAT});
    histos.add("Beauty_pair_mu_DCAT", "Beauty_pair_mu_DCAT", kTH1F, {axisDCAT});
    histos.add("Beauty_pair_DCAT", "Beauty_pair_DCAT", kTH1F, {axisDCAT});
    histos.add("Beauty_vs_mu_pt", "Beauty_vs_mu_pt", kTH2F, {axispT, axispT});
    histos.add("Beauty_mupt_totald", "Beauty_mupt_totald", kTH2F, {axispT, axisAdd});
    histos.add("Beauty_pcad", "Beauty_pcad", kTH1F, {axisPCAD});
    histos.add("Beauty_totald", "Beauty_totald", kTH1F, {axisAdd});

    histos.add("pCharm_mu_DCAT", "pCharm_mu_DCAT", kTH1F, {axisDCAT});
    histos.add("pCharm_mu_DCAT_1_pT_15", "pCharm_mu_DCAT_1_pT_15", kTH1F, {axisDCAT});
    histos.add("pCharm_mu_DCAT_15_pT_2", "pCharm_mu_DCAT_15_pT_2", kTH1F, {axisDCAT});
    histos.add("pCharm_mu_DCAT_2_pT_3", "pCharm_mu_DCAT_2_pT_3", kTH1F, {axisDCAT});
    histos.add("pCharm_mu_DCAT_3_pT_4", "pCharm_mu_DCAT_3_pT_4", kTH1F, {axisDCAT});
    histos.add("pCharm_mu_DCAT_4_pT_6", "pCharm_mu_DCAT_4_pT_6", kTH1F, {axisDCAT});
    histos.add("pCharm_pair_mu_DCAT", "pCharm_pair_mu_DCAT", kTH1F, {axisDCAT});
    histos.add("pCharm_pair_DCAT", "pCharm_pair_DCAT", kTH1F, {axisDCAT});
    histos.add("pCharm_mu_eta", "pCharm_mu_eta", kTH1F, {axisEta});
    histos.add("pCharm_pair_eta", "pCharm_pair_eta", kTH1F, {axisEta});
    histos.add("pCharm_pair_phi", "pCharm_pair_phi", kTH1F, {axisPhi});
    histos.add("pCharm_vs_mu_pt", "pCharm_vs_mu_pt", kTH2F, {axispT, axispT});
    histos.add("pCharm_mupt_totald", "pCharm_mupt_totald", kTH2F, {axispT, axisAdd});
    histos.add("pCharm_pcad", "pCharm_pcad", kTH1F, {axisPCAD});
    histos.add("pCharm_pcar", "pCharm_pcar", kTH1F, {axisPCAD});
    histos.add("pCharm_totald", "pCharm_totald", kTH1F, {axisAdd});
    histos.add("pCharm_DCAXY_Diff", "pCharm_DCAXY_Diff", kTH1F, {axisDCAT});
    histos.add("pCharm_InvMass", "pCharm_InvMass", kTH1F, {axisMass});
    histos.add("pCharm_Mass_wo_neutrino", "pCharm_Mass_wo_neutrino", kTH1F, {axisMass});
    histos.add("pCharm_Mass_wo_neutrino_MFT", "pCharm_Mass_wo_neutrino_MFT", kTH1F, {axisMass});
    histos.add("pCharm_InvMass_Fake", "pCharm_InvMass_Fake", kTH1F, {axisMass});
    histos.add("pCharm_InvMass_Fake_MFT", "pCharm_InvMass_Fake_MFT", kTH1F, {axisMass});
    histos.add("CosSim_true", "CosSim_true", kTH1F, {axisCos});

    histos.add("npCharm_mu_DCAT", "npCharm_mu_DCAT", kTH1F, {axisDCAT});
    histos.add("npCharm_mu_DCAT_1_pT_15", "npCharm_mu_DCAT_1_pT_15", kTH1F, {axisDCAT});
    histos.add("npCharm_mu_DCAT_15_pT_2", "npCharm_mu_DCAT_15_pT_2", kTH1F, {axisDCAT});
    histos.add("npCharm_mu_DCAT_2_pT_3", "npCharm_mu_DCAT_2_pT_3", kTH1F, {axisDCAT});
    histos.add("npCharm_mu_DCAT_3_pT_4", "npCharm_mu_DCAT_3_pT_4", kTH1F, {axisDCAT});
    histos.add("npCharm_mu_DCAT_4_pT_6", "npCharm_mu_DCAT_4_pT_6", kTH1F, {axisDCAT});
    histos.add("npCharm_pair_mu_DCAT", "npCharm_pair_mu_DCAT", kTH1F, {axisDCAT});
    histos.add("npCharm_pair_DCAT", "npCharm_pair_DCAT", kTH1F, {axisDCAT});
    histos.add("npCharm_mupt_totald", "npCharm_mupt_totald", kTH2F, {axispT, axisAdd});
    histos.add("npCharm_pcad", "npCharm_pcad", kTH1F, {axisPCAD});
    histos.add("npCharm_totald", "npCharm_totald", kTH1F, {axisAdd});

    histos.add("Light_mu_DCAT", "Light_mu_DCAT", kTH1F, {axisDCAT});
    histos.add("Light_pair_DCAT", "Light_pair_DCAT", kTH1F, {axisDCAT});
    histos.add("Light_pcad", "Light_pcad", kTH1F, {axisPCAD});
    histos.add("Light_totald", "Light_totald", kTH1F, {axisAdd});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
               soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA, aod::FwdTracksCov> const& fwdtracks, 
               soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
               aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
               aod::McCollisions const&)
  {
    const double k_mass = 0.493677;
    const double mu_mass = 0.105658;
    if(collision.has_mcCollision()){
      if(fabs(collision.posZ())<10){
        const int mcCol_id = collision.mcCollisionId();
        const int Col_id = collision.globalIndex();
        const double Col_x = collision.posX();
        const double Col_y = collision.posY();
        const double Col_z = collision.posZ();
        // Create AmbTrackTable
        vector<uint64_t> ambTrackIds;
        ambTrackIds.clear();
        for(auto& amfttrack : amfttracks){
          ambTrackIds.push_back(amfttrack.mfttrackId());
        }
        // N of Daughters and Muon Source
        for(auto& fwdtrack : fwdtracks){
          if(!fwdtrack.has_collision() || !fwdtrack.has_mcParticle()) continue;
          if(fwdtrack.collisionId()!=Col_id) continue;
          if(fwdtrack.trackType()==1 || fwdtrack.trackType()==3 || fwdtrack.trackType()==4 ) continue;
          if(fwdtrack.eta()<=(-3.6) || fwdtrack.eta()>=(-2.5)) continue;
          auto mcParticle_fwd = fwdtrack.mcParticle();
          if(fabs(mcParticle_fwd.pdgCode())!=13) continue;

          int mucolid = fwdtrack.collisionId();
          int mcmuid = mcParticle_fwd.globalIndex();
          auto mumom = mcParticle_fwd.mothers_first_as<aod::McParticles>();
          auto momid = mumom.globalIndex();
          double mompt = mumom.pt();
          auto Daughters = mumom.daughters_as<aod::McParticles>();

          int muid=0, pid_1=0, pid_2=0, ppdg_1=0, ppdg_2=0;
          double mueta=0, paireta=0;
          double sverx, svery, sverz;
          double px_mu, py_mu, pz_mu, e_mu, px_k, py_k, pz_k, e_k, mass_no_neu, inv_mass;
          double px_mu_mft, py_mu_mft, pz_mu_mft, p_mu_mft;
          double px_neu, py_neu, pz_neu, e_neu;
          bool HasPair = false;
          for(auto& Daughter : Daughters){
            int d_pdg = fabs(Daughter.pdgCode());
            if(d_pdg==13 && Daughter.globalIndex()==mcmuid){
              muid = Daughter.globalIndex();
              mueta = Daughter.eta();
              sverx = Daughter.vx();
              svery = Daughter.vy();
              sverz = Daughter.vz();
            }else if(d_pdg!=13 && d_pdg>20 && (d_pdg%2)!=0 && pid_1==0){
              pid_1 = Daughter.globalIndex();
              paireta = Daughter.eta();
              ppdg_1 = d_pdg;
              HasPair = true;
              px_k = Daughter.px();
              py_k = Daughter.py();
              pz_k = Daughter.pz();
              e_k = Daughter.e();
            }else if(d_pdg!=13 && d_pdg>20 && (d_pdg%2)!=0 && pid_1!=0){
              pid_2 = Daughter.globalIndex();
              ppdg_2 = d_pdg;
              px_neu = Daughter.px();
              py_neu = Daughter.py();
              pz_neu = Daughter.pz();
              e_neu = Daughter.e();
            }
          }
          if(fabs(mumom.pdgCode())==421){
            bool hasmu = false;
            bool hask = false;
            bool hasneu = false;
            for(auto& Daughter : Daughters){
              int d_pdg = fabs(Daughter.pdgCode());
              if(d_pdg==13){
                px_mu = Daughter.px();
                py_mu = Daughter.py();
                pz_mu = Daughter.pz();
                e_mu = Daughter.e();
                px_mu_mft = fwdtrack.px();
                py_mu_mft = fwdtrack.py();
                pz_mu_mft = fwdtrack.pz();
                p_mu_mft = fwdtrack.p();
                hasmu=true;
              }else if(d_pdg==14){
                px_neu = Daughter.px();
                py_neu = Daughter.py();
                pz_neu = Daughter.pz();
                e_neu = Daughter.e();
                hasneu=true;
              }else if(d_pdg==321){
                px_k = Daughter.px();
                py_k = Daughter.py();
                pz_k = Daughter.pz();
                e_k = Daughter.e();
                hask=true;
              }
            }
            if(hasmu && hask && hasneu && Daughters.size()==3){
              inv_mass = sqrt(pow(e_mu+e_k+e_neu,2) - (pow(px_mu+px_k+px_neu,2)+pow(py_mu+py_k+py_neu,2)+pow(pz_mu+pz_k+pz_neu,2)));
              histos.fill(HIST("pCharm_InvMass"), inv_mass);
            }
          }

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
              histos.fill(HIST("pCharm_mu_eta"),mueta);
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
            PartType = 4; // LF Muon
          }
          histos.fill(HIST("ParticleType"), PartType);

          if(!HasPair) continue;

          // Compute the PCA between true pairs
          double mudcaX, mudcaY, mumftX, mumftY, mumftZ, muchi2, muphi, mutgl, musignept, mupt, p1chi2, p1phi, p1tgl, p1signedpt;
          int mcmucolid;
          double p1dcaX, p1dcaY, p1mftX, p1mftY, p1mftZ, p1secverx, p1secvery, p1secverz;
          int p1momid, p1pdg, p1colid, mcp1colid;
          double p2dcaX, p2dcaY, p2mftX, p2mftY, p2mftZ, p2secverx, p2secvery, p2secverz;
          int p2momid, p2pdg;
          vector<vector<double>> muondata{};
          vector<vector<double>> pair1data{};
          bool HasMuonTrack = false;
          bool HasPairTrack = false;

          for(auto& mfttrack : mfttracks){
            if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
            if(mfttrack.collisionId()==Col_id){
              auto mcParticle_mft = mfttrack.mcParticle();
              if(mcParticle_mft.globalIndex()==muid){
                auto mft_mom_mu = mcParticle_mft.mothers_first_as<aod::McParticles>();
                if(mft_mom_mu.globalIndex()==momid && mfttrack.collisionId()==mucolid){
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
                  mcmucolid = mcParticle_fwd.mcCollisionId();
                  muchi2 = mftchi2;
                  muphi = mfttrack.phi();
                  mutgl = mfttrack.tgl();
                  musignept = mfttrack.signed1Pt();
                  mupt = mfttrack.pt();
                  HasMuonTrack = true;
                  for(int i=0; i<10000; i++){
                    auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                    if((cut_z-Col_z)>0) break;
                    mftpars1.propagateToZlinear(cut_z);
                    vector<double> temp = {mftpars1.getX(), mftpars1.getY()};
                    muondata.push_back(temp);
                  }
                  double dcat_mu_all = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
                  if(PartType==1){
                    histos.fill(HIST("Beauty_mu_DCAT"), dcat_mu_all);
                    if(mupt>1 && mupt<1.5){
                      histos.fill(HIST("Beauty_mu_DCAT_1_pT_15"), dcat_mu_all);
                    }else if(mupt>=1.5 && mupt<2){
                      histos.fill(HIST("Beauty_mu_DCAT_15_pT_2"), dcat_mu_all);
                    }else if(mupt>=2 && mupt<3){
                      histos.fill(HIST("Beauty_mu_DCAT_2_pT_3"), dcat_mu_all);
                    }else if(mupt>=3 && mupt<4){
                      histos.fill(HIST("Beauty_mu_DCAT_3_pT_4"), dcat_mu_all);
                    }else if(mupt>=4 && mupt<6){
                      histos.fill(HIST("Beauty_mu_DCAT_4_pT_6"), dcat_mu_all);
                    }
                  }
                  if(PartType==2){
                    histos.fill(HIST("npCharm_mu_DCAT"), dcat_mu_all);
                    if(mupt>1 && mupt<1.5){
                      histos.fill(HIST("npCharm_mu_DCAT_1_pT_15"), dcat_mu_all);
                    }else if(mupt>=1.5 && mupt<2){
                      histos.fill(HIST("npCharm_mu_DCAT_15_pT_2"), dcat_mu_all);
                    }else if(mupt>=2 && mupt<3){
                      histos.fill(HIST("npCharm_mu_DCAT_2_pT_3"), dcat_mu_all);
                    }else if(mupt>=3 && mupt<4){
                      histos.fill(HIST("npCharm_mu_DCAT_3_pT_4"), dcat_mu_all);
                    }else if(mupt>=4 && mupt<6){
                      histos.fill(HIST("npCharm_mu_DCAT_4_pT_6"), dcat_mu_all);
                    }
                  }
                  if(PartType==3){
                    histos.fill(HIST("pCharm_mu_DCAT"), dcat_mu_all);
                    if(mupt>1 && mupt<1.5){
                      histos.fill(HIST("pCharm_mu_DCAT_1_pT_15"), dcat_mu_all);
                    }else if(mupt>=1.5 && mupt<2){
                      histos.fill(HIST("pCharm_mu_DCAT_15_pT_2"), dcat_mu_all);
                    }else if(mupt>=2 && mupt<3){
                      histos.fill(HIST("pCharm_mu_DCAT_2_pT_3"), dcat_mu_all);
                    }else if(mupt>=3 && mupt<4){
                      histos.fill(HIST("pCharm_mu_DCAT_3_pT_4"), dcat_mu_all);
                    }else if(mupt>=4 && mupt<6){
                      histos.fill(HIST("pCharm_mu_DCAT_4_pT_6"), dcat_mu_all);
                    }
                  }
                }
              }

              if(mcParticle_mft.globalIndex()!=muid){
                //if(find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack.globalIndex()) != ambTrackIds.end()) continue;
                if(mcParticle_mft.has_mothers()){
                  auto mft_mom = mcParticle_mft.mothers_first_as<aod::McParticles>();

                  if(mft_mom.globalIndex()==momid && mcParticle_mft.globalIndex()==pid_1 && mfttrack.collisionId()==mucolid){
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
                    p1momid = mft_mom.globalIndex();
                    p1colid = mfttrack.collisionId();
                    p1chi2 = mfttrack.chi2();
                    p1phi = mfttrack.phi();
                    p1tgl = mfttrack.tgl();
                    p1signedpt = mfttrack.signed1Pt();
                    mcp1colid = mcParticle_mft.mcCollisionId();
                    p1secverx = mcParticle_mft.vx();
                    p1secvery = mcParticle_mft.vy();
                    p1secverz = mcParticle_mft.vz();
                    p1pdg = fabs(mcParticle_mft.pdgCode());
                    HasPairTrack = true;
                    for(int i=0; i<=10000; i++){
                      auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                      if((cut_z-Col_z)>0) break;
                      mftpars2.propagateToZlinear(cut_z);
                      vector<double> temp2 = {mftpars2.getX(), mftpars2.getY()};
                      pair1data.push_back(temp2);
                    }
                    if(fabs(mcParticle_mft.pdgCode())==321 && fabs(mumom.pdgCode())==421){
                      mass_no_neu = sqrt(pow(e_mu+e_k,2) - (pow(px_mu+px_k,2)+pow(py_mu+py_k,2)+pow(pz_mu+pz_k,2)));
                      histos.fill(HIST("pCharm_Mass_wo_neutrino"), mass_no_neu);
                      
                      double px_k_mft = mfttrack.px();
                      double py_k_mft = mfttrack.py();
                      double pz_k_mft = mfttrack.pz();
                      double p_k_mft = mfttrack.p();
                      double inv_mass_true_mft = sqrt(pow(sqrt(pow(mu_mass,2)+pow(p_mu_mft,2))+sqrt(pow(k_mass,2)+pow(p_k_mft,2)), 2) - (pow(px_mu_mft+px_k_mft,2)+pow(py_mu_mft+py_k_mft,2)+pow(pz_mu_mft+pz_k_mft,2)));
                      histos.fill(HIST("pCharm_Mass_wo_neutrino_MFT"), inv_mass_true_mft);
                    }
                  }
                  if(mft_mom.globalIndex()==momid && mcParticle_mft.globalIndex()==pid_2 && mfttrack.collisionId()==mucolid){
                    double mftchi2 = mfttrack.chi2();
                    SMatrix5 mftpars_2(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                    vector<double> mftv2;
                    SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
                    o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars_2, mftcovs2, mftchi2};
                    mftpars2.propagateToZlinear(Col_z);

                    p2dcaX = mftpars2.getX();
                    p2dcaY = mftpars2.getY();
                    p2mftX = mfttrack.x();
                    p2mftY = mfttrack.y();
                    p2mftZ = mfttrack.z();
                    p2momid = mft_mom.globalIndex();
                    p2secverx = mcParticle_mft.vx();
                    p2secvery = mcParticle_mft.vy();
                    p2secverz = mcParticle_mft.vz();
                    p2pdg = fabs(mcParticle_mft.pdgCode());
                  }
                }

                // BG info 
                if(mcParticle_mft.globalIndex()!=pid_1 && mcParticle_mft.globalIndex()!=pid_2){
                  double mftchi2 = mfttrack.chi2();
                  SMatrix5 mftpars_2(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                  vector<double> mftv2;
                  SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
                  o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars_2, mftcovs2, mftchi2};
                  mftpars2.propagateToZlinear(Col_z);
                  double bg_delta_xy = sqrt(pow(mftpars2.getX()-mfttrack.x(),2)+pow(mftpars2.getY()-mfttrack.y(),2));
                  histos.fill(HIST("Delta_XY_BG"), bg_delta_xy);

                  if(fabs(mcParticle_mft.pdgCode())==211 && fabs(mumom.pdgCode())==421){
                    double px_pi = mcParticle_mft.px();
                    double py_pi = mcParticle_mft.py();
                    double pz_pi = mcParticle_mft.pz();
                    double p_pi = mcParticle_mft.p();
                    double px_pi_mft = mfttrack.px();
                    double py_pi_mft = mfttrack.py();
                    double pz_pi_mft = mfttrack.pz();
                    double p_pi_mft = mfttrack.p();
                    double inv_mass_fake = sqrt(pow(e_mu+sqrt(pow(k_mass,2)+pow(p_pi,2)),2) - (pow(px_mu+px_pi,2)+pow(py_mu+py_pi,2)+pow(pz_mu+pz_pi,2)));
                    histos.fill(HIST("pCharm_InvMass_Fake"), inv_mass_fake);
                    double inv_mass_fake_mft = sqrt(pow(sqrt(pow(mu_mass,2)+pow(p_mu_mft,2))+sqrt(pow(k_mass,2)+pow(p_pi_mft,2)), 2) - (pow(px_mu_mft+px_pi_mft,2)+pow(py_mu_mft+py_pi_mft,2)+pow(pz_mu_mft+pz_pi_mft,2)));
                    histos.fill(HIST("pCharm_InvMass_Fake_MFT"), inv_mass_fake_mft);
                  }
                }
              }
            }
          }

          if(p1momid==momid && p1secverx==sverx && p1secvery==svery && p1secverz==sverz && p1pdg==ppdg_1){
            if(mcCol_id==mcmucolid && mcCol_id==mcp1colid && HasMuonTrack && HasPairTrack){
              muinfotable(mucolid, muid, mumftX, mumftY, mumftZ, muchi2, muphi, mutgl, musignept, PartType);
              // Compute DCAT
              double dcat_mu = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
              double dcat_pair = sqrt(pow(p1dcaX-Col_x,2)+pow(p1dcaY-Col_y,2));
              histos.fill(HIST("DecayL"), p1secverz-Col_z);
              double dcaxy_diff = sqrt(pow(mudcaX-p1dcaX,2)+pow(mudcaY-p1dcaY,2));
              histos.fill(HIST("pCharm_DCAXY_Diff"), dcaxy_diff);

              //Vector Correlation
              bool MuVec = true;
              bool PairVec = true;
              int muon_vec = 1;
              //muon
              double mdx_mu = mudcaX - mumftX;
              double mdy_mu = mudcaY - mumftY;
              double mu_delta_xy = sqrt(pow(mdx_mu,2)+pow(mdy_mu,2));
              histos.fill(HIST("Delta_XY_mu"), mu_delta_xy);
              if((mumftX-Col_x)>0 && (mumftY-Col_y)>0){
                if(mdx_mu>0 || mdy_mu>0) MuVec=false;
              }else if((mumftX-Col_x)>0 && (mumftY-Col_y)<0){
                if(mdx_mu>0 || mdy_mu<0) MuVec=false;
              }else if((mumftX-Col_x)<0 && (mumftY-Col_y)<0){
                if(mdx_mu<0 || mdy_mu<0) MuVec=false;
              }else if((mumftX-Col_x)<0 && (mumftY-Col_y)>0){
                if(mdx_mu<0 || mdy_mu>0) MuVec=false;
              }
              //pair
              double mdx_p1 = p1dcaX - p1mftX;
              double mdy_p1 = p1dcaY - p1mftY;
              double pair_delta_xy = sqrt(pow(mdx_p1,2)+pow(mdy_p1,2));
              histos.fill(HIST("Delta_XY_pair"), pair_delta_xy);
              if((p1mftX-Col_x)>0 && (p1mftY-Col_y)>0){
                if(mdx_p1>0 || mdy_p1>0) PairVec=false;
              }else if((p1mftX-Col_x)>0 && (p1mftY-Col_y)<0){
                if(mdx_p1>0 || mdy_p1<0) PairVec=false;
              }else if((p1mftX-Col_x)<0 && (p1mftY-Col_y)<0){
                if(mdx_p1<0 || mdy_p1<0) PairVec=false;
              }else if((p1mftX-Col_x)<0 && (p1mftY-Col_y)>0){
                if(mdx_p1<0 || mdy_p1>0) PairVec=false;
              }

              if(MuVec==true && PairVec==true){
                histos.fill(HIST("VectorCorrelation"), 1,1);
              }else{
                if(MuVec==true && PairVec==false) histos.fill(HIST("VectorCorrelation"), 1,0);
                if(MuVec==false && PairVec==true){
                  histos.fill(HIST("VectorCorrelation"), 0,1);
                  muon_vec = -1;
                }
                if(MuVec==false && PairVec==false){
                  histos.fill(HIST("VectorCorrelation"), 0,0);
                  muon_vec = -1;
                }
              }

              // 3D Vector Method
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
              histos.fill(HIST("pCharm_pcar"), predict_r);

              // Cosine Similarity
              auto vecx_mu = mumftX - mudcaX;
              auto vecy_mu = mumftY - mudcaY;
              auto vecz_mu = mumftZ - Col_z;
              auto vecx_p1 = p1mftX - p1dcaX;
              auto vecy_p1 = p1mftY - p1dcaY;
              auto vecz_p1 = p1mftZ - Col_z;
              auto cosxy = (vecx_mu*vecx_p1 + vecy_mu*vecy_p1 + vecz_mu*vecz_p1)/((sqrt(pow(vecx_mu,2)+pow(vecy_mu,2)+pow(vecz_mu,2)))*(sqrt(pow(vecx_p1,2)+pow(vecy_p1,2)+pow(vecz_p1,2))));

              histos.fill(HIST("CosSim_true"), cosxy);
              // xy Plane Total Distance
              double dist_sum = 0;
              if(muondata.size()==pair1data.size()){
                for(int i=0; i<muondata.size(); i++){
                  auto dist = sqrt(pow(pair1data.at(i).at(0)-muondata.at(i).at(0),2)+pow(pair1data.at(i).at(1)-muondata.at(i).at(1),2));
                  dist_sum += dist;
                }
              }
              mcpairtable(Col_id, Col_x, Col_y, Col_z, fabs(mumom.pdgCode()), muid, muon_vec, pid_1, p1pdg, p1dcaX, p1dcaY, pid_2, p2pdg, p2dcaX, p2dcaY, sverx, svery, sverz, pcaX, pcaY, pcaZ, pcaD, predict_r, cosxy, dist_sum, mumftX, mumftY, mumftZ, mupt, px_mu_mft, py_mu_mft, pz_mu_mft, p_mu_mft, muchi2, muphi, mutgl, musignept, p1mftX, p1mftY, p1mftZ, p1chi2, p1phi, p1tgl, p1signedpt, PartType);

              if(PartType==1){ // Besuty Decay
                histos.fill(HIST("Beauty_pair_mu_DCAT"), dcat_mu);
                histos.fill(HIST("Beauty_pair_DCAT"), dcat_pair);
                histos.fill(HIST("Beauty_vs_mu_pt"), mupt, mompt);
                histos.fill(HIST("Beauty_mupt_totald"), mupt, dist_sum);
                histos.fill(HIST("Beauty_pcad"), pcaD);
                histos.fill(HIST("Beauty_totald"), dist_sum);
              }else if(PartType==2){ // Non Prompt Charm
                histos.fill(HIST("npCharm_pair_mu_DCAT"), dcat_mu);
                histos.fill(HIST("npCharm_pair_DCAT"), dcat_pair);
                histos.fill(HIST("npCharm_mupt_totald"), mupt, dist_sum);
                histos.fill(HIST("npCharm_pcad"), pcaD);
                histos.fill(HIST("npCharm_totald"), dist_sum);
              }else if(PartType==3){ // Prompt Charm
                histos.fill(HIST("pCharm_pair_mu_DCAT"), dcat_mu);
                histos.fill(HIST("pCharm_pair_DCAT"), dcat_pair);
                histos.fill(HIST("pCharm_vs_mu_pt"), mupt, mompt);
                histos.fill(HIST("pCharm_mupt_totald"), mupt, dist_sum);
                histos.fill(HIST("pCharm_pcad"), pcaD);
                histos.fill(HIST("pCharm_totald"), dist_sum);
              }else if(PartType==4){ // LF Decay
                histos.fill(HIST("Light_mu_DCAT"), dcat_mu);
                histos.fill(HIST("Light_pair_DCAT"), dcat_pair);
                histos.fill(HIST("Light_pcad"), pcaD);
                histos.fill(HIST("Light_totald"), dist_sum);
              }
            }
          }
        }
      }
    }
  }
};

struct MyAnalysisTask{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&){
    const AxisSpec axisPCAD{15001, -0.0005, 15.0005, "cm"};
    const AxisSpec axispT{501, -0.05, 50.05, "pT(GeV/c)"};
    const AxisSpec axisDCAT{100010, -0.0005, 10.0005, "DCAT"};
    const AxisSpec axisCut{11, -0.5, 10.5, "Cut_DCAT(cm)"};
    const AxisSpec axisRank{5, 0.5, 5.5, "Rank"};
    const AxisSpec axisRankv2{100, 0.5, 100.5, "Rank"};
    const AxisSpec axisMult{500, 0.5, 500.5, "Multiplicity"};
    const AxisSpec axisSigma{501, 0.005, 5.005, "Sigma"};
    const AxisSpec axisPCAD_mc{30001, -0.0005, 30.0005, "cm"};
    const AxisSpec axisMass{10001, -0.0005, 10.0005, "GeV/c^2"};
    const AxisSpec axisCos{1000, 0, 1, ""};
    const AxisSpec axisAddDist{10000, 0, 500, "cm"};

    histos.add("AllTrueMatch_pCharm_pcad", "AllTrueMatch_pCharm_pcad", kTH1F, {axisPCAD});
    histos.add("AllFakeMatch_pCharm_pcad", "AllFakeMatch_pCharm_pcad", kTH1F, {axisPCAD});
    histos.add("AllTrueMatch_Beauty_pcad", "AllTrueMatch_Beauty_pcad", kTH1F, {axisPCAD});
    histos.add("AllFakeMatch_Beauty_pcad", "AllFakeMatch_Beauty_pcad", kTH1F, {axisPCAD});
    histos.add("AllTrueMatch_npCharm_pcad", "AllTrueMatch_npCharm_pcad", kTH1F, {axisPCAD});
    histos.add("AllFakeMatch_npCharm_pcad", "AllFakeMatch_npCharm_pcad", kTH1F, {axisPCAD});
    histos.add("TrueMatch_pCharm_pcad", "TrueMatch_pCharm_pcad", kTH1F, {axisPCAD});
    histos.add("TrueMatch_pCharm_pt", "TrueMatch_pCharm_pt", kTH1F, {axispT});
    histos.add("TrueMatch_pCharm_pcad_pt", "TrueMatch_pCharm_pcad_pt", kTH2F, {axisPCAD, axispT});
    histos.add("TrueMatch_npCharm_pcad", "TrueMatch_npCharm_pcad", kTH1F, {axisPCAD});
    histos.add("TrueMatch_Beauty_pcad", "TrueMatch_Beauty_pcad", kTH1F, {axisPCAD});
    histos.add("FakeMatch_pCharm_pcad", "FakeMatch_pCharm_", kTH1F, {axisPCAD});
    histos.add("FakeMatch_pCharm_pt", "FakeMatch_pCharm_pt", kTH1F, {axispT});
    histos.add("FakeMatch_pCharm_pcad_pt", "FakeMatch_pCharm_pcad_pt", kTH2F, {axisPCAD, axispT});
    histos.add("FakeMatch_npCharm_pcad", "FakeMatch_npCharm_pcad", kTH1F, {axisPCAD});
    histos.add("FakeMatch_Beauty_pcad", "FakeMatch_Beauty_pcad", kTH1F, {axisPCAD});

    histos.add("Signal_Dcatcut_pCharm", "Signal_Dcatcut_pCharm", kTH1F, {axisCut});
    histos.add("TrueClosest_Dcatcut_pCharm", "TrueClosest_Dcatcut_pCharm", kTH1F, {axisCut});
    histos.add("FakeClosest_Dcatcut_pCharm", "FakeClosest_Dcatcut_pCharm", kTH1F, {axisCut});
    histos.add("Signal_Dcatcut_Beauty", "Signal_Dcatcut_Beauty", kTH1F, {axisCut});
    histos.add("TrueClosest_Dcatcut_Beauty", "TrueClosest_Dcatcut_Beauty", kTH1F, {axisCut});
    histos.add("FakeClosest_Dcatcut_Beauty", "FakeClosest_Dcatcut_Beauty", kTH1F, {axisCut});
    histos.add("Signal_Dcatcut_npCharm", "Signal_Dcatcut_npCharm", kTH1F, {axisCut});
    histos.add("TrueClosest_Dcatcut_npCharm", "TrueClosest_Dcatcut_npCharm", kTH1F, {axisCut});
    histos.add("FakeClosest_Dcatcut_npCharm", "FakeClosest_Dcatcut_npCharm", kTH1F, {axisCut});
    histos.add("FakeDCAXY_Diff", "FakeDCAXY_Diff", kTH1F, {axisDCAT});
    histos.add("RankEfficiency_pCharm", "RankEfficiency_pCharm", kTH1F, {axisRank});
    histos.add("RankEfficiency_npCharm", "RankEfficiency_npCharm", kTH1F, {axisRank});

    histos.add("c_DCAT_all", "c_DCAT_all", kTH1F, {axispT});
    histos.add("c_DCAT_cut", "c_DCAT_cut", kTH1F, {axispT});
    histos.add("c_DCAT_PCAD_cut", "c_DCAT_PCAD_cut", kTH1F, {axispT});
    histos.add("BG_DCAT_all", "BG_DCAT_all", kTH1F, {axispT});
    histos.add("BG_DCAT_cut", "BG_DCAT_cut", kTH1F, {axispT});
    histos.add("BG_DCAT_PCAD_cut", "BG_DCAT_PCAD_cut", kTH1F, {axispT});
    histos.add("BG_c_PCAD", "BG_c_PCAD", kTH1F, {axisPCAD_mc});
    histos.add("BG_pcar", "BG_pcar", kTH1F, {axisPCAD_mc});
    histos.add("True_pcar", "True_pcar", kTH1F, {axisPCAD_mc});
    histos.add("Fake_pcar", "Fake_pcar", kTH1F, {axisPCAD_mc});
    histos.add("True_InvMass", "True_InvMass", kTH1F, {axisMass});
    histos.add("Fake_InvMass", "Fake_InvMass", kTH1F, {axisMass});
    histos.add("CosSim_fake", "CosSim_fake", kTH1F, {axisCos});
    histos.add("Add_dist_true", "Add_dist_true", kTH1F, {axisAddDist});
    histos.add("Add_dist_fake", "Add_dist_fake", kTH1F, {axisAddDist});
    histos.add("HasPair_wincut", "HasPair_wincut", kTH1F, {axisRankv2});
    histos.add("True_first", "True_first", kTH1F, {axisMult});
    histos.add("Fake_first", "Fake_first", kTH1F, {axisMult});
    histos.add("BG_InvMass", "BG_InvMass", kTH1F, {axisMass});

    histos.add("Sigma_BG_All", "Sigma_BG_All", kTH2F, {axisSigma, axisSigma});

    auto parcut0 = histos.get<TH1>(HIST("Signal_Dcatcut_pCharm"));
    auto* h0 = parcut0->GetXaxis();
    h0->SetBinLabel(1,"0<DCAT<0.02");
    h0->SetBinLabel(2,"0.02<DCAT<0.04");
    h0->SetBinLabel(3,"0.04<DCAT<0.06");
    h0->SetBinLabel(4,"0.06<DCAT<0.08");
    h0->SetBinLabel(5,"0.08<DCAT<0.1");
    h0->SetBinLabel(6,"0.1<DCAT<0.2");
    h0->SetBinLabel(7,"0.2<DCAT<0.3");
    h0->SetBinLabel(8,"0.3<DCAT<0.4");
    h0->SetBinLabel(9,"0.4<DCAT<0.5");

    auto parcut = histos.get<TH1>(HIST("TrueClosest_Dcatcut_pCharm"));
    auto* h1 = parcut->GetXaxis();
    h1->SetBinLabel(1,"0<DCAT<0.02");
    h1->SetBinLabel(2,"0.02<DCAT<0.04");
    h1->SetBinLabel(3,"0.04<DCAT<0.06");
    h1->SetBinLabel(4,"0.06<DCAT<0.08");
    h1->SetBinLabel(5,"0.08<DCAT<0.1");
    h1->SetBinLabel(6,"0.1<DCAT<0.2");
    h1->SetBinLabel(7,"0.2<DCAT<0.3");
    h1->SetBinLabel(8,"0.3<DCAT<0.4");
    h1->SetBinLabel(9,"0.4<DCAT<0.5");

    auto parcut_f = histos.get<TH1>(HIST("FakeClosest_Dcatcut_pCharm"));
    auto* h2 = parcut_f->GetXaxis();
    h2->SetBinLabel(1,"0<DCAT<0.02");
    h2->SetBinLabel(2,"0.02<DCAT<0.04");
    h2->SetBinLabel(3,"0.04<DCAT<0.06");
    h2->SetBinLabel(4,"0.06<DCAT<0.08");
    h2->SetBinLabel(5,"0.08<DCAT<0.1");
    h2->SetBinLabel(6,"0.1<DCAT<0.2");
    h2->SetBinLabel(7,"0.2<DCAT<0.3");
    h2->SetBinLabel(8,"0.3<DCAT<0.4");
    h2->SetBinLabel(9,"0.4<DCAT<0.5");

    auto parcut1 = histos.get<TH1>(HIST("Signal_Dcatcut_Beauty"));
    auto* h3 = parcut1->GetXaxis();
    h3->SetBinLabel(1,"0<DCAT<0.02");
    h3->SetBinLabel(2,"0.02<DCAT<0.04");
    h3->SetBinLabel(3,"0.04<DCAT<0.06");
    h3->SetBinLabel(4,"0.06<DCAT<0.08");
    h3->SetBinLabel(5,"0.08<DCAT<0.1");
    h3->SetBinLabel(6,"0.1<DCAT<0.2");
    h3->SetBinLabel(7,"0.2<DCAT<0.3");
    h3->SetBinLabel(8,"0.3<DCAT<0.4");
    h3->SetBinLabel(9,"0.4<DCAT<0.5");

    auto parcut2 = histos.get<TH1>(HIST("TrueClosest_Dcatcut_Beauty"));
    auto* h4 = parcut2->GetXaxis();
    h4->SetBinLabel(1,"0<DCAT<0.02");
    h4->SetBinLabel(2,"0.02<DCAT<0.04");
    h4->SetBinLabel(3,"0.04<DCAT<0.06");
    h4->SetBinLabel(4,"0.06<DCAT<0.08");
    h4->SetBinLabel(5,"0.08<DCAT<0.1");
    h4->SetBinLabel(6,"0.1<DCAT<0.2");
    h4->SetBinLabel(7,"0.2<DCAT<0.3");
    h4->SetBinLabel(8,"0.3<DCAT<0.4");
    h4->SetBinLabel(9,"0.4<DCAT<0.5");

    auto parcut_f2 = histos.get<TH1>(HIST("FakeClosest_Dcatcut_Beauty"));
    auto* h5 = parcut_f2->GetXaxis();
    h5->SetBinLabel(1,"0<DCAT<0.02");
    h5->SetBinLabel(2,"0.02<DCAT<0.04");
    h5->SetBinLabel(3,"0.04<DCAT<0.06");
    h5->SetBinLabel(4,"0.06<DCAT<0.08");
    h5->SetBinLabel(5,"0.08<DCAT<0.1");
    h5->SetBinLabel(6,"0.1<DCAT<0.2");
    h5->SetBinLabel(7,"0.2<DCAT<0.3");
    h5->SetBinLabel(8,"0.3<DCAT<0.4");
    h5->SetBinLabel(9,"0.4<DCAT<0.5");

    auto parcut3 = histos.get<TH1>(HIST("Signal_Dcatcut_npCharm"));
    auto* h6 = parcut3->GetXaxis();
    h6->SetBinLabel(1,"0<DCAT<0.02");
    h6->SetBinLabel(2,"0.02<DCAT<0.04");
    h6->SetBinLabel(3,"0.04<DCAT<0.06");
    h6->SetBinLabel(4,"0.06<DCAT<0.08");
    h6->SetBinLabel(5,"0.08<DCAT<0.1");
    h6->SetBinLabel(6,"0.1<DCAT<0.2");
    h6->SetBinLabel(7,"0.2<DCAT<0.3");
    h6->SetBinLabel(8,"0.3<DCAT<0.4");
    h6->SetBinLabel(9,"0.4<DCAT<0.5");

    auto parcut4 = histos.get<TH1>(HIST("TrueClosest_Dcatcut_npCharm"));
    auto* h7 = parcut4->GetXaxis();
    h7->SetBinLabel(1,"0<DCAT<0.02");
    h7->SetBinLabel(2,"0.02<DCAT<0.04");
    h7->SetBinLabel(3,"0.04<DCAT<0.06");
    h7->SetBinLabel(4,"0.06<DCAT<0.08");
    h7->SetBinLabel(5,"0.08<DCAT<0.1");
    h7->SetBinLabel(6,"0.1<DCAT<0.2");
    h7->SetBinLabel(7,"0.2<DCAT<0.3");
    h7->SetBinLabel(8,"0.3<DCAT<0.4");
    h7->SetBinLabel(9,"0.4<DCAT<0.5");

    auto parcut_f3 = histos.get<TH1>(HIST("FakeClosest_Dcatcut_npCharm"));
    auto* h8 = parcut_f3->GetXaxis();
    h8->SetBinLabel(1,"0<DCAT<0.02");
    h8->SetBinLabel(2,"0.02<DCAT<0.04");
    h8->SetBinLabel(3,"0.04<DCAT<0.06");
    h8->SetBinLabel(4,"0.06<DCAT<0.08");
    h8->SetBinLabel(5,"0.08<DCAT<0.1");
    h8->SetBinLabel(6,"0.1<DCAT<0.2");
    h8->SetBinLabel(7,"0.2<DCAT<0.3");
    h8->SetBinLabel(8,"0.3<DCAT<0.4");
    h8->SetBinLabel(9,"0.4<DCAT<0.5");

  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
							 soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
							 soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
							 aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
							 aod::McCollisions const&,
               aod::MCPair const& truepairs,
               aod::MuInfo const& muoninfos)
  {
    const double k_mass = 0.493677;
    const double mu_mass = 0.105658;
    if(collision.has_mcCollision()){
      if(fabs(collision.posZ())<10){
        vector<uint64_t> ambTrackIds;
        ambTrackIds.clear();
        for(auto& amfttrack : amfttracks){
          ambTrackIds.push_back(amfttrack.mfttrackId());
        }
        
        for(auto& truepair : truepairs){
          if(truepair.colid()!=collision.globalIndex()) continue;
          if(truepair.colx()!=collision.posX() || truepair.coly()!=collision.posY() || truepair.colz()!=collision.posZ()) continue;
          double Col_x = collision.posX();
          double Col_y = collision.posY();
          double Col_z = collision.posZ();

          vector<uint64_t> fwdTrackIds;
          fwdTrackIds.clear();
          for(auto& fwdtrack : fwdtracks){
            if(!fwdtrack.has_collision() || fwdtrack.collisionId()!=collision.globalIndex()) continue;
            if(!fwdtrack.has_mcParticle()) continue;
            fwdTrackIds.push_back(fwdtrack.mcParticleId());
            if(fwdtrack.mcParticleId()!=truepair.muid()){
              SMatrix5 fwdpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
              vector<double> fwdv2;
              SMatrix55 fwdcov2(fwdv2.begin(), fwdv2.end());
              o2::track::TrackParCovFwd fwdpars2{fwdtrack.z(), fwdpars, fwdcov2, fwdtrack.chi2()};
              fwdpars2.propagateToZlinear(Col_z);
              double fwddcaX = fwdpars2.getX();
              double fwddcaY = fwdpars2.getY();
              double fwdmftX = fwdtrack.x();
              double fwdmftY = fwdtrack.y();
              double fwdmftZ = fwdtrack.z();
              double fwd_px = fwdtrack.px();
              double fwd_py = fwdtrack.py();
              double fwd_pz = fwdtrack.pz();
              double fwd_p = fwdtrack.p();
              double fwddcat = sqrt(pow(fwddcaX-Col_x,2)+pow(fwddcaY-Col_y,2));
              histos.fill(HIST("BG_DCAT_all"), fwdtrack.pt());
              if(fwddcat>0.01 && fwddcat<0.1){
                histos.fill(HIST("BG_DCAT_cut"), fwdtrack.pt());
              }
              vector<vector<double>> data_r_fwd;

              for(auto& mfttrack : mfttracks){
                if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
                if(mfttrack.collisionId()==truepair.colid() && mfttrack.mcParticleId()!=truepair.muid()){
                  if(find(fwdTrackIds.begin(), fwdTrackIds.end(), mfttrack.mcParticleId()) != fwdTrackIds.end()) continue;

                  auto mcParticle_mft = mfttrack.mcParticle();
                  if(fabs(mcParticle_mft.pdgCode())==2212) continue;
                  //if(truepair.pairid_1()!=mcParticle_mft.globalIndex() && find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack.globalIndex())!=ambTrackIds.end()) continue;

                  double mftchi2 = mfttrack.chi2();
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

                  double canpx = mfttrack.px();
                  double canpy = mfttrack.py();
                  double canpz = mfttrack.pz();
                  double canp = mfttrack.p();

                  auto unit_Na = sqrt(pow(fwdmftX-fwddcaX,2)+pow(fwdmftY-fwddcaY,2)+pow(fwdmftZ-Col_z,2));
                  auto unit_Nc = sqrt(pow(canmftX-candcaX,2)+pow(canmftY-candcaY,2)+pow(canmftZ-Col_z,2));
                  auto Nax = (fwdmftX-fwddcaX)/unit_Na;
                  auto Nay = (fwdmftY-fwddcaY)/unit_Na;
                  auto Naz = (fwdmftZ-Col_z)/unit_Na;
                  auto Ncx = (canmftX-candcaX)/unit_Nc;
                  auto Ncy = (canmftY-candcaY)/unit_Nc;
                  auto Ncz = (canmftZ-Col_z)/unit_Nc;
                  auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
                  auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
                  auto A3 = (fwddcaX-candcaX)*Nax + (fwddcaY-candcaY)*Nay + (Col_z-Col_z)*Naz;
                  auto B1 = A2;
                  auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
                  auto B3 = (candcaX-fwddcaX)*Ncx + (candcaY-fwddcaY)*Ncy + (Col_z-Col_z)*Ncz;
                  auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
                  auto s = -((A2*t+A3)/A1);

                  double predict_mux = fwddcaX + s*Nax;
                  double predict_muy = fwddcaY + s*Nay;
                  double predict_muz = Col_z + s*Naz;
                  double predict_canx = candcaX + t*Ncx;
                  double predict_cany = candcaY + t*Ncy;
                  double predict_canz = Col_z + t*Ncz;
                  double r_xyz = sqrt(pow(predict_canx-predict_mux,2) + pow(predict_cany-predict_muy,2) + pow(predict_canz-predict_muz,2));

                  auto vecx_mu = fwdmftX - fwddcaX;
                  auto vecy_mu = fwdmftY - fwddcaY;
                  auto vecz_mu = fwdmftZ - Col_z;
                  auto vecx_can = canmftX - candcaX;
                  auto vecy_can = canmftY - candcaY;
                  auto vecz_can = canmftZ - Col_z;
                  auto cosxy = (vecx_mu*vecx_can + vecy_mu*vecy_can + vecz_mu*vecz_can)/((sqrt(pow(vecx_mu,2)+pow(vecy_mu,2)+pow(vecz_mu,2)))*(sqrt(pow(vecx_can,2)+pow(vecy_can,2)+pow(vecz_can,2))));
                  
                  auto pcaX = ((predict_mux+predict_canx)/2)-Col_x;
                  auto pcaY = ((predict_muy+predict_cany)/2)-Col_y;
                  auto pcaZ = ((predict_muz+predict_canz)/2)-Col_z;
                  double dcat = sqrt(pow(fwddcaX-Col_x,2)+pow(fwddcaY-Col_y,2));
                  //double dcat = sqrt(pow(candcaX-Col_x,2)+pow(candcaY-Col_y,2));
                  double dcaxy_diff= sqrt(pow(fwddcaX-candcaX,2)+pow(fwddcaY-candcaY,2));

                  double mcid = mfttrack.mcParticleId();
                  double out_colid = mfttrack.collisionId();

                  double add_dist = 0;

                  double mdx_can = candcaX - canmftX;
                  double mdy_can = candcaY - canmftY;
                  bool PairVec = true;
                  if((canmftX-Col_x)>0 && (canmftY-Col_y)>0){
                    if(mdx_can>0 || mdy_can>0) PairVec=false;
                  }else if((canmftX-Col_x)>0 && (canmftY-Col_y)<0){
                    if(mdx_can>0 || mdy_can<0) PairVec=false;
                  }else if((canmftX-Col_x)<0 && (canmftY-Col_y)<0){
                    if(mdx_can<0 || mdy_can<0) PairVec=false;
                  }else if((canmftX-Col_x)<0 && (canmftY-Col_y)>0){
                    if(mdx_can<0 || mdy_can>0) PairVec=false;
                  }
                  int pair_vec = 1;
                  if(PairVec==false){
                    pair_vec = -1;
                  }
                  double vec_cor = truepair.muvec() - pair_vec;

                  double ans = 0;
                  if(mcid==truepair.pairid_1()) ans=1;

                  double inv_mass_candidate = sqrt(pow(sqrt(pow(mu_mass,2)+pow(fwd_p,2))+sqrt(pow(k_mass,2)+pow(canp,2)), 2) - (pow(fwd_px+canpx,2)+pow(fwd_py+canpy,2)+pow(fwd_pz+canpz,2)));

                  vector<double> pardata_r = {r_xyz, cosxy, mcid, fabs(mcParticle_mft.pdgCode()), pcaX, pcaY, pcaZ, add_dist, dcat, dcaxy_diff, out_colid, truepair.pt(), vec_cor, ans, inv_mass_candidate};
                  data_r_fwd.push_back(pardata_r);
                }
              }
              sort(data_r_fwd.begin(), data_r_fwd.end());
              if(data_r_fwd.at(0).at(8)<0.1) histos.fill(HIST("BG_InvMass"), data_r_fwd.at(0).at(14));
              /*if(truepair.particletype()==3){
                int counter_cut=0;
                for(int i=0; i<data_r_fwd.size(); i++){
                  if(data_r_fwd.at(i).at(0)>0.05){
                    cout << endl;
                    break;
                  }
                  if(data_r_fwd.at(i).at(1)>0.99){
                    if(data_r_fwd.at(i).at(7)<30){
                      if(data_r_fwd.at(i).at(8)>0.01 && data_r_fwd.at(i).at(8)<0.1){
                        if(data_r_fwd.at(i).at(9)<0.2 && data_r_fwd.at(i).at(9)>0.02){
                          if(data_r_fwd.at(i).at(14)<2.2){
                            if(counter_cut==0){
                              cout << "BG_Data = { " << data_r_fwd.at(i).at(0) << ", " << data_r_fwd.at(i).at(1) << ", " << data_r_fwd.at(i).at(7) << ", " << data_r_fwd.at(i).at(8) << ", " << data_r_fwd.at(i).at(9) << ", " << data_r_fwd.at(i).at(14) << ", " << data_r_fwd.at(i).at(3) << ", " << data_r_fwd.at(i).at(13) << " }";
                              counter_cut++;
                              if(data_r_fwd.at(i).at(13)==1){
                                histos.fill(HIST("True_first"), data_r_fwd.size());
                              }else{
                                histos.fill(HIST("Fake_first"), data_r_fwd.size());
                              }
                            }else{
                              cout << "," << endl << "{ " << data_r_fwd.at(i).at(0) << ", " << data_r_fwd.at(i).at(1) << ", " << data_r_fwd.at(i).at(7) << ", " << data_r_fwd.at(i).at(8) << ", " << data_r_fwd.at(i).at(9) << ", " << data_r_fwd.at(i).at(14) << ", " << data_r_fwd.at(i).at(3) << ", " << data_r_fwd.at(i).at(13) << " }";
                            }
                            if(data_r_fwd.at(i).at(13)==1){
                              histos.fill(HIST("HasPair_wincut"), i+1);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }*/

              histos.fill(HIST("BG_pcar"), data_r_fwd.at(0).at(0));
              histos.fill(HIST("BG_c_PCAD"), sqrt(pow(data_r_fwd.at(0).at(4),2)+pow(data_r_fwd.at(0).at(5),2)+pow(data_r_fwd.at(0).at(6),2)));
              for(int i=0; i<data_r_fwd.size(); i++){
                auto pcad_fwd_all = sqrt(pow(data_r_fwd.at(i).at(4),2)+pow(data_r_fwd.at(i).at(5),2)+pow(data_r_fwd.at(i).at(6),2));
                if(pcad_fwd_all<1 && data_r_fwd.at(i).at(8)>0.01 && data_r_fwd.at(i).at(8)<0.1){
                  histos.fill(HIST("BG_DCAT_PCAD_cut"), fwdtrack.pt());
                  break;
                }
              }
            }
          }

          SMatrix5 mupars(truepair.x(), truepair.y(), truepair.phi(), truepair.tgl(), truepair.signed1Pt());
          vector<double> muv2;
          SMatrix55 mucovs2(muv2.begin(), muv2.end());
          o2::track::TrackParCovFwd mupars2{truepair.z(), mupars, mucovs2, truepair.chi2()};
          mupars2.propagateToZlinear(Col_z);
          double mudcaX = mupars2.getX();
          double mudcaY = mupars2.getY();
          double mumftX = truepair.x();
          double mumftY = truepair.y();
          double mumftZ = truepair.z();

          //Only DCAT
          double mudcat = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
          if(truepair.particletype()==1){
            if(mudcat>0 && mudcat<0.02){
              histos.fill(HIST("Signal_Dcatcut_Beauty"), 0);
            }else if(mudcat>=0.02 && mudcat<0.04){
              histos.fill(HIST("Signal_Dcatcut_Beauty"), 1);
            }else if(mudcat>=0.04 && mudcat<0.06){
              histos.fill(HIST("Signal_Dcatcut_Beauty"), 2);
            }else if(mudcat>=0.06 && mudcat<0.08){
              histos.fill(HIST("Signal_Dcatcut_Beauty"), 3);
            }else if(mudcat>=0.08 && mudcat<0.1){
              histos.fill(HIST("Signal_Dcatcut_Beauty"), 4);
            }else if(mudcat>=0.1 && mudcat<0.2){
              histos.fill(HIST("Signal_Dcatcut_Beauty"), 5);
            }else if(mudcat>=0.2 && mudcat<0.3){
              histos.fill(HIST("Signal_Dcatcut_Beauty"), 6);
            }else if(mudcat>=0.3 && mudcat<0.4){
              histos.fill(HIST("Signal_Dcatcut_Beauty"), 7);
            }else if(mudcat>=0.4 && mudcat<0.5){
              histos.fill(HIST("Signal_Dcatcut_Beauty"), 8);
            }
          }
          if(truepair.particletype()==2){
            if(mudcat>0 && mudcat<0.02){
              histos.fill(HIST("Signal_Dcatcut_npCharm"), 0);
            }else if(mudcat>=0.02 && mudcat<0.04){
              histos.fill(HIST("Signal_Dcatcut_npCharm"), 1);
            }else if(mudcat>=0.04 && mudcat<0.06){
              histos.fill(HIST("Signal_Dcatcut_npCharm"), 2);
            }else if(mudcat>=0.06 && mudcat<0.08){
              histos.fill(HIST("Signal_Dcatcut_npCharm"), 3);
            }else if(mudcat>=0.08 && mudcat<0.1){
              histos.fill(HIST("Signal_Dcatcut_npCharm"), 4);
            }else if(mudcat>=0.1 && mudcat<0.2){
              histos.fill(HIST("Signal_Dcatcut_npCharm"), 5);
            }else if(mudcat>=0.2 && mudcat<0.3){
              histos.fill(HIST("Signal_Dcatcut_npCharm"), 6);
            }else if(mudcat>=0.3 && mudcat<0.4){
              histos.fill(HIST("Signal_Dcatcut_npCharm"), 7);
            }else if(mudcat>=0.4 && mudcat<0.5){
              histos.fill(HIST("Signal_Dcatcut_npCharm"), 8);
            }
          }
          if(truepair.particletype()==3){
            histos.fill(HIST("c_DCAT_all"), truepair.pt());

            if(mudcat>0.01 && mudcat<0.1){
              histos.fill(HIST("c_DCAT_cut"), truepair.pt());
            }

            if(mudcat>0 && mudcat<0.02){
              histos.fill(HIST("Signal_Dcatcut_pCharm"), 0);
            }else if(mudcat>=0.02 && mudcat<0.04){
              histos.fill(HIST("Signal_Dcatcut_pCharm"), 1);
            }else if(mudcat>=0.04 && mudcat<0.06){
              histos.fill(HIST("Signal_Dcatcut_pCharm"), 2);
            }else if(mudcat>=0.06 && mudcat<0.08){
              histos.fill(HIST("Signal_Dcatcut_pCharm"), 3);
            }else if(mudcat>=0.08 && mudcat<0.1){
              histos.fill(HIST("Signal_Dcatcut_pCharm"), 4);
            }else if(mudcat>=0.1 && mudcat<0.2){
              histos.fill(HIST("Signal_Dcatcut_pCharm"), 5);
            }else if(mudcat>=0.2 && mudcat<0.3){
              histos.fill(HIST("Signal_Dcatcut_pCharm"), 6);
            }else if(mudcat>=0.3 && mudcat<0.4){
              histos.fill(HIST("Signal_Dcatcut_pCharm"), 7);
            }else if(mudcat>=0.4 && mudcat<0.5){
              histos.fill(HIST("Signal_Dcatcut_pCharm"), 8);
            }
          }

          vector<vector<double>> data_r;

          for(auto& mfttrack : mfttracks){
            if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
            if(mfttrack.collisionId()==truepair.colid() && mfttrack.mcParticleId()!=truepair.muid()){
              if(find(fwdTrackIds.begin(), fwdTrackIds.end(), mfttrack.mcParticleId()) != fwdTrackIds.end()) continue;

              auto mcParticle_mft = mfttrack.mcParticle();
              if(fabs(mcParticle_mft.pdgCode())==2212) continue;
              //if(truepair.pairid_1()!=mcParticle_mft.globalIndex() && find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack.globalIndex())!=ambTrackIds.end()) continue;

              double mftchi2 = mfttrack.chi2();
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

              double canpx = mfttrack.px();
              double canpy = mfttrack.py();
              double canpz = mfttrack.pz();
              double canp = mfttrack.p();

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
              double dcat = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
              //double dcat = sqrt(pow(candcaX-Col_x,2)+pow(candcaY-Col_y,2));
              double dcaxy_diff= sqrt(pow(mudcaX-candcaX,2)+pow(mudcaY-candcaY,2));

              double mcid = mfttrack.mcParticleId();
              double out_colid = mfttrack.collisionId();

              double add_dist = 0;
              for(int i=0; i<5000; i++){
                auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                if((cut_z-Col_z)>0) break;
                mupars2.propagateToZlinear(cut_z);
                mftpars2.propagateToZlinear(cut_z);
                auto diff = sqrt(pow(mftpars2.getX()-mupars2.getX(),2)+pow(mftpars2.getY()-mupars2.getY(),2));
                add_dist += diff;
              }

              double mdx_can = candcaX - canmftX;
              double mdy_can = candcaY - canmftY;
              bool PairVec = true;
              if((canmftX-Col_x)>0 && (canmftY-Col_y)>0){
                if(mdx_can>0 || mdy_can>0) PairVec=false;
              }else if((canmftX-Col_x)>0 && (canmftY-Col_y)<0){
                if(mdx_can>0 || mdy_can<0) PairVec=false;
              }else if((canmftX-Col_x)<0 && (canmftY-Col_y)<0){
                if(mdx_can<0 || mdy_can<0) PairVec=false;
              }else if((canmftX-Col_x)<0 && (canmftY-Col_y)>0){
                if(mdx_can<0 || mdy_can>0) PairVec=false;
              }
              int pair_vec = 1;
              if(PairVec==false){
                pair_vec = -1;
              }
              double vec_cor = truepair.muvec() - pair_vec;

              double ans = 0;
              if(mcid==truepair.pairid_1()){
                ans=1;
                histos.fill(HIST("Add_dist_true"), add_dist);
              }else if(mcid!=truepair.pairid_1()){
                histos.fill(HIST("CosSim_fake"), cosxy);
                histos.fill(HIST("Add_dist_fake"), add_dist);
              }

              double inv_mass_candidate = sqrt(pow(sqrt(pow(mu_mass,2)+pow(truepair.p(),2))+sqrt(pow(k_mass,2)+pow(canp,2)), 2) - (pow(truepair.px()+canpx,2)+pow(truepair.py()+canpy,2)+pow(truepair.pz()+canpz,2)));

              vector<double> pardata_r = {r_xyz, cosxy, mcid, fabs(mcParticle_mft.pdgCode()), pcaX, pcaY, pcaZ, add_dist, dcat, dcaxy_diff, out_colid, truepair.pt(), vec_cor, ans, inv_mass_candidate};
              data_r.push_back(pardata_r);
            }
          }
          sort(data_r.begin(), data_r.end());
          double pcad = sqrt(pow(data_r.at(0).at(4),2)+pow(data_r.at(0).at(5),2)+pow(data_r.at(0).at(6),2));

          if(truepair.particletype()==3){
            for(int i=0; i<data_r.size(); i++){
              if(i>5) break;
              int j = i;
              if(data_r.at(i).at(13)==1){
                while(j<5){
                  histos.fill(HIST("RankEfficiency_pCharm"), j+1);
                  j++;
                }
              }
            }
          }
          if(truepair.particletype()==2){
            for(int i=0; i<data_r.size(); i++){
              if(i>5) break;
              int j = i;
              if(data_r.at(i).at(13)==1){
                while(j<5){
                  histos.fill(HIST("RankEfficiency_npCharm"), j+1);
                  j++;
                }
              }
            }
          }

          if(truepair.particletype()==3){
            for(int i=0; i<data_r.size(); i++){
              auto pcad_all = sqrt(pow(data_r.at(i).at(4),2)+pow(data_r.at(i).at(5),2)+pow(data_r.at(i).at(6),2));
              if(data_r.at(i).at(13)==1){
                histos.fill(HIST("AllTrueMatch_pCharm_pcad"), pcad_all);
              }else if(data_r.at(i).at(13)==0){
                histos.fill(HIST("AllFakeMatch_pCharm_pcad"), pcad_all);
              }
            }
          
          }else if(truepair.particletype()==1){
            for(int i=0; i<data_r.size(); i++){
              auto pcad_all = sqrt(pow(data_r.at(i).at(4),2)+pow(data_r.at(i).at(5),2)+pow(data_r.at(i).at(6),2));
              if(data_r.at(i).at(13)==1){
                histos.fill(HIST("AllTrueMatch_Beauty_pcad"), pcad_all);
              }else if(data_r.at(i).at(13)==0){
                histos.fill(HIST("AllFakeMatch_Beauty_pcad"), pcad_all);
              }
            }
          }else if(truepair.particletype()==2){
            for(int i=0; i<data_r.size(); i++){
              auto pcad_all = sqrt(pow(data_r.at(i).at(4),2)+pow(data_r.at(i).at(5),2)+pow(data_r.at(i).at(6),2));
              if(data_r.at(i).at(13)==1){
                histos.fill(HIST("AllTrueMatch_npCharm_pcad"), pcad_all);
              }else if(data_r.at(i).at(13)==0){
                histos.fill(HIST("AllFakeMatch_npCharm_pcad"), pcad_all);
              }
            }
          }
          
          // Use PCA Info.

          //1109~
          if(truepair.particletype()==3){
            int counter_cut=0;
            for(int i=0; i<data_r.size(); i++){
              if(data_r.at(i).at(0)>0.05){
                cout << endl;
                break;
              }
              if(data_r.at(i).at(1)>0.99){
                if(data_r.at(i).at(7)<30){
                  if(data_r.at(i).at(8)>0.01 && data_r.at(i).at(8)<0.1){
                    if(data_r.at(i).at(9)<0.2 && data_r.at(i).at(9)>0.02){
                      if(data_r.at(i).at(14)<2.2){
                        if(counter_cut==0){
                          cout << "Data = { " << data_r.at(i).at(0) << ", " << data_r.at(i).at(1) << ", " << data_r.at(i).at(7) << ", " << data_r.at(i).at(8) << ", " << data_r.at(i).at(9) << ", " << data_r.at(i).at(14) << ", " << data_r.at(i).at(3) << ", " << data_r.at(i).at(13) << " }";
                          counter_cut++;
                          if(data_r.at(i).at(13)==1){
                            histos.fill(HIST("True_first"), data_r.size());
                          }else{
                            histos.fill(HIST("Fake_first"), data_r.size());
                          }
                        }else{
                          cout << "," << endl << "{ " << data_r.at(i).at(0) << ", " << data_r.at(i).at(1) << ", " << data_r.at(i).at(7) << ", " << data_r.at(i).at(8) << ", " << data_r.at(i).at(9) << ", " << data_r.at(i).at(14) << ", " << data_r.at(i).at(3) << ", " << data_r.at(i).at(13) << " }";
                        }
                        if(data_r.at(i).at(13)==1){
                          histos.fill(HIST("HasPair_wincut"), i+1);
                        }
                      }
                    }
                  }
                }
              }
            }
          }


          //~1109
          if(truepair.particletype()==3){
            for(int i=0; i<data_r.size(); i++){
              auto pcad_all = sqrt(pow(data_r.at(i).at(4),2)+pow(data_r.at(i).at(5),2)+pow(data_r.at(i).at(6),2));
              if(pcad_all<1 && data_r.at(i).at(8)>0.01 && data_r.at(i).at(8)<0.1){
                histos.fill(HIST("c_DCAT_PCAD_cut"), truepair.pt());
                break;
              }
            }

            for(int i=0; i<data_r.size(); i++){
              if(i==5) break;
              auto pcad_all = sqrt(pow(data_r.at(i).at(4),2)+pow(data_r.at(i).at(5),2)+pow(data_r.at(i).at(6),2));
              if(pcad_all<1 && data_r.at(i).at(13)==1 && data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.5){
                if(data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.02){
                  histos.fill(HIST("TrueClosest_Dcatcut_pCharm"), 0);
                }else if(data_r.at(i).at(8)>=0.02 && data_r.at(i).at(8)<0.04){
                  histos.fill(HIST("TrueClosest_Dcatcut_pCharm"), 1);
                }else if(data_r.at(i).at(8)>=0.04 && data_r.at(i).at(8)<0.06){
                  histos.fill(HIST("TrueClosest_Dcatcut_pCharm"), 2);
                }else if(data_r.at(i).at(8)>=0.06 && data_r.at(i).at(8)<0.08){
                  histos.fill(HIST("TrueClosest_Dcatcut_pCharm"), 3);
                }else if(data_r.at(i).at(8)>=0.08 && data_r.at(i).at(8)<0.1){
                  histos.fill(HIST("TrueClosest_Dcatcut_pCharm"), 4);
                }else if(data_r.at(i).at(8)>=0.1 && data_r.at(i).at(8)<0.2){
                  histos.fill(HIST("TrueClosest_Dcatcut_pCharm"), 5);
                }else if(data_r.at(i).at(8)>=0.2 && data_r.at(i).at(8)<0.3){
                  histos.fill(HIST("TrueClosest_Dcatcut_pCharm"), 6);
                }else if(data_r.at(i).at(8)>=0.3 && data_r.at(i).at(8)<0.4){
                  histos.fill(HIST("TrueClosest_Dcatcut_pCharm"), 7);
                }else if(data_r.at(i).at(8)>=0.4 && data_r.at(i).at(8)<0.5){
                  histos.fill(HIST("TrueClosest_Dcatcut_pCharm"), 8);
                }
                break;
              }else if(pcad_all<1 && data_r.at(i).at(13)==0 && data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.5){
                if(data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.02){
                  histos.fill(HIST("FakeClosest_Dcatcut_pCharm"), 0);
                }else if(data_r.at(i).at(8)>=0.02 && data_r.at(i).at(8)<0.04){
                  histos.fill(HIST("FakeClosest_Dcatcut_pCharm"), 1);
                }else if(data_r.at(i).at(8)>=0.04 && data_r.at(i).at(8)<0.06){
                  histos.fill(HIST("FakeClosest_Dcatcut_pCharm"), 2);
                }else if(data_r.at(i).at(8)>=0.06 && data_r.at(i).at(8)<0.08){
                  histos.fill(HIST("FakeClosest_Dcatcut_pCharm"), 3);
                }else if(data_r.at(i).at(8)>=0.08 && data_r.at(i).at(8)<0.1){
                  histos.fill(HIST("FakeClosest_Dcatcut_pCharm"), 4);
                }else if(data_r.at(i).at(8)>=0.1 && data_r.at(i).at(8)<0.2){
                  histos.fill(HIST("FakeClosest_Dcatcut_pCharm"), 5);
                }else if(data_r.at(i).at(8)>=0.2 && data_r.at(i).at(8)<0.3){
                  histos.fill(HIST("FakeClosest_Dcatcut_pCharm"), 6);
                }else if(data_r.at(i).at(8)>=0.3 && data_r.at(i).at(8)<0.4){
                  histos.fill(HIST("FakeClosest_Dcatcut_pCharm"), 7);
                }else if(data_r.at(i).at(8)>=0.4 && data_r.at(i).at(8)<0.5){
                  histos.fill(HIST("FakeClosest_Dcatcut_pCharm"), 8);
                }
                break;
              }
            } 
          }else if(truepair.particletype()==2){
            for(int i=0; i<data_r.size(); i++){
              if(i==5) break;
              auto pcad_all = sqrt(pow(data_r.at(i).at(4),2)+pow(data_r.at(i).at(5),2)+pow(data_r.at(i).at(6),2));
              if(pcad_all<1 && data_r.at(i).at(13)==1 && data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.5){
                if(data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.02){
                  histos.fill(HIST("TrueClosest_Dcatcut_npCharm"), 0);
                }else if(data_r.at(i).at(8)>=0.02 && data_r.at(i).at(8)<0.04){
                  histos.fill(HIST("TrueClosest_Dcatcut_npCharm"), 1);
                }else if(data_r.at(i).at(8)>=0.04 && data_r.at(i).at(8)<0.06){
                  histos.fill(HIST("TrueClosest_Dcatcut_npCharm"), 2);
                }else if(data_r.at(i).at(8)>=0.06 && data_r.at(i).at(8)<0.08){
                  histos.fill(HIST("TrueClosest_Dcatcut_npCharm"), 3);
                }else if(data_r.at(i).at(8)>=0.08 && data_r.at(i).at(8)<0.1){
                  histos.fill(HIST("TrueClosest_Dcatcut_npCharm"), 4);
                }else if(data_r.at(i).at(8)>=0.1 && data_r.at(i).at(8)<0.2){
                  histos.fill(HIST("TrueClosest_Dcatcut_npCharm"), 5);
                }else if(data_r.at(i).at(8)>=0.2 && data_r.at(i).at(8)<0.3){
                  histos.fill(HIST("TrueClosest_Dcatcut_npCharm"), 6);
                }else if(data_r.at(i).at(8)>=0.3 && data_r.at(i).at(8)<0.4){
                  histos.fill(HIST("TrueClosest_Dcatcut_npCharm"), 7);
                }else if(data_r.at(i).at(8)>=0.4 && data_r.at(i).at(8)<0.5){
                  histos.fill(HIST("TrueClosest_Dcatcut_npCharm"), 8);
                }
                break;
              }else if(pcad_all<1 && data_r.at(i).at(13)==0){
                if(data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.02 && data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.5){
                  histos.fill(HIST("FakeClosest_Dcatcut_npCharm"), 0);
                }else if(data_r.at(i).at(8)>=0.02 && data_r.at(i).at(8)<0.04){
                  histos.fill(HIST("FakeClosest_Dcatcut_npCharm"), 1);
                }else if(data_r.at(i).at(8)>=0.04 && data_r.at(i).at(8)<0.06){
                  histos.fill(HIST("FakeClosest_Dcatcut_npCharm"), 2);
                }else if(data_r.at(i).at(8)>=0.06 && data_r.at(i).at(8)<0.08){
                  histos.fill(HIST("FakeClosest_Dcatcut_npCharm"), 3);
                }else if(data_r.at(i).at(8)>=0.08 && data_r.at(i).at(8)<0.1){
                  histos.fill(HIST("FakeClosest_Dcatcut_npCharm"), 4);
                }else if(data_r.at(i).at(8)>=0.1 && data_r.at(i).at(8)<0.2){
                  histos.fill(HIST("FakeClosest_Dcatcut_npCharm"), 5);
                }else if(data_r.at(i).at(8)>=0.2 && data_r.at(i).at(8)<0.3){
                  histos.fill(HIST("FakeClosest_Dcatcut_npCharm"), 6);
                }else if(data_r.at(i).at(8)>=0.3 && data_r.at(i).at(8)<0.4){
                  histos.fill(HIST("FakeClosest_Dcatcut_npCharm"), 7);
                }else if(data_r.at(i).at(8)>=0.4 && data_r.at(i).at(8)<0.5){
                  histos.fill(HIST("FakeClosest_Dcatcut_npCharm"), 8);
                }
                break;
              }
            }
          }else if(truepair.particletype()==1){
            for(int i=0; i<data_r.size(); i++){
              if(i==5) break;
              auto pcad_all = sqrt(pow(data_r.at(i).at(4),2)+pow(data_r.at(i).at(5),2)+pow(data_r.at(i).at(6),2));
              if(pcad_all<1 && data_r.at(i).at(13)==1 && data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.5){
                if(data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.02){
                  histos.fill(HIST("TrueClosest_Dcatcut_Beauty"), 0);
                }else if(data_r.at(i).at(8)>=0.02 && data_r.at(i).at(8)<0.04){
                  histos.fill(HIST("TrueClosest_Dcatcut_Beauty"), 1);
                }else if(data_r.at(i).at(8)>=0.04 && data_r.at(i).at(8)<0.06){
                  histos.fill(HIST("TrueClosest_Dcatcut_Beauty"), 2);
                }else if(data_r.at(i).at(8)>=0.06 && data_r.at(i).at(8)<0.08){
                  histos.fill(HIST("TrueClosest_Dcatcut_Beauty"), 3);
                }else if(data_r.at(i).at(8)>=0.08 && data_r.at(i).at(8)<0.1){
                  histos.fill(HIST("TrueClosest_Dcatcut_Beauty"), 4);
                }else if(data_r.at(i).at(8)>=0.1 && data_r.at(i).at(8)<0.2){
                  histos.fill(HIST("TrueClosest_Dcatcut_Beauty"), 5);
                }else if(data_r.at(i).at(8)>=0.2 && data_r.at(i).at(8)<0.3){
                  histos.fill(HIST("TrueClosest_Dcatcut_Beauty"), 6);
                }else if(data_r.at(i).at(8)>=0.3 && data_r.at(i).at(8)<0.4){
                  histos.fill(HIST("TrueClosest_Dcatcut_Beauty"), 7);
                }else if(data_r.at(i).at(8)>=0.4 && data_r.at(i).at(8)<0.5){
                  histos.fill(HIST("TrueClosest_Dcatcut_Beauty"), 8);
                }
                break;
              }else if(pcad_all<1 && data_r.at(i).at(13)==0 && data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.5){
                if(data_r.at(i).at(8)>0 && data_r.at(i).at(8)<0.02){
                  histos.fill(HIST("FakeClosest_Dcatcut_Beauty"), 0);
                }else if(data_r.at(i).at(8)>=0.02 && data_r.at(i).at(8)<0.04){
                  histos.fill(HIST("FakeClosest_Dcatcut_Beauty"), 1);
                }else if(data_r.at(i).at(8)>=0.04 && data_r.at(i).at(8)<0.06){
                  histos.fill(HIST("FakeClosest_Dcatcut_Beauty"), 2);
                }else if(data_r.at(i).at(8)>=0.06 && data_r.at(i).at(8)<0.08){
                  histos.fill(HIST("FakeClosest_Dcatcut_Beauty"), 3);
                }else if(data_r.at(i).at(8)>=0.08 && data_r.at(i).at(8)<0.1){
                  histos.fill(HIST("FakeClosest_Dcatcut_Beauty"), 4);
                }else if(data_r.at(i).at(8)>=0.1 && data_r.at(i).at(8)<0.2){
                  histos.fill(HIST("FakeClosest_Dcatcut_Beauty"), 5);
                }else if(data_r.at(i).at(8)>=0.2 && data_r.at(i).at(8)<0.3){
                  histos.fill(HIST("FakeClosest_Dcatcut_Beauty"), 6);
                }else if(data_r.at(i).at(8)>=0.3 && data_r.at(i).at(8)<0.4){
                  histos.fill(HIST("FakeClosest_Dcatcut_Beauty"), 7);
                }else if(data_r.at(i).at(8)>=0.4 && data_r.at(i).at(8)<0.5){
                  histos.fill(HIST("FakeClosest_Dcatcut_Beauty"), 8);
                }
                break;
              }
            }
          }

          if(data_r.at(0).at(13)==1){
            if(truepair.particletype()==3){
              histos.fill(HIST("TrueMatch_pCharm_pcad"), pcad);
              histos.fill(HIST("TrueMatch_pCharm_pt"), data_r.at(0).at(10));
              histos.fill(HIST("TrueMatch_pCharm_pcad_pt"), pcad, data_r.at(0).at(10));
              histos.fill(HIST("True_pcar"), data_r.at(0).at(0));
              if(truepair.mompdg()==421){
                histos.fill(HIST("True_InvMass"), data_r.at(0).at(14));
              }
            }else if(truepair.particletype()==2){
              histos.fill(HIST("TrueMatch_npCharm_pcad"), pcad);

            }else if(truepair.particletype()==1){
              histos.fill(HIST("TrueMatch_Beauty_pcad"), pcad);

            }
          }else if(data_r.at(0).at(13)==0){
            if(truepair.particletype()==3){
              histos.fill(HIST("FakeMatch_pCharm_pcad"), pcad);
              histos.fill(HIST("FakeMatch_pCharm_pt"), data_r.at(0).at(10));
              histos.fill(HIST("FakeMatch_pCharm_pcad_pt"), pcad, data_r.at(0).at(10));
              histos.fill(HIST("FakeDCAXY_Diff"), data_r.at(0).at(9));
              histos.fill(HIST("Fake_pcar"), data_r.at(0).at(0));
              if(truepair.mompdg()==421){
                histos.fill(HIST("Fake_InvMass"), data_r.at(0).at(14));
              }

            }else if(truepair.particletype()==2){
              histos.fill(HIST("FakeMatch_npCharm_pcad"), pcad);

            }else if(truepair.particletype()==1){
              histos.fill(HIST("FakeMatch_Beauty_pcad"), pcad);

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
    adaptAnalysisTask<McInformation>(cfgc),
    adaptAnalysisTask<MyAnalysisTask>(cfgc)
  };
}