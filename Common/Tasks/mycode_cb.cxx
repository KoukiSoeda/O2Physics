/// \author Koki Soeda
/// \since 13/09/2023

#include <iostream>
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
DECLARE_SOA_COLUMN(CosSim, cossim, double);
DECLARE_SOA_COLUMN(TotalD, totald, double);
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
                  truepair::CosSim,
                  truepair::TotalD,
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

    histos.add("ParticleType", "ParticleType", kTH1F, {axisParType});
    histos.add("DecayL", "DecayL", kTH1F, {axisZ});

    histos.add("Beauty_mu_DCAT", "Beauty_mu_DCAT", kTH1F, {axisDCAT});
    histos.add("Beauty_pair_DCAT", "Beauty_pair_DCAT", kTH1F, {axisDCAT});
    histos.add("Beauty_vs_mu_pt", "Beauty_vs_mu_pt", kTH2F, {axispT, axispT});
    histos.add("Beauty_pcaz", "Beauty_pcaz", kTH1F, {axisZ});
    histos.add("Beauty_totald", "Beauty_totald", kTH1F, {axisAdd});

    histos.add("pCharm_mu_DCAT", "pCharm_mu_DCAT", kTH1F, {axisDCAT});
    histos.add("pCharm_pair_DCAT", "pCharm_pair_DCAT", kTH1F, {axisDCAT});
    histos.add("pCharm_vs_mu_pt", "pCharm_vs_mu_pt", kTH2F, {axispT, axispT});
    histos.add("pCharm_pcaz", "pCharm_pcaz", kTH1F, {axisZ});
    histos.add("pCharm_totald", "pCharm_totald", kTH1F, {axisAdd});

    histos.add("npCharm_mu_DCAT", "npCharm_mu_DCAT", kTH1F, {axisDCAT});
    histos.add("npCharm_pair_DCAT", "npCharm_pair_DCAT", kTH1F, {axisDCAT});
    histos.add("npCharm_pcaz", "npCharm_pcaz", kTH1F, {axisZ});
    histos.add("npCharm_totald", "npCharm_totald", kTH1F, {axisAdd});

    histos.add("Light_mu_DCAT", "Light_mu_DCAT", kTH1F, {axisDCAT});
    histos.add("Light_pair_DCAT", "Light_pair_DCAT", kTH1F, {axisDCAT});
    histos.add("Light_pcaz", "Light_pcaz", kTH1F, {axisZ});
    histos.add("Light_totald", "Light_totald", kTH1F, {axisAdd});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
               soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
               soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
               aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
               aod::McCollisions const&)
  {
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
          if(fwdtrack.trackType()!=0) continue;
          auto mcParticle_fwd = fwdtrack.mcParticle();
          if(fabs(mcParticle_fwd.pdgCode())!=13) continue;

          int mucolid = fwdtrack.collisionId();
          int mcmuid = mcParticle_fwd.globalIndex();
          auto mumom = mcParticle_fwd.mothers_as<aod::McParticles>().back();
          auto momid = mumom.globalIndex();
          double mompt = mumom.pt();
          auto Daughters = mumom.daughters_as<aod::McParticles>();

          int muid=0, pid_1=0, pid_2=0, ppdg_1=0, ppdg_2=0;
          double sverx, svery, sverz;
          bool HasPair = false;
          for(auto& Daughter : Daughters){
            int d_pdg = fabs(Daughter.pdgCode());
            if(d_pdg==13 && Daughter.globalIndex()==mcmuid){
              muid = Daughter.globalIndex();
              sverx = Daughter.vx();
              svery = Daughter.vy();
              sverz = Daughter.vz();
            }else if(d_pdg!=13 && d_pdg>20 && (d_pdg%2)!=0 && pid_1==0){
              pid_1 = Daughter.globalIndex();
              ppdg_1 = d_pdg;
              HasPair = true;
            }else if(d_pdg!=13 && d_pdg>20 && (d_pdg%2)!=0 && pid_1!=0){
              pid_2 = Daughter.globalIndex();
              ppdg_2 = d_pdg;
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
          double mudcaX, mudcaY, mumftX, mumftY, mumftZ, muchi2, muphi, mutgl, musignept, mupt;
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
              if(find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack.globalIndex()) != ambTrackIds.end()) continue;
              auto mcParticle_mft = mfttrack.mcParticle();
              if(mcParticle_mft.globalIndex()==muid){
                auto mft_mom_mu = mcParticle_mft.mothers_as<aod::McParticles>().back();
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
                }
              }

              if(mcParticle_mft.globalIndex()!=muid){
                if(mcParticle_mft.has_mothers()){
                  auto mft_mom = mcParticle_mft.mothers_as<aod::McParticles>().back();

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
                  cout << mfttrack.x()-Col_x << "," << mfttrack.y()-Col_y << "," << mftpars2.getX()-Col_x << "," << mftpars2.getY()-Col_y << endl;
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
              auto pcaX = ((predict_mux+predict_p1x)/2)-Col_x;
              auto pcaY = ((predict_muy+predict_p1y)/2)-Col_y;
              auto pcaZ = ((predict_muz+predict_p1z)/2)-Col_z;
              auto pcaD = sqrt(pow(pcaX,2)+pow(pcaY,2)+pow(pcaZ,2));

              // Cosine Similarity
              auto vecx_mu = mumftX - mudcaX;
              auto vecy_mu = mumftY - mudcaY;
              auto vecz_mu = mumftZ - Col_z;
              auto vecx_p1 = p1mftX - p1dcaX;
              auto vecy_p1 = p1mftY - p1dcaY;
              auto vecz_p1 = p1mftZ - Col_z;
              auto cosxy = (vecx_mu*vecx_p1 + vecy_mu*vecy_p1 + vecz_mu*vecz_p1)/((sqrt(pow(vecx_mu,2)+pow(vecy_mu,2)+pow(vecz_mu,2)))*(sqrt(pow(vecx_p1,2)+pow(vecy_p1,2)+pow(vecz_p1,2))));

              // xy Plane Total Distance
              double dist_sum = 0;
              if(muondata.size()==pair1data.size()){
                for(int i=0; i<muondata.size(); i++){
                  auto dist = sqrt(pow(pair1data.at(i).at(0)-muondata.at(i).at(0),2)+pow(pair1data.at(i).at(1)-muondata.at(i).at(1),2));
                  dist_sum += dist;
                }
              }
              mcpairtable(Col_id, Col_x, Col_y, Col_z, muid, pid_1, p1pdg, p1dcaX, p1dcaY, pid_2, p2pdg, p2dcaX, p2dcaY, sverx, svery, sverz, pcaX, pcaY, pcaZ, pcaD, predict_r, cosxy, dist_sum, PartType);

              if(PartType==1){ // Besuty Decay
                histos.fill(HIST("Beauty_mu_DCAT"), dcat_mu);
                histos.fill(HIST("Beauty_pair_DCAT"), dcat_pair);
                histos.fill(HIST("Beauty_vs_mu_pt"), mupt, mompt);
                histos.fill(HIST("Beauty_pcaz"), pcaZ);
                histos.fill(HIST("Beauty_totald"), dist_sum);
              }else if(PartType==2){ // Non Prompt Charm
                histos.fill(HIST("npCharm_mu_DCAT"), dcat_mu);
                histos.fill(HIST("npCharm_pair_DCAT"), dcat_pair);
                histos.fill(HIST("npCharm_pcaz"), pcaZ);
                histos.fill(HIST("npCharm_totald"), dist_sum);
              }else if(PartType==3){ // Prompt Charm
                histos.fill(HIST("pCharm_mu_DCAT"), dcat_mu);
                histos.fill(HIST("pCharm_pair_DCAT"), dcat_pair);
                histos.fill(HIST("pCharm_vs_mu_pt"), mupt, mompt);
                histos.fill(HIST("pCharm_pcaz"), pcaZ);
                histos.fill(HIST("pCharm_totald"), dist_sum);
              }else if(PartType==4){ // LF Decay
                histos.fill(HIST("Light_mu_DCAT"), dcat_mu);
                histos.fill(HIST("Light_pair_DCAT"), dcat_pair);
                histos.fill(HIST("Light_pcaz"), pcaZ);
                histos.fill(HIST("Light_totald"), dist_sum);
              }
            }
          }
        }
      }
    }
  }
};

// 10/17やること -> Mask状態でのCodeを書く。
/*struct MyAnalysisTask{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&){
    histos.add()
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

  }
};*/

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<McInformation>(cfgc)
    //adaptAnalysisTask<MyAnalysisTask>(cfgc)
  };
}