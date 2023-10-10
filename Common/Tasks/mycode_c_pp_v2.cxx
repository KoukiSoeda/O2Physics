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
/*
namespace o2::aod
{
namespace truepair
{
DECLARE_SOA_COLUMN(ColId, colid, int);
DECLARE_SOA_COLUMN(MuId, muid, int);
DECLARE_SOA_COLUMN(Mupdg, mupdg, int);
DECLARE_SOA_COLUMN(KId, kid, int);
DECLARE_SOA_COLUMN(Kpdg, kpdg, int);
DECLARE_SOA_COLUMN(KdcaX, kdcax, double);
DECLARE_SOA_COLUMN(KdcaY, kdcay, double);
DECLARE_SOA_COLUMN(ColX, colx, double);
DECLARE_SOA_COLUMN(ColY, coly, double);
DECLARE_SOA_COLUMN(ColZ, colz, double);
DECLARE_SOA_COLUMN(SecVerX, secverx, double);
DECLARE_SOA_COLUMN(SecVerY, secvery, double);
DECLARE_SOA_COLUMN(SecVerZ, secverz, double);
DECLARE_SOA_COLUMN(PcaX, pcax, double);
DECLARE_SOA_COLUMN(PcaY, pcay, double);
DECLARE_SOA_COLUMN(PcaZ, pcaz, double);
DECLARE_SOA_COLUMN(CosSim, cossim, double);
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
}
DECLARE_SOA_TABLE(MCPair, "AOD", "MCPAIR",
                  truepair::ColId,
                  truepair::MuId,
                  truepair::Mupdg,
                  truepair::KId,
                  truepair::Kpdg,
                  truepair::KdcaX,
                  truepair::KdcaY,
                  truepair::ColX,
                  truepair::ColY,
                  truepair::ColZ,
                  truepair::SecVerX,
                  truepair::SecVerY,
                  truepair::SecVerZ,
                  truepair::PcaX,
                  truepair::PcaY,
                  truepair::PcaZ,
                  truepair::CosSim);
DECLARE_SOA_TABLE(MuInfo, "AOD", "MUINFO",
                  muinfo::ColId,
                  muinfo::MuId,
                  muinfo::MftX,
                  muinfo::MftY,
                  muinfo::MftZ,
                  muinfo::Chi2,
                  muinfo::Phi,
                  muinfo::Tgl,
                  muinfo::SignedPt);
}

struct McInformation{
  Produces<aod::MCPair> mcpairtable;
  Produces<aod::MuInfo> muinfotable;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsZ{"nBinsZ", 110010, "N bins in Z axis"};

  void init(InitContext const&){
    const AxisSpec axisZ{nBinsZ, -60.0005, 50.0005, "z(cm)"};
    const AxisSpec axisXY{200010, -10.0005, 10.0005, "cm"};
    const AxisSpec axisDCAT{100010, -0.0005, 10.0005, "DCAT"};
    const AxisSpec axisCos{1000, 0, 1, ""};
    const AxisSpec axisAdd{5000, 0, 500, "sum(cm)"};
    const AxisSpec axisCluster{101, -0.5, 100.5, ""};

    histos.add("PredictZ", "PredictZ", kTH1F, {axisZ});
    histos.add("PredictX", "PredictX", kTH1F, {axisXY});
    histos.add("PredictY", "PredictY", kTH1F, {axisXY});
    histos.add("PredictR", "PredictR", kTH1F, {axisXY});
    histos.add("DCAT_Muon", "DCAT_Muon", kTH1F, {axisDCAT});
    histos.add("DCAT_Kaon", "DCAT_Kaon", kTH1F, {axisDCAT});
    histos.add("nClusters", "nClusters", kTH1F, {axisCluster});
    histos.add("CosSimilarity_muk", "CosSimilarity_muk", kTH1F, {axisCos});
    histos.add("Add_distmuk", "Add_distmuk", kTH1F, {axisAdd});
    histos.add("AddDistmuk_Cossim", "AddDistmuk_Cossim", kTH2F, {axisAdd, axisCos});

  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
               soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
               soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
               aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
               aod::McCollisions const&)
  {
    if(collision.has_mcCollision()){
      if(collision.posZ()<10 && collision.posZ()>-10){
        vector<uint64_t> ambTrackIds;
        ambTrackIds.clear();
        for(auto& amfttrack : amfttracks){
          ambTrackIds.push_back(amfttrack.mfttrackId());
        }
        for(auto& fwdtrack : fwdtracks){
          if(!fwdtrack.has_collision() || !fwdtrack.has_mcParticle()) continue;
          if(fwdtrack.trackType()!=0) continue;
          auto mcParticle_fwd = fwdtrack.mcParticle();
          if(fabs(mcParticle_fwd.pdgCode())==13 && mcParticle_fwd.has_mothers()){
            auto fwdColId = fwdtrack.collisionId();
            auto mcMoms = mcParticle_fwd.mothers_as<aod::McParticles>();
            auto mcMom_b = mcMoms.back();
            auto mcMom_f = mcMoms.front();

            if(collision.globalIndex()==fwdColId){
              if(mcMom_f.pdgCode()==mcMom_b.pdgCode() && fabs(mcMom_b.pdgCode())==421){
                auto Daughters = mcMom_b.daughters_as<aod::McParticles>();
                int dc = 0;
                int fc = 0;
                int mcmucolid, mckcolid, mucolid, kcolid, mcmuid, mupdg, mckid, kpdg, kcluster;
                double mudcaX, mudcaY, mumftX, mumftY, mumftZ, kdcaX, kdcaY, kmftX, kmftY, kmftZ, mumcX, mumcY, mumcZ, kmcX, kmcY, kmcZ;
                double muchi2, muphi, mutgl, musignept;
                const double Col_x = collision.posX();
                const double Col_y = collision.posY();
                const double Col_z = collision.posZ();
                vector<vector<double>> muondata{};
                vector<vector<double>> kaondata{};

                for(auto& mfttrack : mfttracks){
                  if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
                  if(find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack.globalIndex()) != ambTrackIds.end()) continue;
                  auto mcParticle_mft = mfttrack.mcParticle();
                  for(auto& Daughter : Daughters){
                    if(fabs(Daughter.pdgCode())==13){
                      if(mfttrack.mcParticleId()==Daughter.globalIndex() && mfttrack.collisionId()==fwdColId){
                        double mftchi2 = mfttrack.chi2();
                        SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                        vector<double> mftv1;
                        SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                        o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
                        mftpars1.propagateToZlinear(Col_z);
                        mudcaX = mftpars1.getX();
                        mudcaY = mftpars1.getY();
                        mumftX = mfttrack.x();
                        mumftY = mfttrack.y();
                        mumftZ = mfttrack.z();
                        mumcX = mcParticle_mft.vx();
                        mumcY = mcParticle_mft.vy();
                        mumcZ = mcParticle_mft.vz();
                        mucolid = mfttrack.collisionId();
                        mcmucolid = mcParticle_mft.mcCollisionId();
                        mcmuid = mfttrack.mcParticleId();
                        mupdg = mcParticle_mft.pdgCode();
                        muchi2 = mftchi2;
                        muphi = mfttrack.phi();
                        mutgl = mfttrack.tgl();
                        musignept = mfttrack.signed1Pt();
                        dc++;
                        for(int i=0; i<10000; i++){
                          auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                          if((cut_z-Col_z)>0) break;
                          mftpars1.propagateToZlinear(cut_z);
                          vector<double> temp = {mftpars1.getX(), mftpars1.getY()};
                          muondata.push_back(temp);
                        }
                      }
                    }
                    if(fabs(Daughter.pdgCode())==321){
                      if(mfttrack.mcParticleId()==Daughter.globalIndex() && mfttrack.collisionId()==fwdColId){
                        double mftchi2 = mfttrack.chi2();
                        SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                        vector<double> mftv1;
                        SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                        o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars, mftcovs, mftchi2};
                        mftpars2.propagateToZlinear(Col_z);
                        kdcaX = mftpars2.getX();
                        kdcaY = mftpars2.getY();
                        kmftX = mfttrack.x();
                        kmftY = mfttrack.y();
                        kmftZ = mfttrack.z();
                        kmcX = mcParticle_mft.vx();
                        kmcY = mcParticle_mft.vy();
                        kmcZ = mcParticle_mft.vz();
                        kcolid = mfttrack.collisionId();
                        mckcolid = mcParticle_mft.mcCollisionId();
                        mckid = mfttrack.mcParticleId();
                        kpdg = mcParticle_mft.pdgCode();
                        kcluster = mfttrack.nClusters();
                        fc++;
                        for(int i=0; i<=10000; i++){
                          auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                          if((cut_z-Col_z)>0) break;
                          mftpars2.propagateToZlinear(cut_z);
                          vector<double> temp2 = {mftpars2.getX(), mftpars2.getY()};
                          kaondata.push_back(temp2);
                        }
                      }
                    }
                  }
                }
                if(mucolid==kcolid && dc+fc==2){
                  if(mumcZ==kmcZ && collision.mcCollisionId()==mcmucolid){
                    histos.fill(HIST("DCAT_Muon"), sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2)));
                    histos.fill(HIST("DCAT_Kaon"), sqrt(pow(kdcaX-Col_x,2)+pow(kdcaY-Col_y,2)));
                    histos.fill(HIST("nClusters"), kcluster);

                    muinfotable(mucolid,mcmuid,mumftX,mumftY,mumftZ,muchi2,muphi,mutgl,musignept);

                    //ベクトル法
                    auto unit_Na = sqrt(pow(mumftX-mudcaX,2)+pow(mumftY-mudcaY,2)+pow(mumftZ-Col_z,2));
                    auto unit_Nc = sqrt(pow(kmftX-kdcaX,2)+pow(kmftY-kdcaY,2)+pow(kmftZ-Col_z,2));
                    auto Nax = (mumftX-mudcaX)/unit_Na;
                    auto Nay = (mumftY-mudcaY)/unit_Na;
                    auto Naz = (mumftZ-Col_z)/unit_Na;
                    auto Ncx = (kmftX-kdcaX)/unit_Nc;
                    auto Ncy = (kmftY-kdcaY)/unit_Nc;
                    auto Ncz = (kmftZ-Col_z)/unit_Nc;
                    auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
                    auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
                    auto A3 = (mudcaX-kdcaX)*Nax + (mudcaY-kdcaY)*Nay + (Col_z-Col_z)*Naz;
                    auto B1 = A2;
                    auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
                    auto B3 = (kdcaX-mudcaX)*Ncx + (kdcaY-mudcaY)*Ncy + (Col_z-Col_z)*Ncz;
                    auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
                    auto s = -((A2*t+A3)/A1);

                    double predict_mux = mudcaX + s*Nax;
                    double predict_muy = mudcaY + s*Nay;
                    double predict_muz = Col_z + s*Naz;
                    double predict_kx = kdcaX + t*Ncx;
                    double predict_ky = kdcaY + t*Ncy;
                    double predict_kz = Col_z + t*Ncz;
                    double predict_r = sqrt(pow(predict_mux-predict_kx,2)+pow(predict_muy-predict_ky,2)+pow(predict_muz-predict_kz,2));
                    auto pcaX = ((predict_mux+predict_kx)/2)-Col_x;
                    auto pcaY = ((predict_muy+predict_ky)/2)-Col_y;
                    auto pcaZ = ((predict_muz+predict_kz)/2)-Col_z;

                    histos.fill(HIST("PredictZ"), fabs(pcaZ));
                    histos.fill(HIST("PredictX"), pcaX);
                    histos.fill(HIST("PredictY"), pcaY);
                    histos.fill(HIST("PredictR"), predict_r);

                    //cos類似度分布
                    auto vecx_mu = mumftX - mudcaX;
                    auto vecy_mu = mumftY - mudcaY;
                    auto vecz_mu = mumftZ - Col_z;
                    auto vecx_k = kmftX - kdcaX;
                    auto vecy_k = kmftY - kdcaY;
                    auto vecz_k = kmftZ - Col_z;
                    auto cosxy = (vecx_mu*vecx_k + vecy_mu*vecy_k + vecz_mu*vecz_k)/((sqrt(pow(vecx_mu,2)+pow(vecy_mu,2)+pow(vecz_mu,2)))*(sqrt(pow(vecx_k,2)+pow(vecy_k,2)+pow(vecz_k,2))));
                    histos.fill(HIST("CosSimilarity_muk"), cosxy);
                    mcpairtable(mucolid,mcmuid,mupdg,mckid,kpdg,kdcaX-Col_x,kdcaY-Col_y,Col_x,Col_y,Col_z,mumcX,mumcY,mumcZ,pcaX-Col_x,pcaY-Col_y,pcaZ-Col_z,cosxy);
                  
                    //x-y平面距離合計
                    double dist_sum = 0;
                    if(muondata.size()==kaondata.size()){
                      for(int i=0; i<muondata.size(); i++){
                        auto dist = sqrt(pow(kaondata.at(i).at(0)-muondata.at(i).at(0),2)+pow(kaondata.at(i).at(1)-muondata.at(i).at(1),2));
                        dist_sum += dist;
                      }
                    }
                    histos.fill(HIST("Add_distmuk"), dist_sum);
                    histos.fill(HIST("AddDistmuk_Cossim"), dist_sum, cosxy);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
};

struct McFakeInfo{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsZ{"nBinsZ", 110010, "N bins in Z axis"};

  void init(InitContext const&){
    const AxisSpec axisZ{nBinsZ, -60.0005, 50.0005, "z(cm)"};
    const AxisSpec axisXY{200010, -10.0005, 10.0005, "cm"};
    const AxisSpec axisDCAT{100010, -0.0005, 10.0005, "DCAT"};
    const AxisSpec axisCos{1000, 0, 1, ""};
    const AxisSpec axisAdd{5000, 0, 500, "sum(cm)"};
    const AxisSpec axisCluster{101, -0.5, 100.5, ""};

    histos.add("FakePCAZ", "FakePCAZ", kTH1F, {axisZ});
    histos.add("FakePCAX", "FakePCAX", kTH1F, {axisXY});
    histos.add("FakePCAY", "FakePCAY", kTH1F, {axisXY});
    histos.add("FakeR", "FakeR", kTH1F, {axisXY});
    histos.add("DCAT_all", "DCAT_all", kTH1F, {axisDCAT});
    histos.add("DCAT_Sec", "DCAT_Sec", kTH1F, {axisDCAT});
    histos.add("nClusters", "nClusters", kTH1F, {axisCluster});
    histos.add("CosSimilarity_fake", "CosSimilarity_fake", kTH1F, {axisCos});
    histos.add("Add_distfake", "Add_distfake", kTH1F, {axisAdd});
    histos.add("AddDistfake_Cossim", "AddDistfake_Cossim", kTH2F, {axisAdd, axisCos});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
							 soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
							 soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
							 aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
							 aod::McCollisions const&,
               aod::MCPair const& truepairs,
               aod::MuInfo const& muinfotable)
  {
    if(collision.has_mcCollision()){
      if(collision.posZ()<10 && collision.posZ()>-10){
        vector<uint64_t> ambTrackIds;
        ambTrackIds.clear();
        for(auto& amfttrack : amfttracks){
          ambTrackIds.push_back(amfttrack.mfttrackId());
        }
        for(auto& truepair : truepairs){
          for(auto& muinfo : muinfotable){
            if(muinfo.muid()==truepair.muid() && muinfo.colid()==truepair.colid()){
              const double Col_x = collision.posX();
              const double Col_y = collision.posY();
              const double Col_z = collision.posZ();
              
              vector<vector<double>> muondata;
              SMatrix5 mftpars(muinfo.x(), muinfo.y(), muinfo.phi(), muinfo.tgl(), muinfo.signed1Pt());
              vector<double> mftv1;
              SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
              o2::track::TrackParCovFwd mftpars1{muinfo.z(), mftpars, mftcovs, muinfo.chi2()};
              
              for(int i=0; i<10000; i++){
                auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                if((cut_z-Col_z)>0) break;
                mftpars1.propagateToZlinear(cut_z);
                vector<double> temp = {mftpars1.getX(), mftpars1.getY()};
                muondata.push_back(temp);
              }
              mftpars1.propagateToZlinear(Col_z);
              double mudcaX = mftpars1.getX();
              double mudcaY = mftpars1.getY();

              for(auto& mfttrack : mfttracks){
                if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
                if(truepair.colid()!=mfttrack.collisionId()) continue;
                if(find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack.globalIndex()) != ambTrackIds.end()) continue;
                double kdcaX, kdcaY, kmftX, kmftY, kmftZ; 
                
                if(mfttrack.mcParticleId()!=truepair.kid() && mfttrack.mcParticleId()!=truepair.muid()){
                  auto mcParticle_mft = mfttrack.mcParticle();
                  double mftchi2 = mfttrack.chi2();
                  SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                  vector<double> mftv2;
                  SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
                  o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars, mftcovs2, mftchi2};
                  mftpars2.propagateToZlinear(Col_z);
                  kdcaX = mftpars2.getX();
                  kdcaY = mftpars2.getY();
                  kmftX = mfttrack.x();
                  kmftY = mfttrack.y();
                  kmftZ = mfttrack.z();

                  histos.fill(HIST("DCAT_all"), sqrt(pow(kdcaX-Col_x,2)+pow(kdcaY-Col_y,2)));
                  if(mcParticle_mft.has_mothers()){
                    auto mcMoms = mcParticle_mft.mothers_as<aod::McParticles>();
                    auto mcMom_b = mcMoms.back();
                    auto mcMom_f = mcMoms.front();
                    if(fabs(mcParticle_mft.pdgCode())==211){
                      if(fabs(mcMom_b.pdgCode())<=37) break;
                      //LOGF(info, "Particle has mothers: %d, %d    MFT PDG: %d", mcMom_f.pdgCode(), mcMom_b.pdgCode(), mcParticle_mft.pdgCode());
                      histos.fill(HIST("DCAT_Sec"), sqrt(pow(kdcaX-Col_x,2)+pow(kdcaY-Col_y,2))); 
                    }
                  }
                  
                  histos.fill(HIST("nClusters"), mfttrack.nClusters());

                  vector<vector<double>> kaondata;
                  for(int i=0; i<=10000; i++){
                    auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                    if((cut_z-Col_z)>0) break;
                    mftpars2.propagateToZlinear(cut_z);
                    vector<double> temp2 = {mftpars2.getX(), mftpars2.getY()};
                    kaondata.push_back(temp2);
                  }

                  auto unit_Na = sqrt(pow(muinfo.x()-mudcaX,2)+pow(muinfo.y()-mudcaY,2)+pow(muinfo.z()-Col_z,2));
                  auto unit_Nc = sqrt(pow(kmftX-kdcaX,2)+pow(kmftY-kdcaY,2)+pow(kmftZ-Col_z,2));
                  auto Nax = (muinfo.x()-mudcaX)/unit_Na;
                  auto Nay = (muinfo.y()-mudcaY)/unit_Na;
                  auto Naz = (muinfo.z()-Col_z)/unit_Na;
                  auto Ncx = (kmftX-kdcaX)/unit_Nc;
                  auto Ncy = (kmftY-kdcaY)/unit_Nc;
                  auto Ncz = (kmftZ-Col_z)/unit_Nc;
                  auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
                  auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
                  auto A3 = (mudcaX-kdcaX)*Nax + (mudcaY-kdcaY)*Nay + (Col_z-Col_z)*Naz;
                  auto B1 = A2;
                  auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
                  auto B3 = (kdcaX-mudcaX)*Ncx + (kdcaY-mudcaY)*Ncy + (Col_z-Col_z)*Ncz;
                  auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
                  auto s = -((A2*t+A3)/A1);

                  double predict_mux = mudcaX + s*Nax;
                  double predict_muy = mudcaY + s*Nay;
                  double predict_muz = Col_z + s*Naz;
                  double predict_kx = kdcaX + t*Ncx;
                  double predict_ky = kdcaY + t*Ncy;
                  double predict_kz = Col_z + t*Ncz;
                  double predict_r = sqrt(pow(predict_mux-predict_kx,2)+pow(predict_muy-predict_ky,2)+pow(predict_muz-predict_kz,2));
                  auto pcaX = ((predict_mux+predict_kx)/2)-Col_x;
                  auto pcaY = ((predict_muy+predict_ky)/2)-Col_y;
                  auto pcaZ = ((predict_muz+predict_kz)/2)-Col_z;
                  
                  histos.fill(HIST("FakePCAZ"), fabs(pcaZ));
                  histos.fill(HIST("FakePCAX"), pcaX);
                  histos.fill(HIST("FakePCAY"), pcaY);
                  histos.fill(HIST("FakeR"), predict_r);

                  auto vecx_mu = muinfo.x() - mudcaX;
                  auto vecy_mu = muinfo.y() - mudcaY;
                  auto vecz_mu = muinfo.z() - Col_z;
                  auto vecx_k = kmftX - kdcaX;
                  auto vecy_k = kmftY - kdcaY;
                  auto vecz_k = kmftZ - Col_z;
                  auto cosxy = (vecx_mu*vecx_k + vecy_mu*vecy_k + vecz_mu*vecz_k)/((sqrt(pow(vecx_mu,2)+pow(vecy_mu,2)+pow(vecz_mu,2)))*(sqrt(pow(vecx_k,2)+pow(vecy_k,2)+pow(vecz_k,2))));
                  histos.fill(HIST("CosSimilarity_fake"), cosxy);

                  double dist_sum = 0;
                  if(muondata.size()==kaondata.size()){
                    for(int i=0; i<muondata.size(); i++){
                      auto dist = sqrt(pow(kaondata.at(i).at(0)-muondata.at(i).at(0),2)+pow(kaondata.at(i).at(1)-muondata.at(i).at(1),2));
                      dist_sum += dist;
                    }
                  }
                  histos.fill(HIST("Add_distfake"), dist_sum);
                  histos.fill(HIST("AddDistfake_Cossim"), dist_sum, cosxy);
                }
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
  Configurable<int> nBinsZ{"nBinsZ", 110010, "N bins in Z axis"};

  void init(InitContext const&){
    const AxisSpec axisZ{nBinsZ, -60.0005, 50.0005, "z(cm)"};
    const AxisSpec axisXY{200010, -10.0005, 10.0005, "cm"};
    const AxisSpec axisAdd{5000, 0, 500, "sum(cm)"};
    histos.add("True_R_closest", "True_R_closest", kTH1F, {axisZ});
    histos.add("Fake_R_closest", "Fake_R_closest", kTH1F, {axisZ});
    histos.add("True_R_closest_dcat_0_01", "True_R_closest_dcat_0_01", kTH1F, {axisZ});
    histos.add("True_R_closest_dcat_0_002", "True_R_closest_dcat_0_002", kTH1F, {axisZ});
    histos.add("True_R_closest_dcat_002_004", "True_R_closest_dcat_002_004", kTH1F, {axisZ});
    histos.add("True_R_closest_dcat_004_006", "True_R_closest_dcat_004_006", kTH1F, {axisZ});
    histos.add("True_R_closest_dcat_006_008", "True_R_closest_dcat_006_008", kTH1F, {axisZ});
    histos.add("True_R_closest_dcat_008_01", "True_R_closest_dcat_008_01", kTH1F, {axisZ});
    histos.add("True_R_closest_dcat_01_02", "True_R_closest_dcat_01_02", kTH1F, {axisZ});
    histos.add("True_R_closest_dcat_02_03", "True_R_closest_dcat_02_03", kTH1F, {axisZ});
    histos.add("True_R_closest_dcat_03_04", "True_R_closest_dcat_03_04", kTH1F, {axisZ});
    histos.add("True_R_closest_dcat_04", "True_R_closest_dcat_04", kTH1F, {axisZ});

    histos.add("True_R_closest_pcaz_0_29", "True_R_closest_pcaz_0_29", kTH1F, {axisZ});
    histos.add("True_R_closest_pcaz_0_05", "True_R_closest_pcaz_0_05", kTH1F, {axisZ});
    histos.add("True_R_closest_pcaz_05_1", "True_R_closest_pcaz_05_1", kTH1F, {axisZ});
    histos.add("True_R_closest_pcaz_1_15", "True_R_closest_pcaz_1_15", kTH1F, {axisZ});
    histos.add("True_R_closest_pcaz_15_2", "True_R_closest_pcaz_15_2", kTH1F, {axisZ});
    histos.add("True_R_closest_pcaz_2_25", "True_R_closest_pcaz_2_25", kTH1F, {axisZ});
    histos.add("True_R_closest_pcaz_25_3", "True_R_closest_pcaz_25_3", kTH1F, {axisZ});

    histos.add("Fake_R_closest_dcat_0_01", "Fake_R_closest_dcat_0_01", kTH1F, {axisZ});
    histos.add("Fake_R_closest_dcat_0_002", "Fake_R_closest_dcat_0_002", kTH1F, {axisZ});
    histos.add("Fake_R_closest_dcat_002_004", "Fake_R_closest_dcat_002_004", kTH1F, {axisZ});
    histos.add("Fake_R_closest_dcat_004_006", "Fake_R_closest_dcat_004_006", kTH1F, {axisZ});
    histos.add("Fake_R_closest_dcat_006_008", "Fake_R_closest_dcat_006_008", kTH1F, {axisZ});
    histos.add("Fake_R_closest_dcat_008_01", "Fake_R_closest_dcat_008_01", kTH1F, {axisZ});
    histos.add("Fake_R_closest_dcat_01_02", "Fake_R_closest_dcat_01_02", kTH1F, {axisZ});
    histos.add("Fake_R_closest_dcat_02_03", "Fake_R_closest_dcat_02_03", kTH1F, {axisZ});
    histos.add("Fake_R_closest_dcat_03_04", "Fake_R_closest_dcat_03_04", kTH1F, {axisZ});
    histos.add("Fake_R_closest_dcat_04", "Fake_R_closest_dcat_04", kTH1F, {axisZ});

    histos.add("Fake_R_closest_pcaz_0_29", "Fake_R_closest_pcaz_0_29", kTH1F, {axisZ});
    histos.add("Fake_R_closest_pcaz_0_05", "Fake_R_closest_pcaz_0_05", kTH1F, {axisZ});
    histos.add("Fake_R_closest_pcaz_05_1", "Fake_R_closest_pcaz_05_1", kTH1F, {axisZ});
    histos.add("Fake_R_closest_pcaz_1_15", "Fake_R_closest_pcaz_1_15", kTH1F, {axisZ});
    histos.add("Fake_R_closest_pcaz_15_2", "Fake_R_closest_pcaz_15_2", kTH1F, {axisZ});
    histos.add("Fake_R_closest_pcaz_2_25", "Fake_R_closest_pcaz_2_25", kTH1F, {axisZ});
    histos.add("Fake_R_closest_pcaz_25_3", "Fake_R_closest_pcaz_25_3", kTH1F, {axisZ});

    histos.add("True_Distance_smallest", "True_Distance_smallest", kTH1F, {axisZ});
    histos.add("Fake_Distance_smallest", "Fake_Distance_smallest", kTH1F, {axisZ});
    histos.add("True_DCAT_0001_05", "True_DCAT_0001_05", kTH1F, {axisZ});
    histos.add("Fake_DCAT_0001_05", "Fake_DCAT_0001_05", kTH1F, {axisZ});
    histos.add("True_Cos_largest", "True_Cos_largest", kTH1F, {axisZ});
    histos.add("Fake_Cos_largest", "Fake_Cos_largest", kTH1F, {axisZ});
    histos.add("True_PCAZ_m29_p24", "True_PCAZ_m29_p24", kTH1F, {axisZ});
    histos.add("Fake_PCAZ_m29_p24", "Fake_PCAZ_m29_p24", kTH1F, {axisZ});

  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
							 soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
							 soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
							 aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
							 aod::McCollisions const&,
               aod::MCPair const& truepairs)
  {
    if(collision.has_mcCollision()){
      if(collision.posZ()<10 && collision.posZ()>-10){
        vector<uint64_t> ambTrackIds;
        ambTrackIds.clear();
        for(auto& amfttrack : amfttracks){
          ambTrackIds.push_back(amfttrack.mfttrackId());
        }
        for(auto& fwdtrack : fwdtracks){
          if(!fwdtrack.has_collision() || !fwdtrack.has_mcParticle()) continue;
          if(fwdtrack.trackType()!=0) continue;
          auto mcParticle_fwd = fwdtrack.mcParticle();
          if(fabs(mcParticle_fwd.pdgCode())!=13 || !mcParticle_fwd.has_mothers()) continue;
          auto fwdColId = fwdtrack.collisionId();
          auto fwdmuId = fwdtrack.mcParticleId();

          for(auto& truepair : truepairs){
            if(truepair.muid()==fwdmuId && truepair.colid()==fwdColId){
              double Col_x = collision.posX();
              double Col_y = collision.posY();
              double Col_z = collision.posZ();

              double mudcaX, mudcaY, mumftX, mumftY, mumftZ;
              for(auto& muontrack : mfttracks){
                if(!muontrack.has_collision() || !muontrack.has_mcParticle()) continue;
                if(muontrack.collisionId()==fwdColId && muontrack.mcParticleId()==truepair.muid()){
                  double muonchi2 = muontrack.chi2();
                  SMatrix5 mftpars(muontrack.x(), muontrack.y(), muontrack.phi(), muontrack.tgl(), muontrack.signed1Pt());
                  vector<double> mftv1;
                  SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                  o2::track::TrackParCovFwd mftpars1{muontrack.z(), mftpars, mftcovs, muonchi2};
                  mftpars1.propagateToZlinear(Col_z);
                  mudcaX = mftpars1.getX();
                  mudcaY = mftpars1.getY();
                  mumftX = muontrack.x();
                  mumftY = muontrack.y();
                  mumftZ = muontrack.z();
                  double candcaX, candcaY, canmftX, canmftY, canmftZ;

                  vector<vector<double>> data_r;
                  vector<vector<double>> data_add;
                  vector<vector<double>> data_cos;

                  for(auto& mfttrack : mfttracks){
                    if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
                    if(mfttrack.collisionId()==fwdColId && mfttrack.mcParticleId()!=truepair.muid()){
                      if(find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack.globalIndex()) != ambTrackIds.end()) continue;
                      auto mcParticle_mft = mfttrack.mcParticle();
                      double mftchi2 = mfttrack.chi2();
                      SMatrix5 mftpars_2(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                      vector<double> mftv2;
                      SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
                      o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars_2, mftcovs2, mftchi2};
                      mftpars2.propagateToZlinear(Col_z);
                      candcaX = mftpars2.getX();
                      candcaY = mftpars2.getY();
                      canmftX = mfttrack.x();
                      canmftY = mfttrack.y();
                      canmftZ = mfttrack.z();

                      //3D
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

                      double dcat = sqrt(pow(candcaX-Col_x,2)+pow(candcaY-Col_y,2));

                      double mcid = mfttrack.mcParticleId();
                      double out_colid = mfttrack.collisionId();

                      double add_dist = 0;
                      for(int i=0; i<5000; i++){
                        auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                        if((cut_z-Col_z)>0) break;
                        mftpars1.propagateToZlinear(cut_z);
                        mftpars2.propagateToZlinear(cut_z);
                        auto diff = sqrt(pow(mftpars2.getX()-mftpars1.getX(),2)+pow(mftpars2.getY()-mftpars2.getY(),2));
                        add_dist += diff;
                      }
                      double ans = 0;
                      if(mcid==truepair.kid()) ans=1;

                      vector<double> pardata_r = {r_xyz, cosxy, mcid, fabs(mcParticle_mft.pdgCode()), pcaX, pcaY, pcaZ, add_dist, dcat, out_colid, ans};
                      data_r.push_back(pardata_r);
                      vector<double> pardata_add = {add_dist, r_xyz, cosxy, mcid, fabs(mcParticle_mft.pdgCode()), pcaX, pcaY, pcaZ, dcat, ans};
                      data_add.push_back(pardata_add);
                      vector<double> pardata_cos = {cosxy, r_xyz, mcid, fabs(mcParticle_mft.pdgCode()), pcaX, pcaY, pcaZ, add_dist, dcat, ans};
                      data_cos.push_back(pardata_cos);
                      
                    }
                  }
                  sort(data_r.begin(), data_r.end());
                  sort(data_add.begin(), data_add.end());
                  sort(data_cos.rbegin(), data_cos.rend());
                  
                  //debug
                  int counter = 0;
                  for(int i=0; i<data_r.size(); i++){
                    for(int j=0; j<data_r.at(0).size(); j++){
                      if(counter==0){
                        cout << "data_r = { " << data_r.at(i).at(j) << ", ";
                        counter++;
                      }else if(j==(data_r.at(0).size()-1)){
                        if(i==(data_r.size()-1)){
                          cout << data_r.at(i).at(j) << " }" << endl;
                        }else{
                          cout << data_r.at(i).at(j) << " }, " << endl << "{ ";
                        }
                      }
                      else{
                        cout << data_r.at(i).at(j) << ", ";
                      }
                    }
                  }

                  //PCAR
                  if(data_r.at(0).at(2)==truepair.kid()){
                    histos.fill(HIST("True_R_closest"), data_r.at(0).at(6));
                    //DCAT Cut
                    if(data_r.at(0).at(8)<0.1) histos.fill(HIST("True_R_closest_dcat_0_01"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0 && data_r.at(0).at(8)<0.02) histos.fill(HIST("True_R_closest_dcat_0_002"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.02 && data_r.at(0).at(8)<0.04) histos.fill(HIST("True_R_closest_dcat_002_004"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.04 && data_r.at(0).at(8)<0.06) histos.fill(HIST("True_R_closest_dcat_004_006"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.06 && data_r.at(0).at(8)<0.08) histos.fill(HIST("True_R_closest_dcat_006_008"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.08 && data_r.at(0).at(8)<0.1) histos.fill(HIST("True_R_closest_dcat_008_01"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.1 && data_r.at(0).at(8)<0.2) histos.fill(HIST("True_R_closest_dcat_01_02"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.2 && data_r.at(0).at(8)<0.3) histos.fill(HIST("True_R_closest_dcat_02_03"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.3 && data_r.at(0).at(8)<0.4) histos.fill(HIST("True_R_closest_dcat_03_04"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.4) histos.fill(HIST("True_R_closest_dcat_04"), data_r.at(0).at(6));

                    //PCAZ Cut
                    if(fabs(data_r.at(0).at(6))>0 && fabs(data_r.at(0).at(6))<2.9) histos.fill(HIST("True_R_closest_pcaz_0_29"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>0 && fabs(data_r.at(0).at(6))<0.5) histos.fill(HIST("True_R_closest_pcaz_0_05"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>=0.5 && fabs(data_r.at(0).at(6))<1) histos.fill(HIST("True_R_closest_pcaz_05_1"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>=1 && fabs(data_r.at(0).at(6))<1.5) histos.fill(HIST("True_R_closest_pcaz_1_15"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>=1.5 && fabs(data_r.at(0).at(6))<2) histos.fill(HIST("True_R_closest_pcaz_15_2"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>=2 && fabs(data_r.at(0).at(6))<2.5) histos.fill(HIST("True_R_closest_pcaz_2_25"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>=2.5 && fabs(data_r.at(0).at(6))<3) histos.fill(HIST("True_R_closest_pcaz_25_3"), fabs(data_r.at(0).at(6)));
                  }else{
                    histos.fill(HIST("Fake_R_closest"), data_r.at(0).at(6));
                    //DCAT Cut
                    if(data_r.at(0).at(8)<0.1) histos.fill(HIST("Fake_R_closest_dcat_0_01"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0 && data_r.at(0).at(8)<0.02) histos.fill(HIST("Fake_R_closest_dcat_0_002"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.02 && data_r.at(0).at(8)<0.04) histos.fill(HIST("Fake_R_closest_dcat_002_004"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.04 && data_r.at(0).at(8)<0.06) histos.fill(HIST("Fake_R_closest_dcat_004_006"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.06 && data_r.at(0).at(8)<0.08) histos.fill(HIST("Fake_R_closest_dcat_006_008"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.08 && data_r.at(0).at(8)<0.1) histos.fill(HIST("Fake_R_closest_dcat_008_01"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.1 && data_r.at(0).at(8)<0.2) histos.fill(HIST("Fake_R_closest_dcat_01_02"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.2 && data_r.at(0).at(8)<0.3) histos.fill(HIST("Fake_R_closest_dcat_02_03"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.3 && data_r.at(0).at(8)<0.4) histos.fill(HIST("Fake_R_closest_dcat_03_04"), data_r.at(0).at(6));
                    if(data_r.at(0).at(8)>=0.4) histos.fill(HIST("Fake_R_closest_dcat_04"), data_r.at(0).at(6));

                    //PCAZ Cut
                    if(fabs(data_r.at(0).at(6))>0 && fabs(data_r.at(0).at(6))<2.9) histos.fill(HIST("Fake_R_closest_pcaz_0_29"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>0 && fabs(data_r.at(0).at(6))<0.5) histos.fill(HIST("Fake_R_closest_pcaz_0_05"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>=0.5 && fabs(data_r.at(0).at(6))<1) histos.fill(HIST("Fake_R_closest_pcaz_05_1"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>=1 && fabs(data_r.at(0).at(6))<1.5) histos.fill(HIST("Fake_R_closest_pcaz_1_15"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>=1.5 && fabs(data_r.at(0).at(6))<2) histos.fill(HIST("Fake_R_closest_pcaz_15_2"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>=2 && fabs(data_r.at(0).at(6))<2.5) histos.fill(HIST("Fake_R_closest_pcaz_2_25"), fabs(data_r.at(0).at(6)));
                    if(fabs(data_r.at(0).at(6))>=2.5 && fabs(data_r.at(0).at(6))<3) histos.fill(HIST("Fake_R_closest_pcaz_25_3"), fabs(data_r.at(0).at(6)));
                  }
                  
                  //Total distance
                  if(data_add.at(0).at(3)==truepair.kid()){
                    histos.fill(HIST("True_Distance_smallest"), data_add.at(0).at(7));
                  }else{
                    histos.fill(HIST("Fake_Distance_smallest"), data_add.at(0).at(7));
                  }

                  //DCAT
                  for(int i=0; i<data_r.size(); i++){
                    if(data_r.at(i).at(8)>0.001 && data_r.at(i).at(8)<0.5){
                      if(data_r.at(i).at(2)==truepair.kid()){
                        histos.fill(HIST("True_DCAT_0001_05"), data_r.at(i).at(6));
                      }else{
                        histos.fill(HIST("Fake_DCAT_0001_05"), data_r.at(i).at(6));
                      }
                    }
                  }

                  //cos
                  if(data_cos.at(0).at(2)==truepair.kid()){
                    histos.fill(HIST("True_Cos_largest"), data_cos.at(0).at(6));
                  }else{
                    histos.fill(HIST("Fake_Cos_largest"), data_cos.at(0).at(6));
                  }

                  //PCAZ
                  for(int i=0; i<data_r.size(); i++){
                    if(data_r.at(i).at(6)<2.4 && data_r.at(i).at(6)>(-2.9)){
                      if(data_r.at(i).at(2)==truepair.kid()){
                        histos.fill(HIST("True_PCAZ_m29_p24"), data_r.at(i).at(6));
                      }else{
                        histos.fill(HIST("Fake_PCAZ_m29_p24"), data_r.at(i).at(6));
                      }
                    }
                  }

                }
              }
            }
          }
        }
      }
    }
  }
};*/

struct Mytest{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&){
    const AxisSpec axistrackType{6, -0.5, 5.5, "TrackType"};
    histos.add("TrackType", "TrackType", kTH1F, {axistrackType});
  }
  void process(aod::Collisions const& collisions,
               soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
							 soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
							 aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks)
  {
    for(auto& collision : collisions){
      if(collision.posZ()<10 && collision.posZ()>-10){
        for(auto& fwdtrack : fwdtracks){
          if(fwdtrack.has_collision()){
            if(fwdtrack.collisionId()!=collision.globalIndex() || !fwdtrack.has_mcParticle()) continue;
            histos.fill(HIST("TrackType"), fwdtrack.trackType());
            if(fwdtrack.trackType()!=2) continue;
            auto mcParticle_fwd = fwdtrack.mcParticle();
            if(fabs(mcParticle_fwd.pdgCode())==13){
              //if(!mcParticle_fwd.has_mothers()){
                LOGF(info, "Collision Vertex: %lf, %lf, %lf  Particle Vertex: %lf, %lf, %lf", collision.posX(), collision.posY(), collision.posZ(), mcParticle_fwd.vx(), mcParticle_fwd.vy(), mcParticle_fwd.vz());
              //}
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
    adaptAnalysisTask<Mytest>(cfgc)
    /*adaptAnalysisTask<McInformation>(cfgc),
    adaptAnalysisTask<McFakeInfo>(cfgc),
    adaptAnalysisTask<MyAnalysisTask>(cfgc)*/
  };
}