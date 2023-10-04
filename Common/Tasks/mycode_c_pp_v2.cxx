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
}

struct McInfomation{
  Produces<aod::MCPair> mcpairtable;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsXY{"nBinsXY", 10010, "N bins in X,Y axis"};
  Configurable<int> nBinsZ{"nBinsZ", 110010, "N bins in Z axis"};
  Configurable<int> nBinsR{"nBinsR", 10010, "N bins in R"};
	Configurable<int> nBinsDist{"nBinsDist", 100010, ""};
	Configurable<int> fwdTrackType{"fwdTrackType", 0, "N TrackType in fwd"};
  Configurable<int> nBinsCut{"nBinsCut", 20010, "N bins in Cut"};

  void init(InitContext const&){
    const AxisSpec axisZ{nBinsZ, -60.0005, 50.0005, "z(cm)"};
    const AxisSpec axisDCAT{100010, -0.0005, 10.0005, "DCAT"};
    const AxisSpec axisCos{1000, 0, 1, ""};
    const AxisSpec axisCutZ{10000, -1.5, 0, "z(cm)"};
    const AxisSpec axisDiff{10000, 0, 10, "(cm)"};
    const AxisSpec axisMult{5000, 0, 500, "multi(cm)"};

    histos.add("SecVtx_Z", "SecVtx_Z", kTH1F, {axisZ});
    histos.add("SecVtx_XY", "SecVtx_XY", kTH1F, {axisDCAT});
    histos.add("PredictZ", "PredictZ", kTH1F, {axisZ});
    histos.add("PredictXY", "PredictXY", kTH1F, {axisDCAT});
    histos.add("DCAT_Muon", "DCAT_Muon", kTH1F, {axisDCAT});
    histos.add("DCAT_Kaon", "DCAT_Kaon", kTH1F, {axisDCAT});
    histos.add("DCAT_Primary", "DCAT_Primary", kTH1F, {axisDCAT});
    histos.add("CosSimilarity_muk", "CosSimilarity_muk", kTH1F, {axisCos});
    histos.add("CosSimilarity_others", "CosSimilarity_others", kTH1F, {axisCos});
    histos.add("Diffmuk_cutz", "Diffmuk_cutz", kTH2F, {axisCutZ, axisDiff});
    histos.add("Add_distmuk", "Add_distmuk", kTH1F, {axisMult});
    histos.add("AddDistmuk_Cossim", "AddDistmuk_Cossim", kTH2F, {axisMult, axisCos});
    histos.add("AddDistmuk_Cossim_cut", "AddDistmuk_Cossim_cut", kTH2F, {axisMult, axisCos});
    histos.add("Add_distmuothers", "Add_distmuothers", kTH1F, {axisMult});
    histos.add("AddDistmuothers_Cossim", "AddDistmuothers_Cossim", kTH2F, {axisMult, axisCos});
    histos.add("AddDistmuothers_Cossim_cut", "AddDistmuothers_Cossim_cut", kTH2F, {axisMult, axisCos});

    //Muon and Kaon 
    const AxisSpec axisXY{10000, -5, 5, "cm"};
    const AxisSpec axisPCAZ{20000, -10, 10, "cm"};
    const AxisSpec axisR{10000, 0, 10, "cm"};
    const AxisSpec axisChi{10001, -0.5, 10000.5, ""};
    histos.add("Exp_muk_x", "Exp_muk_x", kTH1F, {axisXY});
    histos.add("Exp_muk_y", "Exp_muk_y", kTH1F, {axisXY});
    histos.add("Exp_muk_z", "Exp_muk_z", kTH1F, {axisPCAZ});
    histos.add("Exp_muk_r", "Exp_muk_r", kTH1F, {axisR});
    histos.add("Chi2_muk", "Chi2_muk", kTH1F, {axisChi});
    histos.add("Chi2_xyz_add", "Chi2_xyz_add", kTH1F, {axisChi});

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
                int mcmucolid, mckcolid, mucolid, kcolid, mcmuid, mupdg, mckid, kpdg;
                double mudcaX, mudcaY, mumftX, mumftY, mumftZ, kdcaX, kdcaY, kmftX, kmftY, kmftZ, mumcX, mumcY, mumcZ, kmcX, kmcY, kmcZ;
                const double Col_x = collision.posX();
                const double Col_y = collision.posY();
                const double Col_z = collision.posZ();
                vector<double> mux{};
                vector<double> muy{};
                vector<double> kx{};
                vector<double> ky{};
                
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
                        dc++;
                        for(int i=0; i<=10000; i++){
                          auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                          if((cut_z-Col_z)>0) break;
                          mftpars1.propagateToZlinear(cut_z);
                          mux.push_back(mftpars1.getX());
                          muy.push_back(mftpars1.getY());
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
                        fc++;
                        for(int i=0; i<=10000; i++){
                          auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                          if((cut_z-Col_z)>0) break;
                          mftpars2.propagateToZlinear(cut_z);
                          kx.push_back(mftpars2.getX());
                          ky.push_back(mftpars2.getY());
                        }                       
                      }
                    }
                  }
                }
                if(mucolid==kcolid && dc+fc==2){
                  if(mumcZ==kmcZ && collision.mcCollisionId()==mcmucolid){
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
                    auto pcaX = ((predict_mux+predict_kx)/2)-Col_x;
                    auto pcaY = ((predict_muy+predict_ky)/2)-Col_y;
                    auto pcaZ = ((predict_muz+predict_kz)/2)-Col_z;
                    double pre_r = sqrt(pow(predict_mux-predict_kx,2)+pow(predict_muy-predict_ky,2)+pow(predict_muz-predict_kz,2));
                    //LOGF(info, "MC SecVtx: %lf, %lf, %lf   PCA: %lf, %lf, %lf", mumcX,mumcY,mumcZ,pcaX,pcaY,pcaZ);
                    histos.fill(HIST("SecVtx_Z"), mumcZ-collision.mcCollision().posZ());
                    histos.fill(HIST("SecVtx_XY"), sqrt(pow(mumcX-collision.mcCollision().posX(),2)+pow(mumcY-collision.mcCollision().posY(),2)));
                    histos.fill(HIST("PredictZ"), pcaZ);
                    histos.fill(HIST("PredictXY"), sqrt(pow(pcaX-Col_x,2)+pow(pcaY-Col_y,2)));
                    histos.fill(HIST("DCAT_Muon"), sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2)));
                    histos.fill(HIST("DCAT_Kaon"), sqrt(pow(kdcaX-Col_x,2)+pow(kdcaY-Col_y,2)));

                    histos.fill(HIST("Exp_muk_x"), pcaX);
                    histos.fill(HIST("Exp_muk_y"), pcaY);
                    histos.fill(HIST("Exp_muk_z"), pcaZ);
                    histos.fill(HIST("Exp_muk_r"), pre_r);
                    
                    //cos類似度分布
                    auto vecx_mu = mumftX - mudcaX;
                    auto vecy_mu = mumftY - mudcaY;
                    auto vecz_mu = mumftZ - Col_z;
                    auto vecx_k = kmftX - kdcaX;
                    auto vecy_k = kmftY - kdcaY;
                    auto vecz_k = kmftZ - Col_z;
                    auto cosxy = (vecx_mu*vecx_k + vecy_mu*vecy_k + vecz_mu*vecz_k)/((sqrt(pow(vecx_mu,2)+pow(vecy_mu,2)+pow(vecz_mu,2)))*(sqrt(pow(vecx_k,2)+pow(vecy_k,2)+pow(vecz_k,2))));
                    histos.fill(HIST("CosSimilarity_muk"), 1-cosxy);
                    mcpairtable(mucolid,mcmuid,mupdg,mckid,kpdg,kdcaX-Col_x,kdcaY-Col_y,Col_x,Col_y,Col_z,mumcX,mumcY,mumcZ,pcaX-Col_x,pcaY-Col_y,pcaZ-Col_z,cosxy);
                    //LOGF(info, "Collision ID: %d", mucolid);
                    //LOGF(info, "%lf,%lf,%lf,%lf,%lf,%lf", mumftX,mumftY,mumftZ,mudcaX,mudcaY,Col_z+5);
                    //LOGF(info, "%lf,%lf,%lf,%lf,%lf,%lf", kmftX,kmftY,kmftZ,kdcaX,kdcaY,Col_z+5);

                    int nCan = 0;
                    for(auto& mfttrack : mfttracks){
                      if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
                      if(mfttrack.mcParticleId()==mcmuid || mfttrack.mcParticleId()==mckid) continue;
                      if(find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack.globalIndex()) != ambTrackIds.end()) continue;
                      if(mfttrack.collisionId()==mucolid){
                        if(mfttrack.mcParticleId()!=mcmuid){
                          if(mfttrack.sign()!=0) nCan++;
                        }
                        double mftchi2 = mfttrack.chi2();
                        SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                        vector<double> mftv1;
                        SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                        o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
                        mftpars1.propagateToZlinear(Col_z);
                        if(mfttrack.mcParticle().producedByGenerator()==true){
                          histos.fill(HIST("DCAT_Primary"), sqrt(pow(mftpars1.getX()-Col_x,2)+pow(mftpars1.getY()-Col_y,2)));
                        }
                        //LOGF(info, "%lf,%lf,%lf,%lf,%lf,%lf", mfttrack.x(),mfttrack.y(),mfttrack.z(),mftpars1.getX(),mftpars1.getY(),Col_z+5);
                        auto vecx_can = mfttrack.x() - mftpars1.getX();
                        auto vecy_can = mfttrack.y() - mftpars1.getY();
                        auto vecz_can = mfttrack.z() - Col_z;
                        auto cosxy_can = (vecx_mu*vecx_can + vecy_mu*vecy_can + vecz_mu*vecz_can)/((sqrt(pow(vecx_mu,2)+pow(vecy_mu,2)+pow(vecz_mu,2)))*(sqrt(pow(vecx_can,2)+pow(vecy_can,2)+pow(vecz_can,2))));
                        histos.fill(HIST("CosSimilarity_others"), 1-cosxy_can);
                        
                        double addit = 0;
                        for(int i=0; i<mux.size(); i++){
                          auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                          mftpars1.propagateToZlinear(cut_z);
                          auto diff = sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2));
                          addit += diff;
                        }
                        histos.fill(HIST("Add_distmuothers"), addit);
                        histos.fill(HIST("AddDistmuothers_Cossim"), addit, cosxy_can);
                        if(addit<25 && cosxy_can>0.99){
                          histos.fill(HIST("AddDistmuothers_Cossim_cut"), addit, cosxy_can);
                        }
                      }
                    }

                    //x-y平面距離
                    double multi = 0;
                    if(mux.size()==muy.size()){
                      if(kx.size()==ky.size()){
                        if(mux.size()==kx.size()){
                          for(int i=0; i<mux.size(); i++){
                            auto cut_z = (-0.2236-3*0.3081)+0.005*i;
                            auto diff = sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2));
                            histos.fill(HIST("Diffmuk_cutz"), cut_z, diff);
                            multi += diff;
                          }
                        }
                      }
                    }
                    histos.fill(HIST("Add_distmuk"), multi);
                    histos.fill(HIST("AddDistmuk_Cossim"), multi, cosxy);
                    if(multi<25 && cosxy>0.99){
                      histos.fill(HIST("AddDistmuk_Cossim_cut"), multi, cosxy);
                    }

                    //χ二乗適合度検定
                    auto chi_x = pcaX;//10000*pcaX;
                    auto chi_y = pcaY;//10000*pcaY;
                    auto chi_z = -pcaZ;//-10000*pcaZ;
                    auto chi_r = pre_r;//10000*pre_r;
                    auto chi_cos = cosxy;//10*cosxy;
                    //LOGF(info, "Chi2X: %lf, Chi2Y: %lf, Chi2Z: %lf, Chi2R: %lf, Chi2Cos: %lf, Chi2Add: %lf",chi_x,chi_y,chi_z,chi_r,chi_cos,multi);
                    const double mean_xy = 0.003455+0.0003302;//34.55+3.302;
                    const double mean_z = 0.2298;//2298;
                    const double mean_r = 0.02856;//285.6;
                    const double mean_cos = 0.9929;//9.929;
                    const double mean_add = 16.57;

                    double chi_xy = chi_x+chi_y;
                    //LOGF(info, "XY: %lf, Z: %lf, R: %lf, cos: %lf, Add: %lf",pow(chi_xy-mean_xy,2)/mean_xy,pow(chi_z-mean_z,2)/mean_z,pow(chi_r-mean_r,2)/mean_r,pow(chi_cos-mean_cos,2)/mean_cos,pow(multi-mean_add,2)/mean_add);
                    double chi2_muk = (pow(chi_xy-mean_xy,2)/mean_xy)+(pow(chi_z-mean_z,2)/mean_z)+(pow(chi_r-mean_r,2)/mean_r)+(pow(chi_cos-mean_cos,2)/mean_cos)+(pow(multi-mean_add,2)/mean_add);
                    histos.fill(HIST("Chi2_muk"), chi2_muk);
                    double chi2_xyz = (pow(chi_xy-mean_xy,2)/mean_xy)+(pow(chi_z-mean_z,2)/mean_z)+(pow(multi-mean_add,2)/mean_add);
                    histos.fill(HIST("Chi2_xyz_add"), chi2_xyz);
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

struct McInfomation2{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&){
    const AxisSpec axisX{10000, -25, 25, "x(cm)"};
    const AxisSpec axisY{10000, -25, 25, "y(cm)"};
    const AxisSpec axisR{10000, 0, 25, "r(cm)"};
    const AxisSpec axisRatio{500010, -0.0005, 50.0005, ""};
    const AxisSpec axisSumDist{100010, -0.0005, 100.0005, ""};
    const AxisSpec axisTheta{4000, -2, 2, "theta(rad)"};
    /*histos.add("z1_dist_muk", "z1_dist_muk", kTH2F, {axisX, axisY});
    histos.add("z1_dist_bg", "z1_dist_bg", kTH2F, {axisX, axisY});
    histos.add("z2_dist_muk", "z2_dist_muk", kTH2F, {axisX, axisY});
    histos.add("z2_dist_bg", "z2_dist_bg", kTH2F, {axisX, axisY});
    histos.add("z3_dist_muk", "z3_dist_muk", kTH2F, {axisX, axisY});
    histos.add("z3_dist_bg", "z3_dist_bg", kTH2F, {axisX, axisY});
    histos.add("z4_dist_muk", "z4_dist_muk", kTH2F, {axisX, axisY});
    histos.add("z4_dist_bg", "z4_dist_bg", kTH2F, {axisX, axisY});*/
    /*histos.add("z5_dist_muk", "z5_dist_muk", kTH2F, {axisX, axisY});
    histos.add("z5_dist_bg", "z5_dist_bg", kTH2F, {axisX, axisY});*/
    histos.add("z1_dist_muk", "z1_dist_muk", kTH1F, {axisR});
    histos.add("z1_dist_basemuk", "z1_dist_basemuk", kTH1F, {axisR});
    histos.add("z1_dist_bg", "z1_dist_bg", kTH1F, {axisR});
    histos.add("z1_dist_basebg", "z1_dist_basebg", kTH1F, {axisR});
    histos.add("z2_dist_muk", "z2_dist_muk", kTH1F, {axisR});
    histos.add("z2_dist_basemuk", "z2_dist_basemuk", kTH1F, {axisR});
    histos.add("z2_dist_bg", "z2_dist_bg", kTH1F, {axisR});
    histos.add("z2_dist_basebg", "z2_dist_basebg", kTH1F, {axisR});
    histos.add("z3_dist_muk", "z3_dist_muk", kTH1F, {axisR});
    histos.add("z3_dist_bg", "z3_dist_bg", kTH1F, {axisR});
    histos.add("z4_dist_muk", "z4_dist_muk", kTH1F, {axisR});
    histos.add("z4_dist_bg", "z4_dist_bg", kTH1F, {axisR});
    histos.add("SumRatio_prelay_bg", "SumRatio_prelay_bg", kTH1F, {axisRatio});
    histos.add("SumDist_bg", "SumDist_bg", kTH1F, {axisSumDist});
    histos.add("SumRatio_prelay_muk", "SumRatio_prelay_muk", kTH1F, {axisRatio});
    histos.add("SumDist_muk", "SumDist_muk", kTH1F, {axisRatio});
    histos.add("ThetaDistribution_muk", "ThetaDistribution_muk", kTH1F, {axisTheta});
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
                int mcmucolid, mckcolid, mucolid, kcolid, mcmuid, mupdg, mckid, kpdg;
                double mudcaX, mudcaY, mumftX, mumftY, mumftZ, kdcaX, kdcaY, kmftX, kmftY, kmftZ, mumcX, mumcY, mumcZ, kmcX, kmcY, kmcZ;
                const double Col_z = collision.posZ();
                vector<double> mux{};
                vector<double> muy{};
                vector<double> muz{};
                vector<double> kx{};
                vector<double> ky{};
                vector<double> kz{};
                
                for(auto& mfttrack : mfttracks){
                  if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
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
                        dc++;
                        
                        for(int i=1; i<=5; i++){
                          auto propaZ = Col_z - i;
                          mftpars1.propagateToZlinear(propaZ);
                          mux.push_back(mftpars1.getX());
                          muy.push_back(mftpars1.getY());
                          muz.push_back(propaZ);
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
                        fc++;

                        for(int i=1; i<=5; i++){
                          auto propaZ = Col_z - i;
                          mftpars2.propagateToZlinear(propaZ);
                          kx.push_back(mftpars2.getX());
                          ky.push_back(mftpars2.getY());
                          kz.push_back(propaZ);
                        }
                      }
                    }
                  }
                }
                if(mucolid==kcolid && dc+fc==2){
                  auto mu_unit = sqrt(pow(mumftX-mudcaX,2)+pow(mumftY-mudcaY,2)+pow(mumftZ-Col_z,2));
                  auto nx = (mumftX-mudcaX)/mu_unit;
                  auto ny = (mumftY-mudcaY)/mu_unit;
                  auto nz = (mumftZ-Col_z)/mu_unit;
                  auto k_unit = sqrt(pow(kmftX-kdcaX,2)+pow(kmftY-kdcaY,2)+pow(kmftZ-Col_z,2));
                  auto mx = (kmftX-kdcaX)/k_unit;
                  auto my = (kmftY-kdcaY)/k_unit;
                  auto mz = (kmftZ-Col_z)/k_unit;

                  if(mumcZ==kmcZ && collision.mcCollisionId()==mcmucolid){
                    for(auto& mfttrack : mfttracks){
                      if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
                      if(mfttrack.mcParticleId()==mcmuid || mfttrack.mcParticleId()==mckid) continue;
                      if(find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack.globalIndex()) != ambTrackIds.end()) continue;
                      if(mfttrack.collisionId()==mucolid){
                        double mftchi2 = mfttrack.chi2();
                        SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                        vector<double> mftv1;
                        SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                        o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};

                        double first_lay = 0;
                        double sum_ratio = 0;
                        double sum_distance = 0;
                        mftpars1.propagateToZlinear(Col_z);
                        double ldcaX = mftpars1.getX();
                        double ldcaY = mftpars1.getY();
                        auto l_unit = sqrt(pow(mfttrack.x()-mftpars1.getX(),2)+pow(mfttrack.y()-mftpars1.getY(),2)+pow(mfttrack.z()-Col_z,2));
                        auto lx = (mfttrack.x()-mftpars1.getX())/l_unit;
                        auto ly = (mfttrack.y()-mftpars1.getY())/l_unit;
                        auto lz = (mfttrack.z()-Col_z)/l_unit;

                        for(int i=0; i<=4; i++){
                          mftpars1.propagateToZlinear(Col_z-(i+1));
                          if(i==0){
                            //histos.fill(HIST("z1_dist_bg"), mftpars1.getX()-mux[i], mftpars1.getY()-muy[i]);
                            first_lay = sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2));
                            sum_distance += first_lay; 
                            histos.fill(HIST("z1_dist_bg"), first_lay);
                            double s = ((ldcaX-mux[i])*nx+(ldcaY-muy[i])*ny+i*muz[i])/(nx*lx+ny*ly+nz*lz);
                            double dist_basemuk = sqrt(pow((ldcaX+s*lx)-mux[i],2)+pow((ldcaY+s*ly)-muy[i],2)+pow((Col_z+s*lz)-muz[i],2));
                            histos.fill(HIST("z1_dist_basebg"), dist_basemuk);
                          }else if(i==1){
                            //histos.fill(HIST("z2_dist_bg"), mftpars1.getX()-mux[i], mftpars1.getY()-muy[i]);
                            if(first_lay==0){
                              first_lay = sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2));
                              histos.fill(HIST("z2_dist_bg"), first_lay);
                            }else{
                              sum_ratio += sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2))/first_lay;
                              first_lay = sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2));
                              histos.fill(HIST("z2_dist_bg"), first_lay);
                            }
                            sum_distance += first_lay;
                          }else if(i==2){
                            //histos.fill(HIST("z3_dist_bg"), mftpars1.getX()-mux[i], mftpars1.getY()-muy[i]);
                            if(first_lay==0){
                              first_lay = sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2));
                              histos.fill(HIST("z3_dist_bg"), first_lay);
                            }else{
                              sum_ratio += sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2))/first_lay;
                              first_lay = sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2));
                              histos.fill(HIST("z3_dist_bg"), first_lay);
                            }
                            sum_distance += first_lay;
                          }else if(i==3){
                            //histos.fill(HIST("z4_dist_bg"), mftpars1.getX()-mux[i], mftpars1.getY()-muy[i]);
                            if(first_lay==0){
                              first_lay = sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2));
                              histos.fill(HIST("z4_dist_bg"), first_lay);
                            }else{
                              sum_ratio += sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2))/first_lay;
                              first_lay = sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2));
                              histos.fill(HIST("z4_dist_bg"), first_lay);
                            }
                            sum_distance += first_lay;
                          }else if(i==4){
                            //histos.fill(HIST("z5_dist_bg"), mftpars1.getX()-mux[i], mftpars1.getY()-muy[i]);
                            if(first_lay==0){
                              first_lay = sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2));
                            }else{
                              sum_ratio += sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2))/first_lay;
                              first_lay = sqrt(pow(mftpars1.getX()-mux[i],2)+pow(mftpars1.getY()-muy[i],2));
                            }
                            sum_distance += first_lay;
                          }
                        }
                        histos.fill(HIST("SumRatio_prelay_bg"), sum_ratio);
                        histos.fill(HIST("SumDist_bg"), sum_distance);
                      }
                    }
                    double first_lay_muk = 0;
                    double sum_ratio_muk = 0;
                    double sum_distance_muk = 0;
                    double first_theta = 0;
                    double diff_theta = 0;
                    for(int i=0; i<=4; i++){
                      if(i==0){
                        //histos.fill(HIST("z1_dist_muk"), kx[i]-mux[i], ky[i]-muy[i]);
                        first_lay_muk = sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2));
                        sum_distance_muk += first_lay_muk;
                        histos.fill(HIST("z1_dist_muk"), first_lay_muk);
                        first_theta = atan2((ky[i]-muy[i]),(kx[i]-mux[i]));

                        double s = ((kdcaX-mux[i])*nx+(kdcaY-muy[i])*ny+i*muz[i])/(nx*mx+ny*my+nz*mz);
                        double dist_basemuk = sqrt(pow((kdcaX+s*mx)-mux[i],2)+pow((kdcaY+s*my)-muy[i],2)+pow((Col_z+s*mz)-muz[i],2));
                        histos.fill(HIST("z1_dist_basemuk"), dist_basemuk);

                      }else if(i==1){
                        //histos.fill(HIST("z2_dist_muk"), kx[i]-mux[i], ky[i]-muy[i]);
                        if(first_lay_muk==0){
                          first_lay_muk = sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2));
                          histos.fill(HIST("z2_dist_muk"), first_lay_muk);
                        }else{
                          sum_ratio_muk += sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2))/first_lay_muk;
                          first_lay_muk = sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2));
                          histos.fill(HIST("z2_dist_muk"), first_lay_muk);
                        }
                        sum_distance_muk += first_lay_muk;
                        diff_theta += first_theta - atan2((ky[i]-muy[i]),(kx[i]-mux[i]));
                        first_theta = atan2((ky[i]-muy[i]),(kx[i]-mux[i]));
                        
                        double s = ((kdcaX-mux[i])*nx + (kdcaY-muy[i])*ny + (Col_z-muz[i])*nz)/(nx*mx+ny*my+nz*mz);
                        double dist_basemuk = sqrt(pow(((kdcaX+s*mx)-mux[i]),2)+pow(((kdcaY+s*my)-muy[i]),2)+pow(((Col_z+s*mz)-muz[i]),2));
                        histos.fill(HIST("z2_dist_basemuk"), dist_basemuk);

                      }else if(i==2){
                        //histos.fill(HIST("z3_dist_muk"), kx[i]-mux[i], ky[i]-muy[i]);
                        if(first_lay_muk==0){
                          first_lay_muk = sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2));
                          histos.fill(HIST("z3_dist_muk"), first_lay_muk);
                        }else{
                          sum_ratio_muk += sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2))/first_lay_muk;
                          first_lay_muk = sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2));
                          histos.fill(HIST("z3_dist_muk"), first_lay_muk);
                        }
                        sum_distance_muk += first_lay_muk;
                        diff_theta += first_theta - atan2((ky[i]-muy[i]),(kx[i]-mux[i]));
                        first_theta = atan2((ky[i]-muy[i]),(kx[i]-mux[i]));
                      }else if(i==3){
                        //histos.fill(HIST("z4_dist_muk"), kx[i]-mux[i], ky[i]-muy[i]);
                        if(first_lay_muk==0){
                          first_lay_muk = sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2));
                          histos.fill(HIST("z4_dist_muk"), first_lay_muk);
                        }else{
                          sum_ratio_muk += sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2))/first_lay_muk;
                          first_lay_muk = sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2));
                          histos.fill(HIST("z4_dist_muk"), first_lay_muk);
                        }
                        sum_distance_muk += first_lay_muk;
                        diff_theta += first_theta - atan2((ky[i]-muy[i]),(kx[i]-mux[i]));
                        first_theta = atan2((ky[i]-muy[i]),(kx[i]-mux[i]));
                      }else if(i==4){
                        //histos.fill(HIST("z5_dist_muk"), kx[i]-mux[i], ky[i]-muy[i]);
                        if(first_lay_muk==0){
                          first_lay_muk = sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2));
                        }else{
                          sum_ratio_muk += sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2))/first_lay_muk;
                          first_lay_muk = sqrt(pow(kx[i]-mux[i],2)+pow(ky[i]-muy[i],2));
                        }
                        sum_distance_muk += first_lay_muk;
                        diff_theta += first_theta - atan2((ky[i]-muy[i]),(kx[i]-mux[i]));
                        first_theta = atan2((ky[i]-muy[i]),(kx[i]-mux[i]));
                      }
                    }
                    histos.fill(HIST("SumRatio_prelay_muk"), sum_ratio_muk);
                    histos.fill(HIST("SumDist_muk"), sum_distance_muk);
                    histos.fill(HIST("ThetaDistribution_muk"), diff_theta);
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

struct MyAnalysisTask{
  //やること → 普通にベクトル法で最近接粒子を決める。その粒子のcos類似度を確認
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&){
    const AxisSpec axisCos{1000, 0, 1, ""};
    const AxisSpec axisMult{5000, 0, 500, "total(cm)"};
    const AxisSpec axisPDG{21, -0.5, 20.5, ""};
    const AxisSpec axisChi{10001, -0.5, 10000.5};
    histos.add("CosSimilarity_true", "CosSimilarity_true", kTH1F, {axisCos});
    histos.add("CosSimilarity_false", "CosSimilarity_false", kTH1F, {axisCos});
    histos.add("CosSimilarity_runk2", "CosSimilarity_runk2", kTH1F, {axisCos});
    histos.add("CosSimilarity_runk3", "CosSimilarity_runk3", kTH1F, {axisCos});
    histos.add("CosSimilarity_runk4", "CosSimilarity_runk4", kTH1F, {axisCos});
    histos.add("CosSimilarity_runk5", "CosSimilarity_runk5", kTH1F, {axisCos});
    histos.add("CosSimilarity_runk6", "CosSimilarity_runk6", kTH1F, {axisCos});
    histos.add("CosSimilarity_runk7", "CosSimilarity_runk7", kTH1F, {axisCos});
    histos.add("CosSimilarity_runk8", "CosSimilarity_runk8", kTH1F, {axisCos});
    histos.add("CosSimilarity_runk9", "CosSimilarity_runk9", kTH1F, {axisCos});
    histos.add("CosSimilarity_runk10", "CosSimilarity_runk10", kTH1F, {axisCos});

    histos.add("TotalDiff_True", "TotalDiff_True", kTH1F, {axisMult});
    histos.add("TotalDiff_False", "TotalDiff_False", kTH1F, {axisMult});
    histos.add("PDGCode", "PDGCode", kTH1F, {axisPDG});
    histos.add("Chi2_bg", "Chi2_bg", kTH1F, {axisChi});
    histos.add("Chi2_True", "Chi2_True", kTH1F, {axisChi});
    histos.add("Chi2_False", "Chi2_False", kTH1F, {axisChi});

    auto mumompdg_d = histos.get<TH1>(HIST("PDGCode"));
    auto* x3 = mumompdg_d->GetXaxis();
    x3->SetBinLabel(1, "Primary"); // 1-37
    x3->SetBinLabel(2, "Pion0"); // 111
    x3->SetBinLabel(3, "Pion±"); // 211
    x3->SetBinLabel(4, "pho"); // 113
    x3->SetBinLabel(5, "K_l"); // 130
    x3->SetBinLabel(6, "Eta"); // 221
    x3->SetBinLabel(7, "Omega"); // 223
    x3->SetBinLabel(8, "K_s"); // 310
    x3->SetBinLabel(9, "K*0(892)"); // 313
    x3->SetBinLabel(10, "K±"); // 321
    x3->SetBinLabel(11, "K*±(892)"); // 323
    x3->SetBinLabel(12, "Eta_prim"); // 331
    x3->SetBinLabel(13, "Phi"); // 333
    x3->SetBinLabel(14, "D±"); // 411
    x3->SetBinLabel(15, "D*±"); // 413
    x3->SetBinLabel(16, "D0"); // 421
    x3->SetBinLabel(17, "D_s±"); // 431
    x3->SetBinLabel(18, "B0"); // 511
    x3->SetBinLabel(19, "B+"); // 521
    x3->SetBinLabel(20, "B0_s"); // 531
    x3->SetBinLabel(21, "Baryon"); // 1000 -
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
              const double Col_x = collision.posX();
              const double Col_y = collision.posY();
              const double Col_z = collision.posZ();

              double mudcaX, mudcaY, mumftX, mumftY, mumftZ;

              // For muon
              for(auto& muontrack : mfttracks){
                if(!muontrack.has_collision() || !muontrack.has_mcParticle()) continue;
                if(muontrack.collisionId()==fwdColId && muontrack.mcParticleId()==truepair.muid()){
                  //auto mcParticle_muon = muontrack.mcParticle();
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

                  /*double dist_1 = 1000;
                  double dist_2 = 1000;
                  double dist_3 = 1000;
                  double dist_4 = 1000;
                  double dist_5 = 1000;*/
                  double candcaX, candcaY, canmftX, canmftY, canmftZ;
                  /*double cos_sim_1, cos_sim_2, cos_sim_3, cos_sim_4, cos_sim_5;
                  int canid_1, canpdg_1, canid_2, canpdg_2, canid_3, canpdg_3, canid_4, canpdg_4, canid_5, canpdg_5;
                  double canmftparx_1, canmftpary_1, canmftparz_1, canmftparphi_1, canmftpartgl_1, canmftpartpt_1, canmftchi2_1;
                  double canmftparx_2, canmftpary_2, canmftparz_2, canmftparphi_2, canmftpartgl_2, canmftpartpt_2, canmftchi2_2;
                  double canmftparx_3, canmftpary_3, canmftparz_3, canmftparphi_3, canmftpartgl_3, canmftpartpt_3, canmftchi2_3;
                  double canmftparx_4, canmftpary_4, canmftparz_4, canmftparphi_4, canmftpartgl_4, canmftpartpt_4, canmftchi2_4;
                  double canmftparx_5, canmftpary_5, canmftparz_5, canmftparphi_5, canmftpartgl_5, canmftpartpt_5, canmftchi2_5;
                  double chi2_1, chi2_2, chi2_3, chi2_4, chi2_5;*/
                  vector<vector<double>> data;
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

                      auto chi_x = pcaX;//10000*pcaX;
                      auto chi_y = pcaY;//10000*pcaY;
                      auto chi_z = -pcaZ;//-10000*pcaZ;
                      auto chi_r = r_xyz;//10000*pre_r;
                      auto chi_cos = cosxy;//10*cosxy;
                      const double mean_xy = 0.003455+0.0003302;//34.55+3.302;
                      const double mean_z = 0.2298;//2298;
                      const double mean_r = 0.02856;//285.6;
                      const double mean_cos = 0.9929;//9.929;
                      const double mean_add = 16.57;
                      double chi_xy = chi_x+chi_y;

                      double add_can = 0;

                      for(int i=0; i<5000; i++){
                        auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                        if((cut_z-Col_z)>0) break;
                        mftpars1.propagateToZlinear(cut_z);
                        mftpars2.propagateToZlinear(cut_z);
                        auto diff = sqrt(pow(mftpars2.getX()-mftpars1.getX(),2)+pow(mftpars2.getY()-mftpars1.getY(),2));
                        add_can += diff;
                      }
                      double chi2_bg = (pow(chi_xy-mean_xy,2)/mean_xy)+(pow(chi_z-mean_z,2)/mean_z)+(pow(chi_r-mean_r,2)/mean_r)+(pow(chi_cos-mean_cos,2)/mean_cos)+(pow(add_can-mean_add,2)/mean_add);

                      //LOGF(info, "XY: %lf, Z: %lf, R: %lf, cos: %lf, Add: %lf",pow(chi_xy-mean_xy,2)/mean_xy,pow(chi_z-mean_z,2)/mean_z,pow(chi_r-mean_r,2)/mean_r,pow(chi_cos-mean_cos,2)/mean_cos,pow(multi-mean_add,2)/mean_add);
                      if(mfttrack.mcParticleId()!=truepair.kid() && mcParticle_mft.producedByGenerator()==true){
                        histos.fill(HIST("Chi2_bg"), chi2_bg);
                      }

                      //if(cosxy<0.99) continue;

                      /*double add_can = 0;

                      for(int i=0; i<5000; i++){
                        auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                        if((cut_z-Col_z)>0) break;
                        mftpars1.propagateToZlinear(cut_z);
                        mftpars2.propagateToZlinear(cut_z);
                        auto diff = sqrt(pow(mftpars2.getX()-mftpars1.getX(),2)+pow(mftpars2.getY()-mftpars1.getY(),2));
                        add_can += diff;
                      }*/

                      //if(add_can>25) continue;
                      double mcid = mfttrack.mcParticleId();
                      vector<double> pardata = {r_xyz, cosxy, mcid, fabs(mcParticle_mft.pdgCode()), canmftX, canmftY, canmftZ, mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt(), mftchi2, chi2_bg};
                      data.push_back(pardata);
                      /*if(dist_1>=r_xyz){
                        dist_1 = r_xyz;
                        cos_sim_1 = cosxy;
                        canid_1 = mfttrack.mcParticleId();
                        canpdg_1 = fabs(mcParticle_mft.pdgCode());
                        canmftparx_1 = canmftX;
                        canmftpary_1 = canmftY;
                        canmftparz_1 = canmftZ;
                        canmftparphi_1 = mfttrack.phi();
                        canmftpartgl_1 = mfttrack.tgl();
                        canmftpartpt_1 = mfttrack.signed1Pt();
                        canmftchi2_1 = mftchi2;
                        chi2_1 = chi2_bg;
                      }else if(dist_1<r_xyz && dist_2>=r_xyz){
                        dist_2 = r_xyz;
                        cos_sim_2 = cosxy;
                        canid_2 = mfttrack.mcParticleId();
                        canpdg_2 = fabs(mcParticle_mft.pdgCode());
                        canmftparx_2 = canmftX;
                        canmftpary_2 = canmftY;
                        canmftparz_2 = canmftZ;
                        canmftparphi_2 = mfttrack.phi();
                        canmftpartgl_2 = mfttrack.tgl();
                        canmftpartpt_2 = mfttrack.signed1Pt();
                        canmftchi2_2 = mftchi2;
                        chi2_2 = chi2_bg;
                      }else if(dist_2<r_xyz && dist_3>=r_xyz){
                        dist_3 = r_xyz;
                        cos_sim_3 = cosxy;
                        canid_3 = mfttrack.mcParticleId();
                        canpdg_3 = fabs(mcParticle_mft.pdgCode());
                        canmftparx_3 = canmftX;
                        canmftpary_3 = canmftY;
                        canmftparz_3 = canmftZ;
                        canmftparphi_3 = mfttrack.phi();
                        canmftpartgl_3 = mfttrack.tgl();
                        canmftpartpt_3 = mfttrack.signed1Pt();
                        canmftchi2_3 = mftchi2;
                        chi2_3 = chi2_bg;
                      }else if(dist_3<r_xyz && dist_4>=r_xyz){
                        dist_4 = r_xyz;
                        cos_sim_4 = cosxy;
                        canid_4 = mfttrack.mcParticleId();
                        canpdg_4 = fabs(mcParticle_mft.pdgCode());
                        canmftparx_4 = canmftX;
                        canmftpary_4 = canmftY;
                        canmftparz_4 = canmftZ;
                        canmftparphi_4 = mfttrack.phi();
                        canmftpartgl_4 = mfttrack.tgl();
                        canmftpartpt_4 = mfttrack.signed1Pt();
                        canmftchi2_4 = mftchi2;
                        chi2_4 = chi2_bg;
                      }else if(dist_4<r_xyz && dist_5>=r_xyz){
                        dist_5 = r_xyz;
                        cos_sim_5 = cosxy;
                        canid_5 = mfttrack.mcParticleId();
                        canpdg_5 = fabs(mcParticle_mft.pdgCode());
                        canmftparx_5 = canmftX;
                        canmftpary_5 = canmftY;
                        canmftparz_5 = canmftZ;
                        canmftparphi_5 = mfttrack.phi();
                        canmftpartgl_5 = mfttrack.tgl();
                        canmftpartpt_5 = mfttrack.signed1Pt();
                        canmftchi2_5 = mftchi2;
                        chi2_5 = chi2_bg;
                      }*/

                    }
                  }
                  
                  sort(data.begin(), data.end());
                  /*
                  int counter = 0;
                  for(int i=0; i<data.size(); i++){
                    for(int j=0; j<data.at(0).size(); j++){
                      if(counter==0){
                        cout << "data = { " << data.at(i).at(j) << ", ";
                        counter++;
                      }else if(j==(data.at(0).size()-1)){
                        if(i==(data.size()-1)){
                          cout << data.at(i).at(j) << " }" << endl;
                        }else{
                          cout << data.at(i).at(j) << " }, " << endl << "{ ";
                        }
                      }
                      else{
                        cout << data.at(i).at(j) << ", ";
                      }
                    }
                  }*/
                  //LOGF(info, "Vector Size: %d", data.size());
                  SMatrix5 canpars(data.at(0).at(4), data.at(0).at(5), data.at(0).at(7), data.at(0).at(8), data.at(0).at(9));
                  vector<double> mftv3;
                  SMatrix55 mftcovs3(mftv3.begin(), mftv3.end());
                  o2::track::TrackParCovFwd mftpars3{data.at(0).at(6), canpars, mftcovs3, data.at(0).at(10)};
                  double add_true=0;
                  //double add_false=0;
                  //LOGF(info, "TruePairID: %d  Chi2&ID: 1:%lf,%d  2:%lf,%d  3:%lf,%d  4:%lf,%d  5:%lf,%d",truepair.kid(),chi2_1,canid_1,chi2_2,canid_2,chi2_3,canid_3,chi2_4,canid_4,chi2_5,canid_5);
                  if(data.at(0).at(2)==truepair.kid()){
                    //LOGF(info, "CollisionID: %d, correct!", fwdColId);
                    histos.fill(HIST("CosSimilarity_true"), data.at(0).at(1));
                    for(int i=0; i<5000; i++){
                      auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                      if((cut_z-Col_z)>0) break;
                      mftpars1.propagateToZlinear(cut_z);
                      mftpars3.propagateToZlinear(cut_z);
                      auto diff = sqrt(pow(mftpars3.getX()-mftpars1.getX(),2)+pow(mftpars3.getY()-mftpars1.getY(),2));
                      add_true += diff;
                    }
                    histos.fill(HIST("TotalDiff_True"), add_true);
                    histos.fill(HIST("Chi2_True"), data.at(0).at(11));
                  }/*else{
                    histos.fill(HIST("CosSimilarity_false"), cos_sim_1);
                    for(int i=0; i<5000; i++){
                      auto cut_z = Col_z + (-0.2236-3*0.3081)+0.005*i;
                      if((cut_z-Col_z)>0) break;
                      mftpars1.propagateToZlinear(cut_z);
                      mftpars3.propagateToZlinear(cut_z);
                      auto diff = sqrt(pow(mftpars3.getX()-mftpars1.getX(),2)+pow(mftpars3.getY()-mftpars1.getY(),2));
                      add_false += diff;
                    }
                    histos.fill(HIST("TotalDiff_False"), add_false);
                    if(fabs(canpdg_1)<38) histos.fill(HIST("PDGCode"), 0);
                    if(fabs(canpdg_1)==111) histos.fill(HIST("PDGCode"), 1);
                    if(fabs(canpdg_1)==211) histos.fill(HIST("PDGCode"), 2);
                    if(fabs(canpdg_1)==113) histos.fill(HIST("PDGCode"), 3);
                    if(fabs(canpdg_1)==130) histos.fill(HIST("PDGCode"), 4);
                    if(fabs(canpdg_1)==221) histos.fill(HIST("PDGCode"), 5);
                    if(fabs(canpdg_1)==223) histos.fill(HIST("PDGCode"), 6);
                    if(fabs(canpdg_1)==310) histos.fill(HIST("PDGCode"), 7);
                    if(fabs(canpdg_1)==313) histos.fill(HIST("PDGCode"), 8);
                    if(fabs(canpdg_1)==321) histos.fill(HIST("PDGCode"), 9);
                    if(fabs(canpdg_1)==323) histos.fill(HIST("PDGCode"), 10);
                    if(fabs(canpdg_1)==331) histos.fill(HIST("PDGCode"), 11);
                    if(fabs(canpdg_1)==333) histos.fill(HIST("PDGCode"), 12);
                    if(fabs(canpdg_1)==411) histos.fill(HIST("PDGCode"), 13);
                    if(fabs(canpdg_1)==413) histos.fill(HIST("PDGCode"), 14);
                    if(fabs(canpdg_1)==421) histos.fill(HIST("PDGCode"), 15);
                    if(fabs(canpdg_1)==431) histos.fill(HIST("PDGCode"), 16);
                    if(fabs(canpdg_1)==511) histos.fill(HIST("PDGCode"), 17);
                    if(fabs(canpdg_1)==521) histos.fill(HIST("PDGCode"), 18);
                    if(fabs(canpdg_1)==531) histos.fill(HIST("PDGCode"), 19);
                    if(fabs(canpdg_1)>1000) histos.fill(HIST("PDGCode"), 20);
                    histos.fill(HIST("Chi2_False"), chi2_1);
                  }*/
                  for(int i=1; i<data.size(); i++){
                    if(data.at(i).at(2)==truepair.kid()){
                      if(i==1){
                        histos.fill(HIST("CosSimilarity_runk2"), data.at(1).at(1));
                      }else if(i==2){
                        histos.fill(HIST("CosSimilarity_runk3"), data.at(2).at(1));
                      }else if(i==3){
                        histos.fill(HIST("CosSimilarity_runk4"), data.at(3).at(1));
                      }else if(i==4){
                        histos.fill(HIST("CosSimilarity_runk5"), data.at(4).at(1));
                      }else if(i==5){
                        histos.fill(HIST("CosSimilarity_runk6"), data.at(5).at(1));
                      }else if(i==6){
                        histos.fill(HIST("CosSimilarity_runk7"), data.at(6).at(1));
                      }else if(i==7){
                        histos.fill(HIST("CosSimilarity_runk8"), data.at(7).at(1));
                      }else if(i==8){
                        histos.fill(HIST("CosSimilarity_runk9"), data.at(8).at(1));
                      }else if(i==9){
                        histos.fill(HIST("CosSimilarity_runk10"), data.at(9).at(1));
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<McInfomation>(cfgc),
    adaptAnalysisTask<McInfomation2>(cfgc),
    adaptAnalysisTask<MyAnalysisTask>(cfgc)
  };
}