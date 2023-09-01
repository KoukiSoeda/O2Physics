/// \author Koki Soeda
/// \since 24/07/2023

#include <iostream>
#include <array>
#include <cmath>
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
DECLARE_SOA_COLUMN(KdcaX, kdcax, float);
DECLARE_SOA_COLUMN(KdcaY, kdcay, float);
DECLARE_SOA_COLUMN(ColX, colx, float);
DECLARE_SOA_COLUMN(ColY, coly, float);
DECLARE_SOA_COLUMN(ColZ, colz, float);
DECLARE_SOA_COLUMN(SecVerX, secverx, float);
DECLARE_SOA_COLUMN(SecVerY, secvery, float);
DECLARE_SOA_COLUMN(SecVerZ, secverz, float);
DECLARE_SOA_COLUMN(PcaX, pcax, float);
DECLARE_SOA_COLUMN(PcaY, pcay, float);
DECLARE_SOA_COLUMN(PcaZ, pcaz, float);
DECLARE_SOA_COLUMN(R, r, float);
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
                  truepair::R);
}

struct McInfomation{ // PCA estimation via Single muon and Single Kaon (ture pair)
  Produces<aod::MCPair> mcpairtable;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsXY{"nBinsXY", 10010, "N bins in X,Y axis"};
  Configurable<int> nBinsZ{"nBinsZ", 110010, "N bins in Z axis"};
  Configurable<int> nBinsR{"nBinsR", 10010, "N bins in R"};
	Configurable<int> nBinsDist{"nBinsDist", 100010, ""};
	Configurable<int> fwdTrackType{"fwdTrackType", 0, "N TrackType in fwd"};
  Configurable<int> nBinsCut{"nBinsCut", 20010, "N bins in Cut"};

  void init(InitContext const&){
    //const AxisSpec axisVertex{nBinsVer, -30, +30, "cm"};
		const AxisSpec axistrackType{6, -0.5, 5.5, "TrackType"};
		const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisX{nBinsXY, -5.0005, 5.0005, "x(cm)"};
    const AxisSpec axisY{nBinsXY, -5.0005, 5.0005, "y(cm)"};
    const AxisSpec axisZ{nBinsZ, -60.0005, 50.0005, "z(cm)"};
    const AxisSpec axisCut{nBinsCut, -10.005, 10.005, "z(cm)"};
    const AxisSpec axisR{nBinsR, -5.005, 5.005, "r(cm)"};
		const AxisSpec axisDist{nBinsDist, -0.0005, 10.0005, "z(cm)"};
    const AxisSpec axisMFT{510010, -1.0005, 50.0005, "xy(cm)"};
    const AxisSpec axisDCAT{100010, -0.0005, 10.0005, "DCAT"};
    const AxisSpec axispT{10010, -0.005, 10.005, "pT"};
    const AxisSpec axisAccept{101, -0.05, 1.05, ""};

    histos.add("TrackType", "TrackType", kTH1F, {axistrackType});
    histos.add("Distance_D0_z", "Distance_D0_z", kTH1F, {axisDist});
    histos.add("PredictZ", "PredictZ", kTH1F, {axisZ});
    histos.add("PCAZ_Answer", "PCAZ_Answer", kTH1F, {axisZ});
    histos.add("PredictX", "PredictX", kTH1F, {axisZ});
    histos.add("PCAX_Answer", "PCAX_Answer", kTH1F, {axisZ});
    histos.add("PredictY", "PredictY", kTH1F, {axisZ});
    histos.add("PCAY_Answer", "PCAY_Answer", kTH1F, {axisZ});
    histos.add("PredictXY", "PredictXY", kTH1F, {axisR});
    histos.add("MCvertex_Estvertex", "MCvertex_Estvertex", kTH1F, {axisR});
    histos.add("PCAtoEst_R", "PCAtoEst_R", kTH1F, {axisR});
    histos.add("MFT_XY", "MFT_XY", kTH1F, {axisMFT});
    histos.add("MFT_XY_lowpt", "MFT_XY_lowpt", kTH1F, {axisMFT});
    histos.add("MFT_XY_midpt", "MFT_XY_midpt", kTH1F, {axisMFT});
    histos.add("MFT_XY_highpt", "MFT_XY_highpt", kTH1F, {axisMFT});
    histos.add("MuonpT_fromD0", "MuonpT_fromD0", kTH1F, {axispT});
    histos.add("DCAT_Kaon", "DCAT_Kaon", kTH1F, {axisDCAT});
    histos.add("Accept_rate_2sigma", "Accept_rate_2sigma", kTH1F, {axisAccept});
    histos.add("Primary_XY", "Primary_XY", kTH1F, {axisR});
    histos.add("r_muk_plus", "r_muk_plus", kTH1F, {axisR});
    histos.add("r_muk_minus", "r_muk_minus", kTH1F, {axisR});
    histos.add("nCandidate", "nCandidate", kTH1F, {{151, -0.5, 150.5, ""}});
		histos.add("counter1", "counter1", kTH1F, {axisCounter});
		histos.add("counter_1sigma", "counter_1sigma", kTH1F, {axisCounter});
    histos.add("counter_2sigma", "counter_2sigma", kTH1F, {axisCounter});
    histos.add("counter_3sigma", "counter_3sigma", kTH1F, {axisCounter});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
							 soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
							 soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
							 aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
							 aod::McCollisions const&)
  {
    if(collision.has_mcCollision()){
      vector<uint64_t> ambTrackIds;
      ambTrackIds.clear();
      for(auto& amfttrack : amfttracks){
        ambTrackIds.push_back(amfttrack.mfttrackId());
      }
      for(auto& fwdtrack : fwdtracks){
        if(!fwdtrack.has_collision() || !fwdtrack.has_mcParticle()) continue;
        histos.fill(HIST("TrackType"), fwdtrack.trackType());
        if(fwdtrack.trackType()!=0) continue;
        auto mcParticle_fwd = fwdtrack.mcParticle();

        if(fabs(mcParticle_fwd.pdgCode())==13 && mcParticle_fwd.has_mothers()){
          auto fwdColId = fwdtrack.collisionId();
          //auto fwdindex = mcParticle_fwd.globalIndex();
          auto mcMoms = mcParticle_fwd.mothers_as<aod::McParticles>();
          auto mcMom = mcMoms.back();
          auto mcMom_f = mcMoms.front();
          //LOGF(info, "GFT_PDG: %d, Mom_front: %d, Mom_back: %d", mcParticle_fwd.pdgCode(), mcMom.pdgCode(), mcMom_f.pdgCode());
          //auto momindex = mcMom.globalIndex();

          if(collision.globalIndex()==fwdColId){
            if(mcMom_f.pdgCode()==mcMom.pdgCode() && fabs(mcMom.pdgCode())==421){
              auto Daughters = mcMom.daughters_as<aod::McParticles>();
              int dc = 0;
              int fc = 0;
              int mcmucolid, mckcolid, mucolid, mcmuid, mupdg, mckid, kpdg;
              float mcmupt, Col_x, Col_y, Col_z, mudcaX, mudcaY, mumftX, mumftY, mumftZ, kdcaX, kdcaY, kmftX, kmftY, kmftZ, mumcX, mumcY, mumcZ, kmcZ;
              Col_x = collision.posX();
              Col_y = collision.posY();
              Col_z = collision.posZ();
              for(auto& Daughter : Daughters){   
                if(fabs(Daughter.pdgCode())==13){
                  auto muindex = Daughter.globalIndex();
                  for(auto& mfttrack : mfttracks){
                    if(mfttrack.has_collision() && mfttrack.has_mcParticle()){
                      if(mfttrack.collisionId()==fwdColId){
                        auto mcParticle_mft = mfttrack.mcParticle();
                        if(mcParticle_mft.globalIndex()==muindex){
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
                          mcmupt = mcParticle_mft.pt();
                          mucolid = mfttrack.collisionId();
                          mcmucolid = mcParticle_mft.mcCollisionId();
                          mcmuid = mfttrack.mcParticleId();
                          mupdg = mcParticle_mft.pdgCode();
                          dc++;
                        }
                      }
                    }
                  }
                }
                if(fabs(Daughter.pdgCode())==321){
                  auto kindex = Daughter.globalIndex();
                  for(auto& mfttrack2 : mfttracks){
                    if(mfttrack2.has_collision() && mfttrack2.has_mcParticle()){
                      if(mfttrack2.collisionId()==fwdColId){
                        auto mcParticle_mft2 = mfttrack2.mcParticle();
                        if(mcParticle_mft2.globalIndex()==kindex){
                          double mftchi2 = mfttrack2.chi2();
                          SMatrix5 mftpars(mfttrack2.x(), mfttrack2.y(), mfttrack2.phi(), mfttrack2.tgl(), mfttrack2.signed1Pt());
                          vector<double> mftv2;
                          SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
                          o2::track::TrackParCovFwd mftpars2{mfttrack2.z(), mftpars, mftcovs2, mftchi2};
                          mftpars2.propagateToZlinear(Col_z);
                          kdcaX = mftpars2.getX();
                          kdcaY = mftpars2.getY();
                          kmftX = mfttrack2.x();
                          kmftY = mfttrack2.y();
                          kmftZ = mfttrack2.z();
                          kmcZ = mcParticle_mft2.vz();
                          mckcolid = mcParticle_mft2.mcCollisionId();
                          mckid = mfttrack2.mcParticleId();
                          kpdg = mcParticle_mft2.pdgCode();
                          fc++;
                        }
                      }
                    }
                  }
                }
              }
              if(mucolid==fwdColId && collision.mcCollisionId()==mcmucolid){
                if(dc+fc==2 && mumcZ==kmcZ){
                  histos.fill(HIST("Distance_D0_z"), fabs(mcParticle_fwd.vz()-mcMom.vz()));
                  histos.fill(HIST("DCAT_Kaon"), sqrt(pow(kdcaX,2)+pow(kdcaY,2)));
                  auto Nax = mudcaX - mumftX;
                  auto Nay = mudcaY - mumftY;
                  auto Naz = Col_z - mumftZ;
                  auto Ncx = kdcaX - kmftX;
                  auto Ncy = kdcaY - kmftY;
                  auto Ncz = Col_z - kmftZ;
                  auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
                  auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
                  auto A3 = (mumftX-kmftX)*Nax + (mumftY-kmftY)*Nay + (mumftZ-kmftZ)*Naz;
                  auto B1 = A2;
                  auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
                  auto B3 = (kmftX-mumftX)*Ncx + (kmftY-mumftY)*Ncy + (kmftZ-mumftZ)*Ncz;
                  auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
                  auto s = -((A2*t+A3)/A1);
                  float pre_mux = mumftX + s*Nax;
                  float pre_muy = mumftY + s*Nay;
                  float pre_muz = mumftZ + s*Naz;
                  float pre_kx = kmftX + t*Ncx;
                  float pre_ky = kmftY + t*Ncy;
                  float pre_kz = kmftZ + t*Ncz;
                  float pcaX = (pre_mux+pre_kx)/2;
                  float pcaY = (pre_muy+pre_ky)/2;
                  float pcaZ = (pre_muz+pre_kz)/2;

                  float s_org = (mumcZ-mumftZ)/Naz;
                  float t_org = (mumcZ-kmftZ)/Ncz;

                  float r = sqrt(pow((mumftX+s_org*Nax)-(kmftX+t_org*Ncx), 2)+pow((mumftY+s_org*Nay)-(kmftY+t_org*Ncy), 2)+pow((mumftZ+s_org*Naz)-(kmftZ+t_org*Ncz), 2));
                  float r_org = sqrt(pow(pcaX-mumcX, 2)+pow(pcaY-mumcY, 2)+pow(pcaZ-mumcZ, 2));
                  histos.fill(HIST("PredictZ"), pcaZ-Col_z);
                  histos.fill(HIST("PCAtoEst_R"), r);
                  histos.fill(HIST("PCAZ_Answer"), mumcZ-collision.mcCollision().posZ());
                  histos.fill(HIST("PredictX"), pcaX-Col_x);
                  histos.fill(HIST("PCAX_Answer"), mumcX-Col_x);
                  histos.fill(HIST("PredictY"), pcaY-Col_y);
                  histos.fill(HIST("PCAY_Answer"), mumcY-Col_y);
                  float prex = pcaX-Col_x;
                  float prey = pcaY-Col_y;
                  float preR = sqrt(pow(prex,2)+pow(prey,2));
                  histos.fill(HIST("PredictXY"), preR);
                  histos.fill(HIST("MCvertex_Estvertex"), r_org);
                  float MFT_r = sqrt(pow(mumftX-kmftX,2)+pow(mumftY-kmftY,2));
                  histos.fill(HIST("MFT_XY"), MFT_r);
                  histos.fill(HIST("MuonpT_fromD0"), mcmupt);
                  if(mcmupt>0.5 && mcmupt<=1.0) histos.fill(HIST("MFT_XY_lowpt"), MFT_r);
                  if(mcmupt>1.0 && mcmupt<=2.0) histos.fill(HIST("MFT_XY_midpt"), MFT_r);
                  if(mcmupt>2.0 && mcmupt<=4.0) histos.fill(HIST("MFT_XY_highpt"), MFT_r);
                  mcpairtable(mucolid,mcmuid,mupdg,mckid,kpdg,kdcaX,kdcaY,Col_x,Col_y,Col_z,mumcX,mumcY,mumcZ,pcaX-Col_x,pcaY-Col_y,pcaZ-Col_z,r);
                  if(pcaZ>=0){
                    auto r_muk_plus = sqrt(pow((mudcaX-kdcaX)-Col_x,2)+pow((mudcaY-kdcaY)-Col_y,2));
                    histos.fill(HIST("r_muk_plus"), r_muk_plus);
                  }else{
                    auto r_muk_minus = sqrt(pow(pre_mux-pre_kx,2)+pow(pre_muy-pre_ky,2)+pow(pre_muz-pre_kz,2));
                    histos.fill(HIST("r_muk_minus"), r_muk_minus);
                  }
                  /*int ntrack_perCol_int = 0;
                  int ntrack_accept_int = 0;
                  int nCan = 0;
                  for(auto& otmfttrack : mfttracks){
                    if(otmfttrack.has_collision() && otmfttrack.has_mcParticle()){
                      if(otmfttrack.collisionId()==fwdColId){

                        if(otmfttrack.mcParticleId()!=mcmuid){
                          if(find(ambTrackIds.begin(), ambTrackIds.end(), otmfttrack.globalIndex()) != ambTrackIds.end()) continue;
                          if(otmfttrack.sign()!=0) nCan++;
                        }

                        if(otmfttrack.mcParticleId()!=mcmuid && otmfttrack.mcParticleId()!=mckid){
                          if(find(ambTrackIds.begin(), ambTrackIds.end(), otmfttrack.globalIndex()) != ambTrackIds.end()) continue;
                          double mftchi2 = otmfttrack.chi2();
                          SMatrix5 mftpars(otmfttrack.x(), otmfttrack.y(), otmfttrack.phi(), otmfttrack.tgl(), otmfttrack.signed1Pt());
                          vector<double> mftv1;
                          SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                          o2::track::TrackParCovFwd mftpars1{otmfttrack.z(), mftpars, mftcovs, mftchi2};
                          auto z_mean = -0.2236;
                          auto z_std = 0.8555;
                          auto xy_mean = 0.06243;
                          auto xy_std = 0.08104;
                          int count_r_1sigma = 0;
                          int count_r_2sigma = 0;
                          int count_r_3sigma = 0;
                          ntrack_perCol_int++;
                          if(otmfttrack.mcParticle().producedByGenerator()==true){
                            mftpars1.propagateToZlinear(Col_z);
                            histos.fill(HIST("Primary_XY"), sqrt(pow(mftpars1.getX(),2)+pow(mftpars1.getY(),2)));
                          }
                          for(int i=0; i<=1000000; i++){
                            auto cut_z_1sigma = ((z_mean-z_std) + 0.0001*i) + Col_z;
                            auto cut_z_2sigma = ((z_mean-2*z_std) + 0.0001*i) + Col_z;
                            auto cut_z_3sigma = ((z_mean-3*z_std) + 0.0001*i) + Col_z;
                            if(cut_z_3sigma>=(z_mean+3*z_std)+Col_z) break;
                            mftpars1.propagateToZlinear(cut_z_3sigma);
                            auto cut_xaxis = mftpars1.getX() - Col_x;
                            auto cut_yaxis = mftpars1.getY() - Col_y;
                            auto cut_r = sqrt(pow(cut_xaxis,2)+pow(cut_yaxis,2));
                            if(cut_r<=(xy_mean+3*xy_std) && count_r_3sigma==0){
                              histos.fill(HIST("counter_3sigma"), 0.5);
                              count_r_3sigma++;
                            }

                            if(cut_z_2sigma<=(z_mean+2*z_std)+Col_z){
                              mftpars1.propagateToZlinear(cut_z_2sigma);
                              auto cut_xaxis = mftpars1.getX() - Col_x;
                              auto cut_yaxis = mftpars1.getY() - Col_y;
                              auto cut_r = sqrt(pow(cut_xaxis,2)+pow(cut_yaxis,2));
                              if(cut_r<=(xy_mean+2*xy_std) && count_r_2sigma==0){
                                histos.fill(HIST("counter_2sigma"), 0.5);
                                count_r_2sigma++;
                                ntrack_accept_int++;
                              }
                            }
                            
                            if(cut_z_1sigma<=(z_mean+z_std)+Col_z){
                              mftpars1.propagateToZlinear(cut_z_1sigma);
                              auto cut_xaxis = mftpars1.getX() - Col_x;
                              auto cut_yaxis = mftpars1.getY() - Col_y;
                              auto cut_r = sqrt(pow(cut_xaxis,2)+pow(cut_yaxis,2));
                              if(cut_r<=(xy_mean+xy_std) && count_r_1sigma==0){
                                histos.fill(HIST("counter_1sigma"), 0.5);
                                count_r_1sigma++;                            
                              }
                            }
                            //if(cut_z==0)LOGF(info, "R: %f, X: %f, Y: %f", cut_r,cut_xaxis,cut_yaxis);
                          }
                          histos.fill(HIST("counter1"), 0.5);
                        }
                      }
                    }
                  }
                  if(nCan!=0) histos.fill(HIST("nCandidate"), nCan);
                  float ntrack_perCol_float = ntrack_perCol_int;
                  int ntrack_accept_float = ntrack_accept_int;
                  float rate_accept = ntrack_accept_float/ntrack_perCol_float;
                  histos.fill(HIST("Accept_rate_2sigma"), rate_accept);*/
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
  Configurable<int> nBinsPt{"nBinsPt", 25010, "N bins pt axis"};
  Configurable<int> nBinsZ{"nBinsZ", 1100010, "N bins in Z axis"};
  Configurable<int> nBinsR{"nBinsR", 60010, "N bins in R"};
  Configurable<int> nBinsCut{"nBinsCut", 2001, "N bins in cut"};
  Configurable<int> nBinsDCAT{"nBinsDCAT", 50001, "N bins in DCAT"};

  void init(InitContext const&){
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisEvent{1000000, 0, 1000000};
    const AxisSpec axisTest{1, 0, 1, ""};
    const AxisSpec axisPt{nBinsPt, -0.005, 25.005, "p_{T}"};
    const AxisSpec axisR(nBinsR, -1.0005, 5.0005, "r(cm)");
    const AxisSpec axisZ(nBinsZ, -60.0005, 50.0005, "z(cm)");
    const AxisSpec axisCutZ(nBinsCut, -20.005, 0.005, "z(cm)");
    const AxisSpec axisParticle(101, -0.5, 100.5, "N");
    const AxisSpec axisDCAT(nBinsDCAT, -0.0005, 5.0005, "DCAT");
    const AxisSpec axisPDG{21, -0.5, 20.5, ""};
    const AxisSpec axisMFT{510010, -1.0005, 50.0005, "xy(cm)"};

    histos.add("Est_PCAZ", "Est_PCAZ", kTH1F, {axisZ});
    histos.add("Est_PCAZ_Correct", "Est_PCAZ_Correct", kTH1F, {axisZ});
    histos.add("Est_PCAZ_False", "Est_PCAZ_False", kTH1F, {axisZ});
    histos.add("CutZ_Correct", "CutZ_Correct", kTH1F, {axisCutZ});
    histos.add("EstZ_Correct", "EstZ_Correct", kTH1F, {axisZ});
    histos.add("MFT_XY_True", "MFT_XY_True", kTH1F, {axisMFT});
    histos.add("MFT_XY_Fake", "MFT_XY_Fake", kTH1F, {axisMFT});
    histos.add("MFT_XY_Dist", "MFT_XY_Dist", kTH1F, {axisMFT});
    histos.add("pairPDG", "pairPDG", kTH1F, {axisPDG});
    histos.add("Correct_pt", "Correct_pt", kTH1F, {axisPt});
    histos.add("Correct_r_dist", "Correct_r_dist", kTH1F, {axisR});
    histos.add("Incorrect_r_dist", "Incorrect_r_dist", kTH1F, {axisR});
    histos.add("Actual_r", "Actual_r", kTH1F, {axisR});
    histos.add("PosDiff_x", "PosDiff_x", kTH1F, {axisZ});
    histos.add("PosDiff_y", "PosDiff_y", kTH1F, {axisZ});
    histos.add("PosDiff_z", "PosDiff_z", kTH1F, {axisZ});
    histos.add("InnerParticles", "InnerParticles", kTH1F, {axisParticle});
    histos.add("CorDCAXY", "CorDCAXY", kTH1F, {axisDCAT});
    histos.add("KaonDCAXY", "KaonDCAXY", kTH1F, {axisDCAT});
    histos.add("SelDCAXY", "SelDCAXY", kTH1F, {axisDCAT});
    histos.add("Total_pt", "Total_pt", kTH1F, {axisPt});

    auto mumompdg_d = histos.get<TH1>(HIST("pairPDG"));
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
        auto mcMoms = mcParticle_fwd.mothers_as<aod::McParticles>();
        auto mcMom = mcMoms.back();
        auto mcMom_f = mcMoms.front();
        auto Daughters = mcMom.daughters_as<aod::McParticles>();

        if(fabs(mcMom_f.pdgCode())!=421) continue;

        int mucolid, mcmuid, mckid_org, kpdg, estmumomid, estkmomid, mumom, kmom_org, storemuid, estmck, estpdg;
        float mupt, Col_x, Col_y, Col_z, mudcaX, mudcaY, mumftX, mumftY, mumftZ, kdcaX, kdcaY, kmftX, kmftY, kmftZ, mumcX, mumcY, mumcZ, kmcZ, estx, esty, estz;
        float estdcax, estdcay;
        
        for(auto& mfttrack : mfttracks){
          if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
          //muon like
          if(mfttrack.collisionId()==fwdColId && mfttrack.mcParticleId()==fwdtrack.mcParticleId()){
            Col_x = collision.posX();
            Col_y = collision.posY();
            Col_z = collision.posZ();
            auto mcParticle_mft = mfttrack.mcParticle();
            double mftchi2 = mfttrack.chi2();
            SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
            vector<double> mftv1;
            SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
            o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
            mftpars1.propagateToZlinear(Col_z);
            mupt = fwdtrack.pt();
            mudcaX = mftpars1.getX();
            mudcaY = mftpars1.getY();
            mumftX = mfttrack.x();
            mumftY = mfttrack.y();
            mumftZ = mfttrack.z();
            mumcX = mcParticle_mft.vx();
            mumcY = mcParticle_mft.vy();
            mumcZ = mcParticle_mft.vz();
            mupt = mcParticle_mft.pt();
            mucolid = mfttrack.collisionId();
            mcmuid = mfttrack.mcParticleId();
            if(!mcParticle_mft.has_mothers()) continue;
            auto muMoms = mcParticle_mft.mothers_as<aod::McParticles>();
            mumom = muMoms.back().globalIndex();
            int npar = 0;

            for(auto& truepair : truepairs){
              if(truepair.muid()==mcmuid && truepair.colid()==mucolid){
                if(storemuid==mcmuid) LOGF(info, "Analysis MuID: %d", mcmuid);
                //for(int cut=1; cut<=2000; cut++){
                  //float fcut = -5;//-0.01*cut;
                  int colID_ref;
                  float mft_xy;
                  float dist_2track = 10000;
                  for(auto& mfttrack2 : mfttracks){
                    if(mfttrack2.has_collision() && mfttrack2.has_mcParticle()){
                      if(mfttrack2.collisionId()==fwdColId && mfttrack2.mcParticleId()!=mcmuid){
                        if(find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack2.globalIndex()) != ambTrackIds.end()) continue;
                        auto mcParticle_mft2 = mfttrack2.mcParticle();
                        double mftchi2 = mfttrack2.chi2();
                        SMatrix5 mftpars(mfttrack2.x(), mfttrack2.y(), mfttrack2.phi(), mfttrack2.tgl(), mfttrack2.signed1Pt());
                        vector<double> mftv2;
                        SMatrix55 mftcovs2(mftv2.begin(), mftv2.end());
                        o2::track::TrackParCovFwd mftpars2{mfttrack2.z(), mftpars, mftcovs2, mftchi2};
                        mftpars2.propagateToZlinear(Col_z);
                        kdcaX = mftpars2.getX();
                        kdcaY = mftpars2.getY();
                        kmftX = mfttrack2.x();
                        kmftY = mfttrack2.y();
                        kmftZ = mfttrack2.z();
                        kmcZ = mcParticle_mft2.vz();
                        mckid_org = mfttrack2.mcParticleId();
                        kpdg = mcParticle_mft2.pdgCode();
                        if(mcParticle_mft2.has_mothers()){
                          auto kMoms = mcParticle_mft2.mothers_as<aod::McParticles>();
                          kmom_org = kMoms.back().globalIndex();
                        }else{
                          kmom_org = -1;
                        }
                        
                        if(truepair.kid()!=mckid_org) histos.fill(HIST("MFT_XY_Dist"), sqrt(pow(mumftX-kmftX,2)+pow(mumftY-kmftY,2)));
                        auto Nax = mudcaX - mumftX;
                        auto Nay = mudcaY - mumftY;
                        auto Naz = Col_z - mumftZ;
                        auto Ncx = kdcaX - kmftX;
                        auto Ncy = kdcaY - kmftY;
                        auto Ncz = Col_z - kmftZ;
                        auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
                        auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
                        auto A3 = (mumftX-kmftX)*Nax + (mumftY-kmftY)*Nay + (mumftZ-kmftZ)*Naz;
                        auto B1 = A2;
                        auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
                        auto B3 = (kmftX-mumftX)*Ncx + (kmftY-mumftY)*Ncy + (kmftZ-mumftZ)*Ncz;
                        auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
                        auto s = -((A2*t+A3)/A1);
                        float pre_mux = mumftX + s*Nax;
                        float pre_muy = mumftY + s*Nay;
                        float pre_muz = mumftZ + s*Naz;
                        float pre_kx = kmftX + t*Ncx;
                        float pre_ky = kmftY + t*Ncy;
                        float pre_kz = kmftZ + t*Ncz;
                        //float k_unit = sqrt(pow(Ncx,2)+pow(Ncy,2)+pow(Ncz,2));
                        float pcax_org = ((pre_mux+pre_kx)/2) - Col_x;
                        float pcay_org = ((pre_muy+pre_ky)/2) - Col_y;
                        float pcaz_org = ((pre_muz+pre_kz)/2) - Col_z;

                        float pcax, pcay, pcaz, r_xyz;
                        int kmom, mckid;

                        float s_org = (truepair.secverz()-mumftZ)/Naz;
                        float t_org = (truepair.secverz()-kmftZ)/Ncz;
                        float actr_r = sqrt(pow((mumftX+s_org*Nax)-(kmftX+t_org*Ncx), 2)+pow((mumftY+s_org*Nay)-(kmftY+t_org*Ncy), 2)+pow((mumftZ+s_org*Naz)-(kmftZ+t_org*Ncz), 2));
                        if(actr_r<=truepair.r()){
                          npar++;
                        }

                        //if(pcaz_org<0 && pcaz_org>=fcut){
                          r_xyz = sqrt(pow(pre_kx-pre_mux, 2)+pow(pre_ky-pre_muy, 2)+pow(pre_kz-pre_muz, 2));
                          pcax = pcax_org;
                          pcay = pcay_org;
                          pcaz = pcaz_org;
                          kmom = kmom_org;
                          mckid = mckid_org;

                          if(dist_2track>=r_xyz){
                            dist_2track = r_xyz;
                            estx = pcax;
                            esty = pcay;
                            estz = pcaz;
                            estmumomid = mumom;
                            estkmomid = kmom;
                            estmck = mckid;
                            estdcax = kdcaX;
                            estdcay = kdcaY;
                            estpdg = kpdg;
                            mft_xy = sqrt(pow(mumftX-kmftX,2)+pow(mumftY-kmftY,2));
                            colID_ref = mfttrack2.collisionId();
                          }
                        //}
                      }
                    }
                  }

                  //if(estmumomid==estkmomid && truepair.kid()==estmck) histos.fill(HIST("CutZ_Correct"), fcut);

                  //if(fcut==(-5)){
                    histos.fill(HIST("Est_PCAZ"), estz);
                    histos.fill(HIST("Total_pt"), mupt);
                    
                    if(fabs(estpdg)<38) histos.fill(HIST("pairPDG"), 0);
                    if(fabs(estpdg)==111) histos.fill(HIST("pairPDG"), 1);
                    if(fabs(estpdg)==211) histos.fill(HIST("pairPDG"), 2);
                    if(fabs(estpdg)==113) histos.fill(HIST("pairPDG"), 3);
                    if(fabs(estpdg)==130) histos.fill(HIST("pairPDG"), 4);
                    if(fabs(estpdg)==221) histos.fill(HIST("pairPDG"), 5);
                    if(fabs(estpdg)==223) histos.fill(HIST("pairPDG"), 6);
                    if(fabs(estpdg)==310) histos.fill(HIST("pairPDG"), 7);
                    if(fabs(estpdg)==313) histos.fill(HIST("pairPDG"), 8);
                    if(fabs(estpdg)==321) histos.fill(HIST("pairPDG"), 9);
                    if(fabs(estpdg)==323) histos.fill(HIST("pairPDG"), 10);
                    if(fabs(estpdg)==331) histos.fill(HIST("pairPDG"), 11);
                    if(fabs(estpdg)==333) histos.fill(HIST("pairPDG"), 12);
                    if(fabs(estpdg)==411) histos.fill(HIST("pairPDG"), 13);
                    if(fabs(estpdg)==413) histos.fill(HIST("pairPDG"), 14);
                    if(fabs(estpdg)==421) histos.fill(HIST("pairPDG"), 15);
                    if(fabs(estpdg)==431) histos.fill(HIST("pairPDG"), 16);
                    if(fabs(estpdg)==511) histos.fill(HIST("pairPDG"), 17);
                    if(fabs(estpdg)==521) histos.fill(HIST("pairPDG"), 18);
                    if(fabs(estpdg)==531) histos.fill(HIST("pairPDG"), 19);
                    if(fabs(estpdg)>1000) histos.fill(HIST("pairPDG"), 20);
                    //histos.fill(HIST("InnerParticles"), npar);

                    if(estmumomid==estkmomid && truepair.kid()==estmck){
                      histos.fill(HIST("MFT_XY_True"), mft_xy);
                      histos.fill(HIST("Correct_pt"), mupt);
                      histos.fill(HIST("Correct_r_dist"), dist_2track);
                      histos.fill(HIST("Est_PCAZ_Correct"), estz);
                      histos.fill(HIST("EstZ_Correct"), estz-(truepair.secverz()-Col_z));
                      histos.fill(HIST("InnerParticles"), npar);
                      //histos.fill(HIST("CutZ_Correct"), fcut);
                      if(estdcax==truepair.kdcax() && estdcay==truepair.kdcay()){
                        float CorDCAT = sqrt(pow(estdcax-Col_x,2)+pow(estdcay-Col_y,2));
                        histos.fill(HIST("CorDCAXY"), CorDCAT);
                      }
                    }else{
                      histos.fill(HIST("Est_PCAZ_False"), estz);
                      histos.fill(HIST("MFT_XY_Fake"), mft_xy);
                      histos.fill(HIST("Actual_r"), truepair.r());
                      histos.fill(HIST("PosDiff_x"), estx-truepair.pcax());
                      histos.fill(HIST("PosDiff_y"), esty-truepair.pcay());
                      histos.fill(HIST("PosDiff_z"), estz-truepair.pcaz());
                      float trueDCAT = sqrt(pow(truepair.kdcax()-Col_x,2)+pow(truepair.kdcay()-Col_y,2));
                      float estDCAT = sqrt(pow(estdcax-Col_x,2)+pow(estdcay-Col_y,2));
                      histos.fill(HIST("KaonDCAXY"), trueDCAT);
                      histos.fill(HIST("SelDCAXY"), estDCAT);
                      histos.fill(HIST("Incorrect_r_dist"), dist_2track);

                    }
                  //}
                  storemuid = mcmuid;
                //}
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
    adaptAnalysisTask<MyAnalysisTask>(cfgc)};
}