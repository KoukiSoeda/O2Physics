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
  Configurable<int> nBinsZ{"nBinsZ", 1100010, "N bins in Z axis"};
  Configurable<int> nBinsR{"nBinsR", 110010, "N bins in R"};
	Configurable<int> nBinsDist{"nBinsDist", 100010, ""};
	Configurable<int> fwdTrackType{"fwdTrackType", 0, "N TrackType in fwd"};

  void init(InitContext const&){
    //const AxisSpec axisVertex{nBinsVer, -30, +30, "cm"};
		const AxisSpec axistrackType{6, -0.5, 5.5, "TrackType"};
		const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisZ{nBinsZ, -60.0005, 50.0005, "z(cm)"};
    const AxisSpec axisR{nBinsR, -1.0005, 10.0005, "r(cm)"};
		const AxisSpec axisDist{nBinsDist, -0.0005, 10.0005, "z(cm)"};
    const AxisSpec axisMFT{510010, -1.0005, 50.0005, "xy(cm)"};
    const AxisSpec axisDCAT{100010, -0.0005, 10.0005, "DCAT"};
    const AxisSpec axispT{10010, -0.005, 10.005, "pT"};

    histos.add("TrackType", "TrackType", kTH1F, {axistrackType});
    histos.add("Distance_D0_z", "Distance_D0_z", kTH1F, {axisDist});
    histos.add("PredictZ", "PredictZ", kTH1F, {axisZ});
    histos.add("PCAZ_Answer", "PCAZ_Answer", kTH1F, {axisZ});
    histos.add("PredictX", "PredictX", kTH1F, {axisZ});
    histos.add("PCAX_Answer", "PCAX_Answer", kTH1F, {axisZ});
    histos.add("PredictY", "PredictY", kTH1F, {axisZ});
    histos.add("PCAY_Answer", "PCAY_Answer", kTH1F, {axisZ});
    histos.add("MCvertex_Estvertex", "MCvertex_Estvertex", kTH1F, {axisR});
    histos.add("PCAtoEst_R", "PCAtoEst_R", kTH1F, {axisR});
    histos.add("MFT_XY", "MFT_XY", kTH1F, {axisMFT});
    histos.add("MFT_XY_lowpt", "MFT_XY_lowpt", kTH1F, {axisMFT});
    histos.add("MFT_XY_midpt", "MFT_XY_midpt", kTH1F, {axisMFT});
    histos.add("MFT_XY_highpt", "MFT_XY_highpt", kTH1F, {axisMFT});
    histos.add("MuonpT_fromD0", "MuonpT_fromD0", kTH1F, {axispT});
    histos.add("DCAT_Kaon", "DCAT_Kaon", kTH1F, {axisDCAT});
		histos.add("counter1", "counter1", kTH1F, {axisCounter});
		histos.add("counter2", "counter2", kTH1F, {axisCounter});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
							 soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
							 soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
							 aod::McParticles const& MCparticle,
							 aod::McCollisions const&)
  {
    if(collision.has_mcCollision()){
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
                  histos.fill(HIST("MCvertex_Estvertex"), r_org);
                  float MFT_r = sqrt(pow(mumftX-kmftX,2)+pow(mumftY-kmftY,2));
                  histos.fill(HIST("MFT_XY"), MFT_r);
                  histos.fill(HIST("MuonpT_fromD0"), mcmupt);
                  if(mcmupt>0.5 && mcmupt<=1.0) histos.fill(HIST("MFT_XY_lowpt"), MFT_r);
                  if(mcmupt>1.0 && mcmupt<=2.0) histos.fill(HIST("MFT_XY_midpt"), MFT_r);
                  if(mcmupt>2.0 && mcmupt<=4.0) histos.fill(HIST("MFT_XY_highpt"), MFT_r);
                  mcpairtable(mucolid,mcmuid,mupdg,mckid,kpdg,kdcaX,kdcaY,Col_x,Col_y,Col_z,mumcX,mumcY,mumcZ,pcaX-Col_x,pcaY-Col_y,pcaZ-Col_z,r);
                }
              }
            }
          }
        }
      }
    }
  }
};

struct Test{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&){
    const AxisSpec axisTest{1, 12.5, 13.5, ""};
    histos.add("test", "test", kTH1F, {axisTest});
  }
  void process(aod::MCPair const& truepairs){
    for(auto& truepair : truepairs){
      //LOGF(info, "MuID: %d", truepair.muid());
      auto test = truepair.mupdg();
      histos.fill(HIST("test"), fabs(test));
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

                  //Topological vertex finder (v0813)
                  float vc_sum = 0.0;
                  float vc2_sum = 0.0;
                  float ave_x = 0.0;
                  float ave_y = 0.0;
                  float ave_z = 0.0;
                  int ntrack_int = 0;
                  float ntrack_float = 0.0;
                  float org_x = 0.0;
                  float org_y = 0.0;
                  float org_z = 0.0;

                  for(auto& mfttrack_ave : mfttracks){
                    if(mfttrack_ave.has_collision() && mfttrack_ave.has_mcParticle()){
                      if(mfttrack_ave.collisionId()==fwdColId){
                        if(find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack_ave.globalIndex()) != ambTrackIds.end()) continue;
                        double mftchi2 = mfttrack_ave.chi2();
                        SMatrix5 mftpars(mfttrack_ave.x(), mfttrack_ave.y(), mfttrack_ave.phi(), mfttrack_ave.tgl(), mfttrack_ave.signed1Pt());
                        vector<double> mftv_ave;
                        SMatrix55 mftcovs(mftv_ave.begin(), mftv_ave.end());
                        o2::track::TrackParCovFwd mftpars_ave{mfttrack_ave.z(), mftpars, mftcovs, mftchi2};
                        mftpars_ave.propagateToZlinear(Col_z);

                        auto a = mfttrack_ave.x() - mftpars_ave.getX();
                        auto b = mfttrack_ave.y() - mftpars_ave.getY();
                        auto c = mfttrack_ave.z() - Col_z;
                        auto t = (a*estx + b*esty + c*estz)/(a*a + b*b + c*c);

                        org_x += t*a - Col_x;
                        org_y += t*b - Col_y;
                        org_z += t*c - Col_z;
                        ntrack_int++;
                      }
                    }
                  }

                  if(ntrack_int==0) continue;
                  ntrack_float = ntrack_int;
                  ave_x = org_x/ntrack_float;
                  ave_y = org_y/ntrack_float;
                  ave_z = org_z/ntrack_float;
                  /*
                  double a11_0 = 0;
                  double a22_0 = 0;
                  double a33_0 = 0;
                  double a12_0 = 0;
                  double a13_0 = 0;
                  double a21_0 = 0;
                  double a23_0 = 0;
                  double a31_0 = 0;
                  double a32_0 = 0;

                  for(auto& mfttrack3 : mfttracks){
                    if(mfttrack3.has_collision() && mfttrack3.has_mcParticle()){
                      if(mfttrack3.collisionId()==fwdColId && mfttrack3.mcParticleId()!=mcmuid){
                        if(find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack3.globalIndex()) != ambTrackIds.end()) continue;
                        double mftchi2 = mfttrack3.chi2();
                        SMatrix5 mftpars(mfttrack3.x(), mfttrack3.y(), mfttrack3.phi(), mfttrack3.tgl(), mfttrack3.signed1Pt());
                        vector<double> mftv3;
                        SMatrix55 mftcovs(mftv3.begin(), mftv3.end());
                        o2::track::TrackParCovFwd mftpars3{mfttrack3.z(), mftpars, mftcovs, mftchi2};
                        mftpars3.propagateToZlinear(Col_z);

                        auto a0 = mfttrack3.x() - mftpars3.getX();
                        auto b0 = mfttrack3.y() - mftpars3.getY();
                        auto c0 = mfttrack3.z() - Col_z;
                        auto t = (a0*estx + b0*esty + c0*estz)/(a0*a0 + b0*b0 + c0*c0);
                        //最近接点
                        auto closest_x = t*a0 - Col_x;
                        auto closest_y = t*b0 - Col_y;
                        auto closest_z = t*c0 - Col_z;
                        //偏差ベクトル
                        double a = closest_x - ave_x;
                        double b = closest_y - ave_y;
                        double c = closest_z - ave_z;

                        //分散
                         a11_0 += a*a;
                         a22_0 += b*b;
                         a33_0 += c*c;
                        //共分散
                         a12_0 += a*b;
                         a13_0 += a*c;
                         a21_0 += b*a;
                         a23_0 += b*c;
                         a31_0 += c*a;
                         a32_0 += c*b;
                        //LOGF(info,"%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", a11,a12,a13,a21,a22,a23,a31,a32,a33);
                      }
                    }
                  }
                  double a11 = a11_0/ntrack_float;
                  double a12 = (a12_0/ntrack_float);//((sqrt(a11_0)*sqrt(a22_0)));
                  double a13 = (a13_0/ntrack_float);//((sqrt(a11_0)*sqrt(a33_0)));
                  double a21 = (a21_0/ntrack_float);//((sqrt(a11_0)*sqrt(a22_0)));
                  double a22 = a22_0/ntrack_float;
                  double a23 = (a23_0/ntrack_float);//((sqrt(a22_0)*sqrt(a33_0)));
                  double a31 = (a31_0/ntrack_float);//((sqrt(a11_0)*sqrt(a33_0)));
                  double a32 = (a32_0/ntrack_float);//((sqrt(a22_0)*sqrt(a33_0)));
                  double a33 = a33_0/ntrack_float;

                  double datv = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a31*a22*a13 - a21*a12*a33 - a32*a23*a11;
                  if(datv!=0) LOGF(info, "datv is not 0!: %f", datv);
                  if(datv==0) LOGF(info, "datv is 0");
                  
                  double b11 = a22*a33 - a32*a23;
                  double b12 = -(a12*a33 - a13*a32);
                  double b13 = a12*a23 - a22*a13;
                  double b21 = -(a21*a33 - a31*a23);
                  double b22 = a11*a33 - a31*a13;
                  double b23 = -(a11*a23 - a21*a13);
                  double b31 = a21*a32 - a31*a22;
                  double b32 = -(a11*a32 - a31*a12);
                  double b33 = a11*a22 - a21*a12;
                  //LOGF(info,"%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", b11/datv,b12/datv,b13/datv,b21/datv,b22/datv,b23/datv,b31/datv,b32/datv,b33/datv);
                  */

                  for(auto& mfttrack4 : mfttracks){
                    if(mfttrack4.has_collision() && mfttrack4.has_mcParticle()){
                      if(mfttrack4.collisionId()==fwdColId && mfttrack4.mcParticleId()!=mcmuid){
                        if(find(ambTrackIds.begin(), ambTrackIds.end(), mfttrack4.globalIndex()) != ambTrackIds.end()) continue;
                        double mftchi2 = mfttrack4.chi2();
                        SMatrix5 mftpars(mfttrack4.x(), mfttrack4.y(), mfttrack4.phi(), mfttrack4.tgl(), mfttrack4.signed1Pt());
                        vector<double> mftv4;
                        SMatrix55 mftcovs(mftv4.begin(), mftv4.end());
                        o2::track::TrackParCovFwd mftpars4{mfttrack4.z(), mftpars, mftcovs, mftchi2};
                        mftpars4.propagateToZlinear(Col_z);
                        auto a1 = mfttrack4.x() - mftpars4.getX();
                        auto b1 = mfttrack4.y() - mftpars4.getY();
                        auto c1 = mfttrack4.z() - Col_z;
                        auto t = (a1*estx + b1*esty + c1*estz)/(a1*a1 + b1*b1 + c1*c1);
                        //最近接点
                        auto closest_x = t*a1 - Col_x;
                        auto closest_y = t*b1 - Col_y;
                        auto closest_z = t*c1 - Col_z;
                        double l_x = estx - closest_x;
                        double l_y = esty - closest_y;
                        double l_z = estz - closest_z;
                        /*ave_x = (estx+closest_x)/2;
                        ave_y = (esty+closest_y)/2;
                        ave_z = (estz+closest_z)/2;*/
                        double a0 = estx-ave_x;
                        double a = closest_x-ave_x;
                        double b0 = esty-ave_y;
                        double b = closest_y-ave_y;
                        double c0 = estz-ave_z;
                        double c = closest_z-ave_z;
                        double a11 = (pow(a0,2)+pow(a,2))/2;
                        double a22 = (pow(b0,2)+pow(b,2))/2;
                        double a33 = (pow(c0,2)+pow(c,2))/2;
                        double a12 = (a0*b0+a*b)/2;//(2*sqrt(a11)*sqrt(a22));
                        double a13 = (a0*c0+a*c)/2;//(2*sqrt(a11)*sqrt(a33));
                        double a21 = (b0*a0+b*a)/2;//(2*sqrt(a22)*sqrt(a11));
                        double a23 = (b0*c0+b*c)/2;//(2*sqrt(a22)*sqrt(a33));
                        double a31 = (c0*a0+c*a)/2;//(2*sqrt(a33)*sqrt(a11));
                        double a32 = (c0*b0+c*b)/2;//(2*sqrt(a33)*sqrt(a22));
                        
                        double datv = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a31*a22*a13 - a21*a12*a33 - a32*a23*a11;
                        double b11 = (a22*a33 - a32*a23)/datv;
                        double b12 = -(a12*a33 - a13*a32)/datv;
                        double b13 = (a12*a23 - a22*a13)/datv;
                        double b21 = -(a21*a33 - a31*a23)/datv;
                        double b22 = (a11*a33 - a31*a13)/datv;
                        double b23 = -(a11*a23 - a21*a13)/datv;
                        double b31 = (a21*a32 - a31*a22)/datv;
                        double b32 = -(a11*a32 - a31*a12)/datv;
                        double b33 = (a11*a22 - a21*a12)/datv;
                        double vc_matrix = (l_x*(b11*l_x+b21*l_y+b31*l_z) + l_y*(b12*l_x+b22*l_y+b32*l_z) + l_z*(b13*l_x+b23*l_y+b33*l_z));
                        double f_pcar = exp(-0.5*vc_matrix);
                        //LOGF(info, "PCAX-CloX: %lf, PCAY-CloY: %lf, PCAZ-CloZ: %lf", l_x,l_y,l_z);
                        //if(datv!=0) LOGF(info, "datv is not 0!: %f", datv);
                        if(datv==0){
                          LOGF(info, "datv is 0");
                          LOGF(info, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", a11,a12,a13,a21,a22,a23,a31,a32,a33);
                        }else{
                          LOGF(info, "vc_matrix: %lf", vc_matrix);
                          vc_sum += f_pcar;
                          vc2_sum += f_pcar*f_pcar;
                        }
                      }
                    }
                  }
                  double p_pcar = vc_sum - vc2_sum/vc_sum;
                  //LOGF(info, "vc2: %lf, vc: %lf", vc2_sum, vc_sum);
                  LOGF(info, "p_PCAr: %lf", p_pcar);

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
    adaptAnalysisTask<Test>(cfgc),
    adaptAnalysisTask<MyAnalysisTask>(cfgc)};
}