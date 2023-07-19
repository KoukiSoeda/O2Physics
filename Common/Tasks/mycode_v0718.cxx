/// \author Koki Soeda
/// \since 04/07/2023

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
    const AxisSpec axisZ(nBinsZ, -60.0005, 50.0005, "z(cm)");
    const AxisSpec axisR(nBinsR, -1.0005, 10.0005, "r(cm)");
		const AxisSpec axisDist(nBinsDist, -0.0005, 10.0005, "z(cm)");

    histos.add("TrackType", "TrackType", kTH1F, {axistrackType});
    histos.add("Distance_D0_z", "Distance_D0_z", kTH1F, {axisDist});
    histos.add("PredictZ", "PredictZ", kTH1F, {axisZ});
    histos.add("PCAZ_Answer", "PCAZ_Answer", kTH1F, {axisZ});
    histos.add("Diff_D0Ver_ColVer", "Diff_D0Ver_ColVer", kTH1F, {axisZ});
    histos.add("Diff_D0Ver_ColVer_5over", "Diff_D0Ver_ColVer_5over", kTH1F, {axisZ});
    histos.add("PredictX", "PredictX", kTH1F, {axisZ});
    histos.add("PCAX_Answer", "PCAX_Answer", kTH1F, {axisZ});
    histos.add("PredictY", "PredictY", kTH1F, {axisZ});
    histos.add("PCAY_Answer", "PCAY_Answer", kTH1F, {axisZ});
    histos.add("MCvertex_Estvertex", "MCvertex_Estvertex", kTH1F, {axisR});
    histos.add("PCAtoEst_R", "PCAtoEst_R", kTH1F, {axisR});
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
        if(fabs(mcParticle_fwd.pdgCode())!=13 || !mcParticle_fwd.has_mothers()) continue;
        auto fwdColId = fwdtrack.collisionId();
        auto mcMoms = mcParticle_fwd.mothers_as<aod::McParticles>();
        auto mcMom = mcMoms.back();
        auto mcMom_f = mcMoms.front();
        //auto momindex = mcMom.globalIndex();

        if(mcMom_f.pdgCode()==mcMom.pdgCode() && fabs(mcMom.pdgCode())==421){
          auto Daughters = mcMom.daughters_as<aod::McParticles>();
          int dc = 0;
          int mccolid, mucolid, mcmuid, mupdg, mckid, kpdg;
          float Col_x, Col_y, Col_z, mudcaX, mudcaY, mumftX, mumftY, mumftZ, kdcaX, kdcaY, kmftX, kmftY, kmftZ, mumcX, mumcY, mumcZ, kmcZ;
          Col_x = collision.posX();
          Col_y = collision.posY();
          Col_z = collision.posZ();
          for(auto& Daughter : Daughters){
            if(fabs(Daughter.pdgCode())==13){
              auto muindex = Daughter.globalIndex();
              for(auto& mfttrack : mfttracks){
                if(mfttrack.has_collision() && mfttrack.has_mcParticle()){
                  if(mfttrack.collisionId()!=fwdColId) continue;
                  auto mcParticle_mft = mfttrack.mcParticle();
                  if(mfttrack.mcParticleId()==muindex){
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
                    mccolid = mcParticle_mft.mcCollisionId();
                    mcmuid = mfttrack.mcParticleId();
                    mupdg = mcParticle_mft.pdgCode();
                    dc++;
                  }
                }
              }
            }else if(fabs(Daughter.pdgCode())==321){
              auto kindex = Daughter.globalIndex();
              for(auto& mfttrack2 : mfttracks){
                if(mfttrack2.has_collision() && mfttrack2.has_mcParticle()){
                  if(mfttrack2.collisionId()!=fwdColId) continue;
                  auto mcParticle_mft2 = mfttrack2.mcParticle();
                  if(mfttrack2.mcParticleId()==kindex){
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
                    mckid = mfttrack2.mcParticleId();
                    kpdg = mcParticle_mft2.pdgCode();
                    dc++;
                  }
                }
              }
            }
          }
          if(collision.mcCollisionId()==mccolid){
            if(dc==2 && mumcZ==kmcZ){
              histos.fill(HIST("Distance_D0_z"), fabs(mcParticle_fwd.vz()-mcMom.vz()));
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
              if(fabs(mumcZ-collision.mcCollision().posZ())>5){
                histos.fill(HIST("Diff_D0Ver_ColVer_5over"), mcMom.vz()-collision.mcCollision().posZ());
                //LOGF(info, "Collision GlobaIndex: %d, Track CollisionID: %d", collision.globalIndex(), mfttrack.collisionId());
              }else{
                histos.fill(HIST("Diff_D0Ver_ColVer"), mcMom.vz()-collision.mcCollision().posZ());
              }
              histos.fill(HIST("PredictX"), pcaX-Col_x);
              histos.fill(HIST("PCAX_Answer"), mumcX-Col_x);
              histos.fill(HIST("PredictY"), pcaY-Col_y);
              histos.fill(HIST("PCAY_Answer"), mumcY-Col_y);
              histos.fill(HIST("MCvertex_Estvertex"), r_org);
              mcpairtable(mucolid,mcmuid,mupdg,mckid,kpdg,Col_x,Col_y,Col_z,mumcX,mumcY,mumcZ,pcaX-Col_x,pcaY-Col_y,pcaZ-Col_z,r);
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
  Configurable<int> nBinsCut{"nBinsCut", 2011, "N bins in cut"};

  void init(InitContext const&){
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisEvent{1000000, 0, 1000000};
    const AxisSpec axisTest{1, 0, 1, ""};
    const AxisSpec axisPt{nBinsPt, -0.005, 25.005, "p_{T}"};
    const AxisSpec axisR(nBinsR, -1.0005, 5.0005, "r(cm)");
    const AxisSpec axisZ(nBinsZ, -60.0005, 50.0005, "z(cm)");
    const AxisSpec axisCutZ(nBinsCut, -20.05, 0.05, "z(cm)");
    const AxisSpec axisParticle(101, -0.5, 100.5, "N");

    histos.add("Est_PCAZ", "Est_PCAZ", kTH1F, {axisZ});
    histos.add("Est_PCAZ_Correct", "Est_PCAZ_Correct", kTH1F, {axisZ});
    histos.add("CutZ_Correct", "CutZ_Correct", kTH1F, {axisCutZ});
    histos.add("EstZ_Correct", "EstZ_Correct", kTH1F, {axisZ});
    histos.add("Correct_pt", "Correct_pt", kTH1F, {axisPt});
    histos.add("Correct_r_dist", "Correct_r_dist", kTH1F, {axisR});
    histos.add("Incorrect_r_dist", "Incorrect_r_dist", kTH1F, {axisR});
    histos.add("Actual_r", "Actual_r", kTH1F, {axisR});
    histos.add("PosDiff_x", "PosDiff_x", kTH1F, {axisZ});
    histos.add("PosDiff_y", "PosDiff_y", kTH1F, {axisZ});
    histos.add("PosDiff_z", "PosDiff_z", kTH1F, {axisZ});
    histos.add("PosDiff", "PosDiff", kTH1F, {axisZ});
    histos.add("InnerParticles", "InnerParticles", kTH1F, {axisParticle});
    histos.add("Total_pt", "Total_pt", kTH1F, {axisPt});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
							 soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
							 soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
							 aod::McParticles const& MCparticle,
							 aod::McCollisions const&,
               aod::MCPair const& truepairs)
  {
    if(collision.has_mcCollision()){
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

        int mucolid, mcmuid, mckid_org, kpdg, estmumomid, estkmomid, mumom, kmom_org, storemuid, estmck;
        float mupt, Col_x, Col_y, Col_z, mudcaX, mudcaY, mumftX, mumftY, mumftZ, kdcaX, kdcaY, kmftX, kmftY, kmftZ, mumcX, mumcY, mumcZ, kmcZ, estx, esty, estz;
        float dist_2track = 10000;
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
            //mupt = fwdtrack.pt();
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
                
                float fcut = 10;
                int colID_ref;
                for(auto& mfttrack2 : mfttracks){
                  if(mfttrack2.has_collision() && mfttrack2.has_mcParticle()){
                    if(mfttrack2.collisionId()==fwdColId && mfttrack2.mcParticleId()!=mcmuid){
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

                      if(pcaz_org<=0 && pcaz_org>=-fcut){
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
                          colID_ref = mfttrack2.collisionId();
                        }
                      }/*else{
                        for(int i=0; i<100001; i++){
                          float h = 1 - 0.00001*i;
                          auto x0 = mumftX + h*Nax;
                          auto y0 = mumftY + h*Nay;
                          auto z0 = mumftZ + h*Naz;
                          auto kx0 = kmftX + ((x0-kmftX)*(Ncx/k_unit)+(y0-kmftY)*(Ncy/k_unit)+(z0-kmftZ)*(Ncz/k_unit))*(Ncx/k_unit);
                          auto ky0 = kmftY + ((y0-kmftY)*(Ncx/k_unit)+(y0-kmftY)*(Ncy/k_unit)+(z0-kmftZ)*(Ncz/k_unit))*(Ncy/k_unit);
                          auto kz0 = kmftZ + ((z0-kmftZ)*(Ncx/k_unit)+(y0-kmftY)*(Ncy/k_unit)+(z0-kmftZ)*(Ncz/k_unit))*(Ncz/k_unit);

                          if((((z0+kz0)/2)-Col_z)<=0 && (((z0+kz0)/2)-Col_z)>=-fcut){
                            auto r_xyz_org = sqrt(pow(x0-kx0,2)+pow(y0-ky0,2)+pow(z0-kz0,2));
                            if(r_xyz>=r_xyz_org){
                              r_xyz = r_xyz_org;
                              pcax = (x0+kx0)/2 - Col_x;
                              pcay = (y0+ky0)/2 - Col_y;
                              pcaz = (z0+kz0)/2 - Col_z;
                              kmom = kmom_org;
                              mckid = mckid_org;
                            }
                          }
                        }
                      }*/
                    }
                  }
                }

                histos.fill(HIST("Est_PCAZ"), estz);
                histos.fill(HIST("Total_pt"), mupt);
                //histos.fill(HIST("InnerParticles"), npar);

                if(estmumomid==estkmomid && truepair.kid()==estmck){
                  histos.fill(HIST("Correct_pt"), mupt);
                  histos.fill(HIST("Correct_r_dist"), dist_2track);
                  histos.fill(HIST("Est_PCAZ_Correct"), estz);
                  histos.fill(HIST("EstZ_Correct"), estz-(truepair.secverz()-Col_z));
                  histos.fill(HIST("CutZ_Correct"), -fcut);
                  histos.fill(HIST("InnerParticles"), npar);
                }else{
                  histos.fill(HIST("Actual_r"), truepair.r());
                  histos.fill(HIST("PosDiff_x"), estx-truepair.pcax());
                  histos.fill(HIST("PosDiff_y"), esty-truepair.pcay());
                  histos.fill(HIST("PosDiff_z"), estz-truepair.pcaz());
                  if(fabs(estx-truepair.pcax())<0.1 && fabs(esty-truepair.pcay())<0.1){
                    if(fabs(estz-truepair.pcaz())<0.1) histos.fill(HIST("PosDiff"), estz-truepair.pcaz());
                  }
                  histos.fill(HIST("Incorrect_r_dist"), dist_2track);
                }
                storemuid = mcmuid;
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