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

struct MyAnalysisTask{

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsZ{"nBinsZ", 5000, "N bins in Z axis"};
  Configurable<int> nBinsR{"nBinsR", 60010, "N bins in R"};
	Configurable<int> nBinsDCA{"nBinsDCA", 10001, "N bins in DCA"};
	Configurable<int> nBinsDCAT{"nBinsDCAT", 1000, "N bins in DCAT"};
	Configurable<int> nBinsDist{"nBinsDist", 400010, ""};
	Configurable<int> fwdTrackType{"fwdTrackType", 0, "N TrackType in fwd"};

  void init(InitContext const&){
    //const AxisSpec axisVertex{nBinsVer, -30, +30, "cm"};
		const AxisSpec axistrackType{6, -0.5, 5.5, "TrackType"};
		const AxisSpec axisCounter{1, 0, +1, ""};
		const AxisSpec axisChi2(1005, -4.5, +1000.5, "");
		const AxisSpec axisPDG(19, -0.5, +18.5, "");
    const AxisSpec axisZ(nBinsZ, -60.0005, 50.0005, "z(cm)");
    const AxisSpec axisMFT(80010, -80.005, 0.005, "z(cm)");
    const AxisSpec axisR(nBinsR, -1.0005, 5.0005, "r(cm)");
		const AxisSpec axisDist(nBinsDist, -0.0005, 40.0005, "z(cm)");
		const AxisSpec axisDCAX(nBinsDCA, -0.50005, 0.50005, "x(cm)");
		const AxisSpec axisDCAY(nBinsDCA, -0.50005, 0.50005, "y(cm)");
    const AxisSpec axisPt{nBinsPt, 0, 25, "p_{T}"};

    histos.add("TrackType", "TrackType", kTH1F, {axistrackType});
    histos.add("Distance_D0_z", "Distance_D0_z", kTH1F, {axisDist});
    histos.add("PredictZ", "PredictZ", kTH1F, {axisZ});
    histos.add("PCAZ_Answer", "PCAZ_Answer", kTH1F, {axisZ});
    histos.add("Z_diff", "Z_diff", kTH1F, {axisZ});
    histos.add("xy_pca", "xy_pca", kTH1F, {axisR});
    histos.add("DCAXY_correct", "DCAXY_correct", kTH2F, {axisDCAX, axisDCAY});
    histos.add("DCAXY_est", "DCAXY_est", kTH2F, {axisDCAX, axisDCAY});
		histos.add("PredictZ_bg", "PredictZ_bg", kTH1F, {axisZ});
    histos.add("xy_incorrect", "xy_incorrect", kTH1F, {axisR});
    histos.add("xy_actual", "xy_actual", kTH1F, {axisR});
    histos.add("DCAXY_incorrect", "DCAXY_incorrect", kTH2F, {axisDCAX, axisDCAY});
    histos.add("predictPDG", "predictPDG", kTH1F, {axisPDG});
		histos.add("counter1", "counter1", kTH1F, {axisCounter});
		histos.add("counter2", "counter2", kTH1F, {axisCounter});

    auto pdgcode = histos.get<TH1>(HIST("predictPDG"));
		auto* x2 = pdgcode->GetXaxis();
		x2->SetBinLabel(1, "Primary"); // 1 - 37
		x2->SetBinLabel(2, "Pion0"); // 111
		x2->SetBinLabel(3, "Pion±"); // 211
		x2->SetBinLabel(4, "pho"); // 113
		x2->SetBinLabel(5, "K_l"); // 130
		x2->SetBinLabel(6, "Eta"); // 221
		x2->SetBinLabel(7, "Omega"); // 223
		x2->SetBinLabel(8, "K_s"); // 310
		x2->SetBinLabel(9, "K*0(892)"); // 313
		x2->SetBinLabel(10, "K±"); // 321
		x2->SetBinLabel(11, "K*±(892)"); // 323
		x2->SetBinLabel(12, "Eta_prim"); // 331
		x2->SetBinLabel(13, "Phi"); // 333
		x2->SetBinLabel(14, "D±"); // 411
		x2->SetBinLabel(15, "D*±"); // 413
		x2->SetBinLabel(16, "D0"); // 421
		x2->SetBinLabel(17, "D_s±"); // 431
		x2->SetBinLabel(18, "Beauty"); // 500-599
		x2->SetBinLabel(19, "Baryon"); // 1000 -
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
        auto momindex = mcMom.globalIndex();

        int dc = 0;
        int mucolid;
        float Col_z, mudcaX, mudcaY, mumftX, mumftY, mumftZ, kdcaX, kdcaY, kmftX, kmftY, kmftZ, mumcZ, kmcZ;
        if(fabs(mcMom.pdgCode())!=421) continue;
        for(auto& mfttrack : mfttracks){
          if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
          if(mfttrack.collisionId()!=fwdColId) continue;
          auto mcParticle_mft = mfttrack.mcParticle();
          if(!mcParticle_mft.has_mothers()) continue;
          auto mft_Moms = mcParticle_mft.mothers_as<aod::McParticles>();
          auto mft_Mom = mft_Moms.back();
          Col_z = collision.posZ();

          if(fwdtrack.mcParticleId()==mfttrack.mcParticleId()){
            if(mft_Mom.globalIndex()!=momindex) continue;
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
            mumcZ = mcParticle_mft.vz();
            mucolid = mfttrack.collisionId();
            dc++;
          }
          if(mft_Mom.globalIndex()==momindex && fabs(mcParticle_mft.pdgCode())==321){
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
            kmcZ = mcParticle_mft.vz();
            dc++;
          }
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
            float pcaZ = (pre_muz+pre_kz)/2;
            float r_actual = sqrt(pow(pre_kx-pre_mux, 2)+pow(pre_ky-pre_muy, 2)+pow(pre_kz-pre_muz, 2))/2;
            if(pcaZ-Col_z>0) cout << "Mother: " << mcMom.pdgCode() << " vz: " << mcMom.vz() << ", ColZ: " << Col_z << ", D1: " << mumcZ << ", D2: " << kmcZ << ", PCAZ: " << pcaZ << endl; 
            histos.fill(HIST("PredictZ"), pcaZ-Col_z);
            histos.fill(HIST("PCAZ_Answer"), mumcZ-Col_z);
            histos.fill(HIST("Z_diff"), (pcaZ-Col_z)-(mumcZ-Col_z));
            histos.fill(HIST("xy_pca"), r_actual);
            dc = 0;
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyAnalysisTask>(cfgc)};
}