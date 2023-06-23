/// \author Koki Soeda
/// \since 12/06/2023

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

struct MyAnalysisTask {
	// Histogram registry: an object to hold your histograms
	HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

	Configurable<int> nBinsVer{"nBinsVer", 300, "N bins in Collisino Vertex histo"};
	Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
	Configurable<int> nBinsEta{"nBinsEta", 100, ""};
	Configurable<int> nBinsPhi{"nBinsPhi", 3000, ""};
	Configurable<int> nBinsXY{"nBinsXY", 10000, "N bins in X and Y axis"};
	Configurable<int> nBinsZ{"nBinsZ", 5000, "N bins in Z axis"};
  Configurable<int> nBinsR{"nBinsR", 20010, "N bins in R"};
	Configurable<int> nBinsDCA{"nBinsDCA", 200010, "N bins in DCA"};
	Configurable<int> nBinsDCAT{"nBinsDCAT", 1000, "N bins in DCAT"};
	Configurable<int> nBinsDist{"nBinsDist", 800010, ""};
	Configurable<int> fwdTrackType{"fwdTrackType", 0, "N TrackType in fwd"};
	Configurable<int> fwdTrackType2{"fwdTrackType2", 2, "N TrackType in fwd2"};
	Configurable<float> zaxisMaxCut{"zaxisMaxCut", 40, "MaxCutForPCA"};
	Configurable<float> zaxisMinCut{"zaxisMinCut", -40, "MinCutForPCA"};

	void init(InitContext const&)
	{
		// define axes you want to use
		const AxisSpec axisVertex{nBinsVer, -30, +30, "cm"};
		const AxisSpec axistrackType{6, -0.5, 5.5, "TrackType"};
		const AxisSpec axisCounter{1, 0, +1, ""};
		const AxisSpec axisChi2(1005, -4.5, +1000.5, "");
		const AxisSpec axisPDG(19, -0.5, +18.5, "");
		const AxisSpec axisPDGv2(1001, -0.5, +1000.5, "");
		const AxisSpec axisX(nBinsXY, -10, 10, "x(cm)");
		const AxisSpec axisY(nBinsXY, -10, 10, "y(cm)");
		const AxisSpec axisZ(nBinsZ, -50.0005, 50.0005, "z(cm)");
    const AxisSpec axisR(nBinsR, -0.0005, 2.0005, "r(cm)");
		const AxisSpec axisDist(nBinsDist, -40.0005, 40.0005, "z(cm)");
		const AxisSpec axisDCAX(nBinsDCA, -1.00005, 1.00005, "x(cm)");
		const AxisSpec axisDCAY(nBinsDCA, -1.00005, 1.00005, "y(cm)");
		const AxisSpec axisDCAT(nBinsDCAT, -0.01, +1, "cm");
		const AxisSpec axisEta{nBinsEta, -10, +10, "#eta"};
		const AxisSpec axisPhi{nBinsPhi, -6.5, +6.5, "#phi"};
		const AxisSpec axisPt{nBinsPt, 0, 25, "p_{T}"};
		const AxisSpec axisPro{200, -0.5, 1.5, ""};
		const AxisSpec axisClu{101, -0.5, 100.5, ""};

		// create histograms
		histos.add("mcCollisionPos", "CollisionPos", kTH1F, {axisVertex});
		histos.add("TrackType", "TrackType", kTH1F, {axistrackType});
		histos.add("muonMomPDG", "muonMomPDG", kTH1F, {axisPDG});
    histos.add("Distance_D0_z", "Distance_D0_z", kTH1F, {axisDist});
    histos.add("D0Decay_muK", "D0Decay_muK", kTH1F, {axisZ});
		histos.add("PredictZ", "PredictZ", kTH1F, {axisZ});
    histos.add("chi2MCHMFT_correct", "chi2MCHMFT_correct", kTH1F, {axisChi2});
    histos.add("xy_correct", "xy_correct", kTH1F, {axisR});
    histos.add("DCAXY_correct", "DCAXY_correct", kTH2F, {axisDCAX, axisDCAY});
		histos.add("PredictZ_bg", "PredictZ_bg", kTH1F, {axisZ});
    histos.add("chi2MCHMFT_incorrect", "chi2MCHMFT_incorrect", kTH1F, {axisChi2});
    histos.add("xy_incorrect", "xy_incorrect", kTH1F, {axisR});
    histos.add("xy_actual", "xy_actual", kTH1F, {axisR});
    histos.add("DCAXY_incorrect", "DCAXY_incorrect", kTH2F, {axisDCAX, axisDCAY});
		histos.add("predictPDG", "predictPDG", kTH1F, {axisPDG});
		histos.add("counter1", "counter1", kTH1F, {axisCounter});
		histos.add("counter2", "counter2", kTH1F, {axisCounter});

		auto hmpdgcode = histos.get<TH1>(HIST("muonMomPDG"));
		auto* x1 = hmpdgcode->GetXaxis();
		x1->SetBinLabel(1, "Primary"); // 1 - 37
		x1->SetBinLabel(2, "Pion0"); // 111
		x1->SetBinLabel(3, "Pion±"); // 211
		x1->SetBinLabel(4, "pho"); // 113
		x1->SetBinLabel(5, "K_l"); // 130
		x1->SetBinLabel(6, "Eta"); // 221
		x1->SetBinLabel(7, "Omega"); // 223
		x1->SetBinLabel(8, "K_s"); // 310
		x1->SetBinLabel(9, "K*0(892)"); // 313
		x1->SetBinLabel(10, "K±"); // 321
		x1->SetBinLabel(11, "K*±(892)"); // 323
		x1->SetBinLabel(12, "Eta_prim"); // 331
		x1->SetBinLabel(13, "Phi"); // 333
		x1->SetBinLabel(14, "D±"); // 411
		x1->SetBinLabel(15, "D*±"); // 413
		x1->SetBinLabel(16, "D0"); // 421
		x1->SetBinLabel(17, "D_s±"); // 431
		x1->SetBinLabel(18, "Beauty"); // 500-599
		x1->SetBinLabel(19, "Baryon"); // 1000 -

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
			//Informatin of MC collision
			auto mcCol_z = collision.mcCollision().posZ();
			histos.fill(HIST("mcCollisionPos"), mcCol_z);

			//Information of FwdTracks
			for(auto& fwdtrack : fwdtracks){
        if(collision.posZ()>10 && collision.posZ()<-10) continue;
				histos.fill(HIST("TrackType"), fwdtrack.trackType());

				if(!fwdtrack.has_mcParticle()) continue;
				auto mcParticle_fwd = fwdtrack.mcParticle();
				auto fwdTrackPDG = mcParticle_fwd.pdgCode();
        auto fwdID = mcParticle_fwd.globalIndex();

				//Informaion of D0 meson
				float mudcaX, mudcaY, mumftX, mumftY, mumftZ, kdcaX, kdcaY, kmftX, kmftY, kmftZ;
				float mu_phi, mu_eta, k_phi, k_eta;
				if(fwdtrack.trackType()==fwdTrackType){
					auto chi2GMT = fwdtrack.chi2MatchMCHMFT();

					if(chi2GMT==-1) continue;
					if(fabs(fwdTrackPDG)==13){
						auto a_ID = mcParticle_fwd.globalIndex();

						if(mcParticle_fwd.has_mothers()){
              auto mcMoms = mcParticle_fwd.mothers_as<aod::McParticles>();
              auto mcMom = mcMoms.back(); //The last mother of global muon track
							auto mcMomPDG = mcMom.pdgCode();
							auto Daughters = mcMom.daughters_as<aod::McParticles>();
              auto fwdcollID = fwdtrack.collisionId();

							if(fabs(mcMomPDG)==421){
                if(mcMom.producedByGenerator()==true){ //Selected the particle produced by the generator
                  float pcaZ, pcaX, pcaY;
                  int preidPDG, idPDG, k_mom, mu_mom;
                  int64_t mftID;
                  int daughter_count = 0;
                  float secver_z;
                  float predictdist = 1000;
                  float predictX, predictY, predictZ;
                  
                  for(auto Daughter : Daughters){

                    if(fabs(Daughter.pdgCode())==13){//muon
                      auto mu_ID = Daughter.globalIndex();
                      if(mu_ID==a_ID){ //Check the index
                        Double_t verZ = collision.posZ();
                        double mftchi2 = fwdtrack.chi2();
                        SMatrix5 mftpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
                        vector<double> mftv1;
                        SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                        o2::track::TrackParCovFwd mftpars1{fwdtrack.z(), mftpars, mftcovs, mftchi2};
                        mftpars1.propagateToZlinear(verZ);

                        mudcaX = mftpars1.getX() - collision.posX();
                        mudcaY = mftpars1.getY() - collision.posY();
                        mumftX = fwdtrack.x();
                        mumftY = fwdtrack.y();
                        mumftZ = fwdtrack.z();
                        mu_phi = fwdtrack.phi();
                        mu_eta = fwdtrack.eta();
                        secver_z = Daughter.vz() - collision.posZ();

                        histos.fill(HIST("Distance_D0_z"), Daughter.vz()-mcMom.vz()); //fabs(Daughter.vz()-mcMom.vz()));
                        mu_mom = Daughter.mothers_as<aod::McParticles>().back().globalIndex();
                        if(secver_z<0) daughter_count++;
                      }
                    }
                    if(fabs(Daughter.pdgCode())==321){//kaon
                      auto k_ID = Daughter.globalIndex();
                      for(auto& mfttrack : mfttracks){
                        if(!mfttrack.has_collision() && !mfttrack.has_mcParticle()) continue;
                        auto mcParticle_mft = mfttrack.mcParticle();
                        mftID = mcParticle_mft.globalIndex();
                        if(mftID==k_ID){
                          Double_t verZ = collision.posZ();
                          double mftchi2 = mfttrack.chi2();
                          SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                          vector<double> mftv1;
                          SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                          o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
                          mftpars1.propagateToZlinear(verZ);
                          
                          k_phi = mfttrack.phi();
                          k_eta = mfttrack.eta();
                          kmftX = mfttrack.x();
                          kmftY = mfttrack.y();
                          kmftZ = mfttrack.z();
                          kdcaX = mftpars1.getX() - collision.posX();
                          kdcaY = mftpars1.getY() - collision.posY();
                          k_mom = Daughter.mothers_as<aod::McParticles>().back().globalIndex();
                          daughter_count++;
                        }
                      }
                    }

                  }

                  if(daughter_count==2){
                    histos.fill(HIST("D0Decay_muK"), secver_z);
                  
                    //Actual PCA
                    auto Nax = mudcaX - mumftX;
                    auto Nay = mudcaY - mumftY;
                    auto Naz = collision.posZ() - mumftZ;
                    auto Ncx = kdcaX - kmftX;
                    auto Ncy = kdcaY - kmftY;
                    auto Ncz = collision.posZ() - kmftZ;
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
                    float r_actual = sqrt(pow(pre_kx-pre_mux, 2)+pow(pre_ky-pre_muy, 2)+pow(pre_kz-pre_muz, 2))/2;

                    //Scan all mfttrack for calculate the PCA
                    int64_t predictID, preID;
                    float dist_2track, r_can, r_xy, actual_dcaX, actual_dcaY;
                    for(auto& mfttrack : mfttracks){
                      if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
                      auto mcParticle_mft = mfttrack.mcParticle();
                      auto mftcollID = mfttrack.collisionId();
                      if(mftcollID!=fwdcollID) continue;
                      if(mcParticle_mft.globalIndex()==fwdID) continue;
                      
                      Double_t verZ = collision.posZ();
                      double mftchi2 = mfttrack.chi2();
                      SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                      vector<double> mftv1;
                      SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                      o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
                      mftpars1.propagateToZlinear(verZ);

                      auto mftX = mfttrack.x();
                      auto mftY = mfttrack.y();
                      auto mftZ = mfttrack.z();
                      auto dcaX = mftpars1.getX() - collision.posX();
                      auto dcaY = mftpars1.getY() - collision.posY();

                      auto Nax = mudcaX - mumftX;
                      auto Nay = mudcaY - mumftY;
                      auto Naz = collision.posZ() - mumftZ;
                      auto Ncx = dcaX - mftX;
                      auto Ncy = dcaY - mftY;
                      auto Ncz = collision.posZ() - mftZ;
                      auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
                      auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
                      auto A3 = (mumftX-mftX)*Nax + (mumftY-mftY)*Nay + (mumftZ-mftZ)*Naz;
                      auto B1 = A2;
                      auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
                      auto B3 = (mftX-mumftX)*Ncx + (mftY-mumftY)*Ncy + (mftZ-mumftZ)*Ncz;
                      auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
                      auto s = -((A2*t+A3)/A1);
                      float pre_mux = mumftX + s*Nax;
                      float pre_muy = mumftY + s*Nay;
                      float pre_muz = mumftZ + s*Naz;
                      float pre_kx = mftX + t*Ncx;
                      float pre_ky = mftY + t*Ncy;
                      float pre_kz = mftZ + t*Ncz;

                      if(((pre_muz+pre_kz)/2)>=mumftZ && ((pre_muz+pre_kz)/2)<=collision.posZ()){
                        dist_2track = sqrt(pow(pre_kx-pre_mux, 2)+pow(pre_ky-pre_muy, 2)+pow(pre_kz-pre_muz, 2));
                        pcaX = (pre_mux+pre_kx)/2;
                        pcaY = (pre_muy+pre_ky)/2;
                        pcaZ = (pre_muz+pre_kz)/2;
                        r_can = sqrt(pow(pcaX-pre_mux,2)+pow(pcaY-pre_muy,2)+pow(pcaZ-pre_muz,2));
                        preidPDG = mcParticle_mft.pdgCode();
                        preID = mcParticle_mft.globalIndex();
                      }else{
                        for(int i=0; i<=500000 ; i++){
                          s = 1 - 0.000001*i;
                          auto x0 = mumftX + s*Nax;
                          auto y0 = mumftY + s*Nay;
                          auto z0 = mumftZ + s*Naz;
                          auto paraT = (Ncx*(x0-mftX)+Ncy*(y0-mftY)+Ncz*(z0-mftZ))/B2;
                          auto d = sqrt(pow(x0-(mftX+paraT*Ncx),2)+pow(y0-(mftY+paraT*Ncy),2)+pow(z0-(mftZ+paraT*Ncz),2));
                          if(d<=dist_2track){
                            dist_2track = d;
                            pcaX = (x0+(mftX+paraT*Ncx))/2;
                            pcaY = (y0+(mftY+paraT*Ncy))/2;
                            pcaZ = (z0+(mftZ+paraT*Ncz))/2;
                            r_can = sqrt(pow(pcaX-x0,2)+pow(pcaY-y0,2));
                            preidPDG = mcParticle_mft.pdgCode();
                            preID = mcParticle_mft.globalIndex();
                          }
                        }
                      }

                      if(predictdist>dist_2track){
                        predictdist = dist_2track;
                        predictX = pcaX;
                        predictY = pcaY;
                        predictZ = pcaZ;
                        r_xy = r_can;
                        idPDG = preidPDG;
                        predictID = preID;
                        actual_dcaX = dcaX;
                        actual_dcaY = dcaY;
                      }
                    }

                    //histos.fill(HIST("predictPDG"), fabs(idPDG));
                    if(fabs(idPDG) < 38) histos.fill(HIST("predictPDG"), 0.0, 1);
                    if(fabs(idPDG) == 111) histos.fill(HIST("predictPDG"), 1.0, 1);
                    if(fabs(idPDG) == 211) histos.fill(HIST("predictPDG"), 2.0, 1);
                    if(fabs(idPDG) == 113 ) histos.fill(HIST("predictPDG"), 3.0, 1);
                    if(fabs(idPDG) == 130 ) histos.fill(HIST("predictPDG"), 4.0, 1);
                    if(fabs(idPDG) == 221 ) histos.fill(HIST("predictPDG"), 5.0, 1);
                    if(fabs(idPDG) == 223 ) histos.fill(HIST("predictPDG"), 6.0, 1);
                    if(fabs(idPDG) == 310 ) histos.fill(HIST("predictPDG"), 7.0, 1);
                    if(fabs(idPDG) == 313 ) histos.fill(HIST("predictPDG"), 8.0, 1);
                    if(fabs(idPDG) == 321 ) histos.fill(HIST("predictPDG"), 9.0, 1);
                    if(fabs(idPDG) == 323 ) histos.fill(HIST("predictPDG"), 10.0, 1);
                    if(fabs(idPDG) == 333 ) histos.fill(HIST("predictPDG"), 12.0, 1);
                    if(fabs(idPDG) == 411 ) histos.fill(HIST("predictPDG"), 13.0, 1);
                    if(fabs(idPDG) == 413 ) histos.fill(HIST("predictPDG"), 14.0, 1);
                    if(fabs(idPDG) == 431 ) histos.fill(HIST("predictPDG"), 16.0, 1);
                    if(fabs(idPDG) > 499 && fabs(idPDG) < 600) histos.fill(HIST("predictPDG"), 17.0, 1);
                    if(fabs(idPDG) > 999 ) histos.fill(HIST("predictPDG"), 18.0, 1);

                    if(k_mom==mu_mom){
                      histos.fill(HIST("counter1"), 0.5);
                      if(mftID==predictID){ // Correct
                        histos.fill(HIST("counter2"), 0.5);
                        histos.fill(HIST("PredictZ"), predictZ-collision.posZ());
                        histos.fill(HIST("chi2MCHMFT_correct"), chi2GMT);
                        histos.fill(HIST("xy_correct"), r_xy);
                        histos.fill(HIST("DCAXY_correct"), actual_dcaX, actual_dcaY);

                      }else{ // Incorrect
                        histos.fill(HIST("PredictZ_bg"), predictZ-collision.posZ());
                        histos.fill(HIST("chi2MCHMFT_incorrect"), chi2GMT);
                        histos.fill(HIST("xy_incorrect"), r_xy);
                        histos.fill(HIST("xy_actual"), r_actual);
                        histos.fill(HIST("DCAXY_incorrect"), actual_dcaX, actual_dcaY);

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
    adaptAnalysisTask<MyAnalysisTask>(cfgc)};
}
