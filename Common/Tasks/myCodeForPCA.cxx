/// \author Koki Soeda
/// \since 07/06/2023

#include <iostream>
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
	Configurable<int> nBinsDCA{"nBinsDCA", 5000, "N bins in DCA"};
	Configurable<int> nBinsDCAT{"nBinsDCAT", 1000, "N bins in DCAT"};
	Configurable<int> nBinsDist{"nBinsDist", 500, ""};
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
		const AxisSpec axisChi2(505, -4.5, +500.5, "");
		const AxisSpec axisPDG(19, -0.5, +18.5, "");
		const AxisSpec axisPDGv2(1001, -0.5, +1000.5, "");
		const AxisSpec axisX(nBinsXY, -10, 10, "x(cm)");
		const AxisSpec axisY(nBinsXY, -10, 10, "y(cm)");
		const AxisSpec axisZ(nBinsZ, -50.0005, 50.0005, "z(cm)");
		const AxisSpec axisDist(nBinsDist, -0.0005, 20.0005, "z(cm)");
		const AxisSpec axisDCAX(nBinsDCA, -1, 1, "x(cm)");
		const AxisSpec axisDCAY(nBinsDCA, -1, 1, "y(cm)");
		const AxisSpec axisDCAT(nBinsDCAT, -0.01, +1, "cm");
		const AxisSpec axisEta{nBinsEta, -10, +10, "#eta"};
		const AxisSpec axisPhi{nBinsPhi, -6.5, +6.5, "#phi"};
		const AxisSpec axisPt{nBinsPt, 0, 25, "p_{T}"};
		const AxisSpec axisPro{200, -0.5, 1.5, ""};
		const AxisSpec axisClu{101, -0.5, 100.5, ""};

		// create histograms
		histos.add("CollisionPos", "CollisionPos", kTH1F, {axisVertex});
		histos.add("TrackType", "TrackType", kTH1F, {axistrackType});
		histos.add("FwdParticleCounter", "FwdParticleCounter", kTH1F, {axisCounter});
		histos.add("FwdMuonCounter", "FwdMuonCounter", kTH1F, {axisCounter});
		histos.add("GMT_pT", "GMT_pT", kTH1F, {axisPt});
		histos.add("muonMomPDG", "muonMomPDG", kTH1F, {axisPDG});
		histos.add("D0DecayX", "D0DecayX", kTH1F, {axisX});
		histos.add("D0DecayY", "D0DecayY", kTH1F, {axisY});
		histos.add("D0DecayZ", "D0DecayZ", kTH1F, {axisZ});
		histos.add("Distance_D0_z", "Distance_D0_z", kTH1F, {axisDist});
		histos.add("D02mu_DCA_X", "DCA_X", kTH1F, {axisDCAX});
		histos.add("D02mu_DCA_Y", "DCA_Y", kTH1F, {axisDCAY});
		histos.add("D02mu_DCAT", "D02mu_DCAT", kTH1F, {axisDCAT});
		histos.add("D0_pT", "D0_pT", kTH1F, {axisPt});
		histos.add("D02mu_eta", "D02mu_eta", kTH1F, {axisEta});
		histos.add("D02mu_phi", "D02mu_phi", kTH1F, {axisPhi});
		histos.add("D02Kp_DCA_X", "DCA_X", kTH1F, {axisDCAX});
		histos.add("D02Kp_DCA_Y", "DCA_Y", kTH1F, {axisDCAY});
		histos.add("MFT_Kaon_pT", "MFT_Kaon_pT", kTH1F, {axisPt});
		histos.add("MFT_bg_pT", "MFT_bg_pT", kTH1F, {axisPt});
		histos.add("PredictZ", "PredictZ", kTH1F, {axisZ});
		histos.add("PredictZ_bg", "PredictZ_bg", kTH1F, {axisZ});
		histos.add("PredictX", "PredictX", kTH1F, {axisZ});
		histos.add("PredictX_bg", "PredictX_bg", kTH1F, {axisZ});
		histos.add("PredictY", "PredictY", kTH1F, {axisZ});
		histos.add("PredictY_bg", "PredictY_bg", kTH1F, {axisZ});
		histos.add("Cluster_K", "Cluster_K", kTH1F, {axisClu});
		histos.add("Cluster_bg", "Cluster_bg", kTH1F, {axisClu});
		histos.add("predictPDG", "predictPDG", kTH1F, {axisPDG});
		histos.add("D02k_etaphi", "D02k_etaphi", kTH2F, {axisPhi, axisEta});
		//histos.add("PredictZandPDG", "PredictZandPDG", kTH2F, {axisZ, axisPDGv2});
		histos.add("D0Decay_muK", "D0Decay_muK", kTH1F, {axisZ});
		histos.add("PredictZ_muK", "PredictZ_muK", kTH1F, {axisZ});
		histos.add("Diff_preZ_muK", "Diff_preZ_muK", kTH1F, {{2010, -10.5, 10.5, "z(cm)"}});
		histos.add("PCAZ_muK", "PCAZ_muK", kTH1F, {axisZ});
		histos.add("PCA_Quality", "PCA_Quality", kTH1F, {axisPro});
		histos.add("MFTParticleCounter", "MFTParticleCounter", kTH1F, {axisCounter});
		histos.add("MFTMuonCounter", "MFTMuonCounter", kTH1F, {axisCounter});
		histos.add("MFT_muon_pT", "MFT_muon_pT", kTH1F, {axisPt});
		histos.add("MFT_KaonDCAXY", "MFT_KaonDCAXY", kTH2F, {axisDCAX, axisDCAY});
		histos.add("MFT_PionDCAXY", "MFT_PionDCAXY", kTH2F, {axisDCAX, axisDCAY});
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
		int g = 0;
		if(collision.has_mcCollision()){
			//Informatin of MC collision
			auto mcCol_z = collision.mcCollision().posZ();
			histos.fill(HIST("CollisionPos"), mcCol_z);

			//Information of FwdTracks
			for(auto& fwdtrack : fwdtracks){
				histos.fill(HIST("TrackType"), fwdtrack.trackType());
				auto fwdID = fwdtrack.globalIndex();

				if(!fwdtrack.has_mcParticle()) continue;
				auto mcParticle_fwd = fwdtrack.mcParticle();
				histos.fill(HIST("FwdParticleCounter"), 0.5);
				auto fwdTrackPDG = mcParticle_fwd.pdgCode();

				if(fabs(fwdTrackPDG)==13) histos.fill(HIST("FwdMuonCounter"), 0.5);

				//Informaion of D0 meson
				float mudcaX, mudcaY, mumftX, mumftY, mumftZ, kdcaX, kdcaY, kmftX, kmftY, kmftZ;
				float mu_phi, mu_eta, k_phi, k_eta;
				float secver_z;
				if(fwdtrack.trackType()==fwdTrackType){
					auto chi2GMT = fwdtrack.chi2MatchMCHMFT();
					auto GMT_pT = mcParticle_fwd.pt();
					histos.fill(HIST("GMT_pT"), GMT_pT);

					if(chi2GMT==-1) continue;
					if(fabs(fwdTrackPDG)==13){
						auto a_ID = mcParticle_fwd.globalIndex();

						if(mcParticle_fwd.has_mothers()){
							auto mcMom = MCparticle.rawIteratorAt(mcParticle_fwd.mothersIds()[0]);
							auto mcMomPDG = mcMom.pdgCode();
							auto Daughters = mcMom.daughters_as<aod::McParticles>();
							
							if(fabs(mcMomPDG) < 38) histos.fill(HIST("muonMomPDG"), 0.0, 1);
							if(fabs(mcMomPDG) == 111) histos.fill(HIST("muonMomPDG"), 1.0, 1);
							if(fabs(mcMomPDG) == 113 ) histos.fill(HIST("muonMomPDG"), 3.0, 1);
							if(fabs(mcMomPDG) == 130 ) histos.fill(HIST("muonMomPDG"), 4.0, 1);
							if(fabs(mcMomPDG) == 221 ) histos.fill(HIST("muonMomPDG"), 5.0, 1);
							if(fabs(mcMomPDG) == 223 ) histos.fill(HIST("muonMomPDG"), 6.0, 1);
							if(fabs(mcMomPDG) == 310 ) histos.fill(HIST("muonMomPDG"), 7.0, 1);
							if(fabs(mcMomPDG) == 313 ) histos.fill(HIST("muonMomPDG"), 8.0, 1);
							if(fabs(mcMomPDG) == 321 ) histos.fill(HIST("muonMomPDG"), 9.0, 1);
							if(fabs(mcMomPDG) == 323 ) histos.fill(HIST("muonMomPDG"), 10.0, 1);
							if(fabs(mcMomPDG) == 333 ) histos.fill(HIST("muonMomPDG"), 12.0, 1);
							if(fabs(mcMomPDG) == 411 ) histos.fill(HIST("muonMomPDG"), 13.0, 1);
							if(fabs(mcMomPDG) == 413 ) histos.fill(HIST("muonMomPDG"), 14.0, 1);
							if(fabs(mcMomPDG) == 431 ) histos.fill(HIST("muonMomPDG"), 16.0, 1);
							if(fabs(mcMomPDG) > 499 && fabs(mcMomPDG) < 600) histos.fill(HIST("muonMomPDG"), 17.0, 1);
							if(fabs(mcMomPDG) > 999 ) histos.fill(HIST("muonMomPDG"), 18.0, 1);
							

							if(fabs(mcMomPDG)==211){
								histos.fill(HIST("muonMomPDG"), 2.0, 1);
							}

							if(fabs(mcMomPDG)==421){
								//auto momID = mcMom.globalIndex();
								histos.fill(HIST("muonMomPDG"), 15.0, 1);
								histos.fill(HIST("D0DecayX"), mcParticle_fwd.vx());
								histos.fill(HIST("D0DecayY"), mcParticle_fwd.vy());
								histos.fill(HIST("D0DecayZ"), mcParticle_fwd.vz());
								histos.fill(HIST("D0_pT"), mcMom.pt());
								histos.fill(HIST("Distance_D0_z"), fabs(mcCol_z-mcParticle_fwd.vz()));

								int daughter_count = 0;
								float pcaCan = 1000;
								float pcaZ, pcaX, pcaY, nClu, mftpT;
								int idPDG, k_mom, mu_mom;
								for(auto Daughter : Daughters){
									if(fabs(Daughter.pdgCode())==13){//muon
										auto mu_ID = Daughter.globalIndex();
										if(mu_ID==a_ID){
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
											secver_z = Daughter.vz();
											histos.fill(HIST("D02mu_eta"), Daughter.eta());
											histos.fill(HIST("D02mu_phi"), mu_phi);//Daughter.phi());
											histos.fill(HIST("D02mu_DCA_X"), mudcaX);
											histos.fill(HIST("D02mu_DCA_Y"), mudcaY);
											histos.fill(HIST("D02mu_DCAT"), sqrt(pow(mudcaX, 2.0)*pow(mudcaY, 2.0)));
											mu_mom = MCparticle.rawIteratorAt(Daughter.mothersIds()[0]).globalIndex();
											daughter_count++;
										}
									}
									if(fabs(Daughter.pdgCode())==321){//kaon
										auto k_ID = Daughter.globalIndex();
										for(auto& mfttrack : mfttracks){
											if(!mfttrack.has_collision()) continue;
											if(!mfttrack.has_mcParticle()) continue;
											auto mcParticle_mft = mfttrack.mcParticle();
											auto mftID = mcParticle_mft.globalIndex();
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
												histos.fill(HIST("D02Kp_DCA_X"), kdcaX);
												histos.fill(HIST("D02Kp_DCA_Y"), kdcaY);
												histos.fill(HIST("D02k_etaphi"), k_phi, k_eta);
												k_mom = MCparticle.rawIteratorAt(Daughter.mothersIds()[0]).globalIndex();
												daughter_count++;
											}
										}
									}
								}
								if(daughter_count==2){
									g++;
									histos.fill(HIST("counter1"), 0.5);
									histos.fill(HIST("D0Decay_muK"), secver_z);
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
									auto pre_muz = mumftZ + s*Naz;
									auto pre_kz = kmftZ + t*Ncz;
									histos.fill(HIST("PredictZ_muK"), pre_muz);
									histos.fill(HIST("Diff_preZ_muK"), pre_muz-pre_kz);
									histos.fill(HIST("PCAZ_muK"), (pre_muz+pre_kz)/2);
									
									//Scan all mfttrack for calculate the PCA
									for(auto& mfttrack : mfttracks){
										if(!mfttrack.has_collision() || !mfttrack.has_mcParticle()) continue;
										if(mfttrack.globalIndex()==fwdID) continue;
										
										auto mcParticle_mft = mfttrack.mcParticle();
										//auto mftID = mcParticle_mft.globalIndex();
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

										//問題ありの可能性大
										if(((pre_kz+pre_muz)/2)>zaxisMaxCut || ((pre_kz+pre_muz)/2)<zaxisMinCut) continue;
										//if(fabs(pre_kz+pre_muz)/2<0.05 || (pre_kz+pre_muz/2)==0) continue;
										float dist_2track = sqrt(pow(pre_kx-pre_mux, 2)+pow(pre_ky-pre_muy, 2)+pow(pre_kz-pre_muz, 2));
										if(pcaCan>dist_2track){
											pcaCan = dist_2track;
											idPDG = mcParticle_mft.pdgCode();
											nClu = mfttrack.nClusters();
											mftpT = mfttrack.p();
											pcaZ = (pre_muz+pre_kz)/2;
											pcaX = (pre_mux+pre_kx)/2;
											pcaY = (pre_muy+pre_kx)/2;
										}else if(pcaCan==dist_2track && fabs(mcParticle_mft.pdgCode())==321){
											idPDG = mcParticle_mft.pdgCode();
											pcaZ = (pre_muz+pre_kz)/2;
											pcaX = (pre_mux+pre_kx)/2;
											pcaY = (pre_muy+pre_kx)/2;
											mftpT = mfttrack.p();
											nClu = mfttrack.nClusters();
										}
									}
									if(pcaZ<=0.005 && pcaZ>=-0.005) continue;

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
									if(fabs(idPDG) == 321 ){
										histos.fill(HIST("predictPDG"), 9.0, 1);
										histos.fill(HIST("PredictZ"), pcaZ);
										histos.fill(HIST("PredictX"), pcaX);
										histos.fill(HIST("PredictY"), pcaY);
										histos.fill(HIST("Cluster_K"), nClu);
										histos.fill(HIST("MFT_Kaon_pT"), mftpT);
									}else{
										histos.fill(HIST("PredictZ_bg"), pcaZ);
										histos.fill(HIST("PredictX_bg"), pcaX);
										histos.fill(HIST("PredictY_bg"), pcaY);
										histos.fill(HIST("Cluster_bg"), nClu);
										histos.fill(HIST("MFT_bg_pT"), mftpT);
									}
									if(fabs(idPDG) == 323 ) histos.fill(HIST("predictPDG"), 10.0, 1);
									if(fabs(idPDG) == 333 ) histos.fill(HIST("predictPDG"), 12.0, 1);
									if(fabs(idPDG) == 411 ) histos.fill(HIST("predictPDG"), 13.0, 1);
									if(fabs(idPDG) == 413 ) histos.fill(HIST("predictPDG"), 14.0, 1);
									if(fabs(idPDG) == 431 ) histos.fill(HIST("predictPDG"), 16.0, 1);
									if(fabs(idPDG) > 499 && fabs(idPDG) < 600) histos.fill(HIST("predictPDG"), 17.0, 1);
									if(fabs(idPDG) > 999 ) histos.fill(HIST("predictPDG"), 18.0, 1);

									if(k_mom==mu_mom){
										histos.fill(HIST("counter2"), 0.5);

									}
								}
							}
							if(fabs(mcMomPDG)==321){
								histos.fill(HIST("muonMomPDG"), 11.0, 1);
							}
						}
					}
				}
			}
			//Information of MFTTracks
			for(auto& mfttrack : mfttracks){
				if(!mfttrack.has_mcParticle()) continue;
				auto mcParticle_mft = mfttrack.mcParticle();
				histos.fill(HIST("MFTParticleCounter"), 0.5);
				if(fabs(mcParticle_mft.pdgCode())==13){
					histos.fill(HIST("MFTMuonCounter"), 0.5);
					auto mft_mupT = mcParticle_mft.pt();
					histos.fill(HIST("MFT_muon_pT"), mft_mupT);
				}
				//Hadron DCA
				Double_t verZ = collision.posZ();
				double mftchi2 = mfttrack.chi2();
				SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
				vector<double> mftv1;
				SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
				o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
				mftpars1.propagateToZlinear(verZ);

				if(fabs(mcParticle_mft.pdgCode())==211) histos.fill(HIST("MFT_PionDCAXY"), mftpars1.getX() - collision.posX(), mftpars1.getY() - collision.posY());
				if(fabs(mcParticle_mft.pdgCode())==321) histos.fill(HIST("MFT_KaonDCAXY"), mftpars1.getX() - collision.posX(), mftpars1.getY() - collision.posY());
			}
		}
	}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyAnalysisTask>(cfgc)};
}
