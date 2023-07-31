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

struct CharmInfo{ // MC dateset LHC22f5b
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axistrackType{6, -0.5, 5.5, ""};
    const AxisSpec axisPDG{21, -0.5, 20.5, ""};
    const AxisSpec axisDCAT{100010, -0.0005, 10.0005, "DCAT"};
    const AxisSpec axisPCA{1100010, -60.0005, 50.0005, "z(cm)"};
    const AxisSpec axisCounter{1, 0.5, 1.5, ""};

    histos.add("TrackType", "TrackType", kTH1F, {axistrackType});
    histos.add("MuonMomPDG_Prompt", "MuonMomPDG_Prompt", kTH1F, {axisPDG});
    histos.add("DCAT_all", "DCAT_all", kTH1F, {axisDCAT});
    histos.add("DCAT_Charm", "DCAT_Charm", kTH1F, {axisDCAT});
    histos.add("DCAT_Beauty", "DCAT_Beauty", kTH1F, {axisDCAT});
    histos.add("DCAT_Strange", "DCAT_Strange", kTH1F, {axisDCAT});
    histos.add("DCAT_pi_prompt", "DCAT_pi_prompt", kTH1F, {axisDCAT});
    histos.add("EstPCAZ_all", "EstPCAZ_all", kTH1F, {axisPCA});
    histos.add("pairPDG", "pairPDG", kTH1F, {axisPDG});
    histos.add("EstPCAZ_Charm", "EstPCAZ_Charm", kTH1F, {axisPCA});
    histos.add("EstPCAZ_Beauty", "EstPCAZ_Beauty", kTH1F, {axisPCA});
    histos.add("EstPCAZ_Strange", "EstPCAZ_Strange", kTH1F, {axisPCA});
    histos.add("EstPCAZ_pi", "EstPCAZ_pi", kTH1F, {axisPCA});
    histos.add("MuonMomPDG_Second", "MuonMomPDG_Second", kTH1F, {axisPDG});
    histos.add("DCAT_sec", "DCAT_sec", kTH1F, {axisDCAT});
    histos.add("DCAT_pi_sec", "DCAT_pi_sec", kTH1F, {axisDCAT});
    histos.add("D0_Distance", "D0_Distance", kTH1F, {axisPCA});
    histos.add("Counter", "Counter", kTH1F, {axisCounter});

    auto mumompdg = histos.get<TH1>(HIST("MuonMomPDG_Prompt"));
    auto* x1 = mumompdg->GetXaxis();
    x1->SetBinLabel(1, "Primary"); // 1-37
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
    x1->SetBinLabel(18, "B0"); // 511
    x1->SetBinLabel(19, "B+"); // 521
    x1->SetBinLabel(20, "B0_s"); // 531
    x1->SetBinLabel(21, "Baryon"); // 1000 -

    auto mumompdg_s = histos.get<TH1>(HIST("MuonMomPDG_Second"));
    auto* x2 = mumompdg_s->GetXaxis();
    x2->SetBinLabel(1, "Primary"); // 1-37
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
    x2->SetBinLabel(18, "B0"); // 511
    x2->SetBinLabel(19, "B+"); // 521
    x2->SetBinLabel(20, "B0_s"); // 531
    x2->SetBinLabel(21, "Baryon"); // 1000 -

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
							 aod::McCollisions const&)
  {
    if(collision.has_mcCollision()){
      for(auto& fwdtrack : fwdtracks){
        if(fwdtrack.has_collision() && fwdtrack.has_mcParticle()){
          histos.fill(HIST("TrackType"), fwdtrack.trackType());
          auto mcParticle_fwd = fwdtrack.mcParticle();
          //if(fwdtrack.pt()<=3) continue;
          if(fabs(mcParticle_fwd.pdgCode())==13 && mcParticle_fwd.has_mothers()){
            auto mcMoms = mcParticle_fwd.mothers_as<aod::McParticles>();
            auto mcMom_b = mcMoms.back();
            auto mcMom_f = mcMoms.front();
            /*if(mcMom_b.producedByGenerator()==false){
              LOGF(info, "FrontPDG & GlobalIndex: %d : %d, BackPDG & GlobalIndex: %d : %d", mcMom_f.pdgCode(),mcMom_f.globalIndex(),mcMom_b.pdgCode(),mcMom_b.globalIndex());
            }*/
            if(mcParticle_fwd.mcCollisionId()==collision.mcCollision().globalIndex()){
              auto Col_x = collision.posX();
              auto Col_y = collision.posY();
              auto Col_z = collision.posZ();
              auto mcCol_z = collision.mcCollision().posZ();
              int mucolid, mcmuid, mu_momId;
              float mudcaX, mudcaY, mumftX, mumftY, mumftZ;
              //Muon Mother PDG(Prompt)
              if(mcMom_b.globalIndex()==mcMom_f.globalIndex() && mcMom_b.producedByGenerator()==true){
                if(fabs(mcMom_b.pdgCode())<38) histos.fill(HIST("MuonMomPDG_Prompt"), 0);
                if(fabs(mcMom_b.pdgCode())==111) histos.fill(HIST("MuonMomPDG_Prompt"), 1);
                if(fabs(mcMom_b.pdgCode())==211) histos.fill(HIST("MuonMomPDG_Prompt"), 2);
                if(fabs(mcMom_b.pdgCode())==113) histos.fill(HIST("MuonMomPDG_Prompt"), 3);
                if(fabs(mcMom_b.pdgCode())==130) histos.fill(HIST("MuonMomPDG_Prompt"), 4);
                if(fabs(mcMom_b.pdgCode())==221) histos.fill(HIST("MuonMomPDG_Prompt"), 5);
                if(fabs(mcMom_b.pdgCode())==223) histos.fill(HIST("MuonMomPDG_Prompt"), 6);
                if(fabs(mcMom_b.pdgCode())==310) histos.fill(HIST("MuonMomPDG_Prompt"), 7);
                if(fabs(mcMom_b.pdgCode())==313) histos.fill(HIST("MuonMomPDG_Prompt"), 8);
                if(fabs(mcMom_b.pdgCode())==321) histos.fill(HIST("MuonMomPDG_Prompt"), 9);
                if(fabs(mcMom_b.pdgCode())==323) histos.fill(HIST("MuonMomPDG_Prompt"), 10);
                if(fabs(mcMom_b.pdgCode())==331) histos.fill(HIST("MuonMomPDG_Prompt"), 11);
                if(fabs(mcMom_b.pdgCode())==333) histos.fill(HIST("MuonMomPDG_Prompt"), 12);
                if(fabs(mcMom_b.pdgCode())==411) histos.fill(HIST("MuonMomPDG_Prompt"), 13);
                if(fabs(mcMom_b.pdgCode())==413) histos.fill(HIST("MuonMomPDG_Prompt"), 14);
                if(fabs(mcMom_b.pdgCode())==421) histos.fill(HIST("MuonMomPDG_Prompt"), 15);
                if(fabs(mcMom_b.pdgCode())==431) histos.fill(HIST("MuonMomPDG_Prompt"), 16);
                if(fabs(mcMom_b.pdgCode())==511) histos.fill(HIST("MuonMomPDG_Prompt"), 17);
                if(fabs(mcMom_b.pdgCode())==521) histos.fill(HIST("MuonMomPDG_Prompt"), 18);
                if(fabs(mcMom_b.pdgCode())==531) histos.fill(HIST("MuonMomPDG_Prompt"), 19);
                if(fabs(mcMom_b.pdgCode())>1000) histos.fill(HIST("MuonMomPDG_Prompt"), 20);

                if(fwdtrack.trackType()==0 || fwdtrack.trackType()==2){ //Global Muon Track (MFT-MCH-MID) or Global Track (MFT-MCH)
                  double fwdchi2 = fwdtrack.chi2();
                  SMatrix5 fwdpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
                  vector<double> fwdv1;
                  SMatrix55 fwdcovs(fwdv1.begin(), fwdv1.end());
                  o2::track::TrackParCovFwd fwdpars1{fwdtrack.z(), fwdpars, fwdcovs, fwdchi2};
                  fwdpars1.propagateToZlinear(Col_z);

                  mudcaX = fwdpars1.getX();
                  mudcaY = fwdpars1.getY();
                  mumftX = fwdtrack.x();
                  mumftY = fwdtrack.y();
                  mumftZ = fwdtrack.z();
                  mucolid = fwdtrack.collisionId();
                  mcmuid = fwdtrack.mcParticleId();
                  mu_momId = mcMom_b.globalIndex();
                  auto mftmchId = fwdtrack.matchMFTTrackId();

                  //Calcurate PCA
                  float candcaX, candcaY, canmftX, canmftY, canmftZ, estpcax, estpcay, estpcaz;
                  int mcMom_b_canId, estMomId, estpairPDG;
                  float dist_2track = 10000;
                  for(auto& mfttrack : mfttracks){
                    if(mfttrack.has_collision() && mfttrack.has_mcParticle()){
                      if(mfttrack.collisionId()==mucolid && mfttrack.globalIndex()!=mftmchId){
                        auto mcParticle_mft = mfttrack.mcParticle();
                        double mftchi2 = mfttrack.chi2();
                        SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                        vector<double> mftv1;
                        SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                        o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
                        mftpars1.propagateToZlinear(Col_z);

                        candcaX = mftpars1.getX();
                        candcaY = mftpars1.getY();
                        canmftX = mfttrack.x();
                        canmftY = mfttrack.y();
                        canmftZ = mfttrack.z();
                        if(mcParticle_mft.has_mothers()){
                          auto mcMoms_can = mcParticle_mft.mothers_as<aod::McParticles>();
                          mcMom_b_canId = mcMoms_can.back().globalIndex();
                        }else{
                          mcMom_b_canId = -10;
                        }

                        auto Nax = mudcaX - mumftX;
                        auto Nay = mudcaY - mumftY;
                        auto Naz = Col_z - mumftZ;
                        auto Ncx = candcaX - canmftX;
                        auto Ncy = candcaY - canmftY;
                        auto Ncz = Col_z - canmftZ;
                        auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
                        auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
                        auto A3 = (mumftX-canmftX)*Nax + (mumftY-canmftY)*Nay + (mumftZ-canmftZ)*Naz;
                        auto B1 = A2;
                        auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
                        auto B3 = (canmftX-mumftX)*Ncx + (canmftY-mumftY)*Ncy + (canmftZ-mumftZ)*Ncz;
                        auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
                        auto s = -((A2*t+A3)/A1);
                        float pre_mux = mumftX + s*Nax;
                        float pre_muy = mumftY + s*Nay;
                        float pre_muz = mumftZ + s*Naz;
                        float pre_canx = canmftX + t*Ncx;
                        float pre_cany = canmftY + t*Ncy;
                        float pre_canz = canmftZ + t*Ncz;
                        float pcax_org = ((pre_mux+pre_canx)/2) - Col_x;
                        float pcay_org = ((pre_muy+pre_cany)/2) - Col_y;
                        float pcaz_org = ((pre_muz+pre_canz)/2) - Col_z;
                        auto r_xyz = sqrt(pow(pre_mux-pre_canx,2)+pow(pre_muy-pre_cany,2)+pow(pre_muz-pre_canz,2));

                        if(dist_2track>=r_xyz){
                          dist_2track = r_xyz;
                          estpcax = pcax_org;
                          estpcay = pcay_org;
                          estpcaz = pcaz_org;
                          estMomId = mcMom_b_canId;
                          estpairPDG = mcParticle_mft.pdgCode();
                        }
                      }
                    }
                  }
                  if(mu_momId==estMomId){// Correct pairs
                    histos.fill(HIST("Counter"), 1);
                  }

                  //DCAT Infomation
                  auto DCAT = sqrt(pow(mudcaX-Col_x,2)+pow(mudcaY-Col_y,2));
                  histos.fill(HIST("DCAT_all"), DCAT);
                  histos.fill(HIST("EstPCAZ_all"), estpcaz);
                  if(fabs(estpairPDG)<38) histos.fill(HIST("pairPDG"), 0);
                  if(fabs(estpairPDG)==111) histos.fill(HIST("pairPDG"), 1);
                  if(fabs(estpairPDG)==211) histos.fill(HIST("pairPDG"), 2);
                  if(fabs(estpairPDG)==113) histos.fill(HIST("pairPDG"), 3);
                  if(fabs(estpairPDG)==130) histos.fill(HIST("pairPDG"), 4);
                  if(fabs(estpairPDG)==221) histos.fill(HIST("pairPDG"), 5);
                  if(fabs(estpairPDG)==223) histos.fill(HIST("pairPDG"), 6);
                  if(fabs(estpairPDG)==310) histos.fill(HIST("pairPDG"), 7);
                  if(fabs(estpairPDG)==313) histos.fill(HIST("pairPDG"), 8);
                  if(fabs(estpairPDG)==321) histos.fill(HIST("pairPDG"), 9);
                  if(fabs(estpairPDG)==323) histos.fill(HIST("pairPDG"), 10);
                  if(fabs(estpairPDG)==331) histos.fill(HIST("pairPDG"), 11);
                  if(fabs(estpairPDG)==333) histos.fill(HIST("pairPDG"), 12);
                  if(fabs(estpairPDG)==411) histos.fill(HIST("pairPDG"), 13);
                  if(fabs(estpairPDG)==413) histos.fill(HIST("pairPDG"), 14);
                  if(fabs(estpairPDG)==421) histos.fill(HIST("pairPDG"), 15);
                  if(fabs(estpairPDG)==431) histos.fill(HIST("pairPDG"), 16);
                  if(fabs(estpairPDG)==511) histos.fill(HIST("pairPDG"), 17);
                  if(fabs(estpairPDG)==521) histos.fill(HIST("pairPDG"), 18);
                  if(fabs(estpairPDG)==531) histos.fill(HIST("pairPDG"), 19);
                  if(fabs(estpairPDG)>1000) histos.fill(HIST("pairPDG"), 20);


                  if(fabs(mcMom_b.pdgCode())==421){
                    auto muvz = mcParticle_fwd.vz();
                    auto Daughters = mcMom_b.daughters_as<aod::McParticles>();
                    int coun = 0;
                    for(auto& Daughter : Daughters){
                      if(fabs(Daughter.pdgCode())==13 || fabs(Daughter.pdgCode())==321){
                        coun++;
                      }
                    }
                    if(coun==2){
                      histos.fill(HIST("D0_Distance"), muvz-mcCol_z);
                      coun = 0;
                    }
                  }
                  if(fabs(mcMom_b.pdgCode())>=400 && fabs(mcMom_b.pdgCode())<500){
                    histos.fill(HIST("DCAT_Charm"), DCAT);
                    histos.fill(HIST("EstPCAZ_Charm"), estpcaz);
                  }
                  if(fabs(mcMom_b.pdgCode())>=500 && fabs(mcMom_b.pdgCode())<600){
                    histos.fill(HIST("DCAT_Beauty"), DCAT);
                    histos.fill(HIST("EstPCAZ_Beauty"), estpcaz);
                  }
                  if(fabs(mcMom_b.pdgCode())>=300 & fabs(mcMom_b.pdgCode())<400){
                    histos.fill(HIST("DCAT_Strange"), DCAT);
                    histos.fill(HIST("EstPCAZ_Strange"), estpcaz);
                  }
                  if(fabs(mcMom_b.pdgCode())==211){
                    histos.fill(HIST("DCAT_pi_prompt"), DCAT);
                    histos.fill(HIST("EstPCAZ_pi"), estpcaz);
                  }
                }
              }
              if(mcMom_b.producedByGenerator()==false){
                if(fabs(mcMom_b.pdgCode())<38) histos.fill(HIST("MuonMomPDG_Second"), 0);
                if(fabs(mcMom_b.pdgCode())==111) histos.fill(HIST("MuonMomPDG_Second"), 1);
                if(fabs(mcMom_b.pdgCode())==211) histos.fill(HIST("MuonMomPDG_Second"), 2);
                if(fabs(mcMom_b.pdgCode())==113) histos.fill(HIST("MuonMomPDG_Second"), 3);
                if(fabs(mcMom_b.pdgCode())==130) histos.fill(HIST("MuonMomPDG_Second"), 4);
                if(fabs(mcMom_b.pdgCode())==221) histos.fill(HIST("MuonMomPDG_Second"), 5);
                if(fabs(mcMom_b.pdgCode())==223) histos.fill(HIST("MuonMomPDG_Second"), 6);
                if(fabs(mcMom_b.pdgCode())==310) histos.fill(HIST("MuonMomPDG_Second"), 7);
                if(fabs(mcMom_b.pdgCode())==313) histos.fill(HIST("MuonMomPDG_Second"), 8);
                if(fabs(mcMom_b.pdgCode())==321) histos.fill(HIST("MuonMomPDG_Second"), 9);
                if(fabs(mcMom_b.pdgCode())==323) histos.fill(HIST("MuonMomPDG_Second"), 10);
                if(fabs(mcMom_b.pdgCode())==331) histos.fill(HIST("MuonMomPDG_Second"), 11);
                if(fabs(mcMom_b.pdgCode())==333) histos.fill(HIST("MuonMomPDG_Second"), 12);
                if(fabs(mcMom_b.pdgCode())==411) histos.fill(HIST("MuonMomPDG_Second"), 13);
                if(fabs(mcMom_b.pdgCode())==413) histos.fill(HIST("MuonMomPDG_Second"), 14);
                if(fabs(mcMom_b.pdgCode())==421) histos.fill(HIST("MuonMomPDG_Second"), 15);
                if(fabs(mcMom_b.pdgCode())==431) histos.fill(HIST("MuonMomPDG_Second"), 16);
                if(fabs(mcMom_b.pdgCode())==511) histos.fill(HIST("MuonMomPDG_Second"), 17);
                if(fabs(mcMom_b.pdgCode())==521) histos.fill(HIST("MuonMomPDG_Second"), 18);
                if(fabs(mcMom_b.pdgCode())==531) histos.fill(HIST("MuonMomPDG_Second"), 19);
                if(fabs(mcMom_b.pdgCode())>1000) histos.fill(HIST("MuonMomPDG_Second"), 20);

                if(fwdtrack.trackType()==0 || fwdtrack.trackType()==2){ //Global Muon Track (MFT-MCH-MID) or Global Track (MFT-MCH)
                  double fwdchi2 = fwdtrack.chi2();
                  SMatrix5 fwdpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
                  vector<double> fwdv1;
                  SMatrix55 fwdcovs(fwdv1.begin(), fwdv1.end());
                  o2::track::TrackParCovFwd fwdpars1{fwdtrack.z(), fwdpars, fwdcovs, fwdchi2};
                  fwdpars1.propagateToZlinear(Col_z);

                  mudcaX = fwdpars1.getX() - Col_x;
                  mudcaY = fwdpars1.getY() - Col_y;
                  mumftX = fwdtrack.x();
                  mumftY = fwdtrack.y();
                  mumftZ = fwdtrack.z();
                  mucolid = fwdtrack.collisionId();
                  mcmuid = fwdtrack.mcParticleId();

                  //DCAT Infomation
                  auto DCAT = sqrt(pow(mudcaX,2)+pow(mudcaY,2));
                  histos.fill(HIST("DCAT_all"), DCAT);
                  if(fabs(mcMom_b.pdgCode())==211){
                    histos.fill(HIST("DCAT_pi_sec"), DCAT);
                  }else{
                    histos.fill(HIST("DCAT_sec"), DCAT);
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

struct MFTMacthingInfo{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&){
    const AxisSpec axischi2{501, -0.5, 500.5, ""};
    const AxisSpec axispT{3000010, -0.00005, 30.00005, "pT"};

    histos.add("MFTMCH_chi2", "MFTMCH_chi2", kTH1F, {axischi2});
    histos.add("TrueMatch_pT", "TrueMatch_pT", kTH1F, {axispT});
    histos.add("FakeMatch_pT", "FakeMatch_pT", kTH1F, {axispT});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
							 soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA> const& fwdtracks, 
							 soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
							 aod::McParticles const& MCparticle,
							 aod::McCollisions const&)
  {
    if(collision.has_mcCollision()){
      for(auto& fwdtrack : fwdtracks){
        if(fwdtrack.has_collision() && fwdtrack.has_mcParticle()){
          if(fwdtrack.trackType()==0 || fwdtrack.trackType()==2){
            auto fwd_colId = fwdtrack.collisionId();
            auto mftmch_chi2 = fwdtrack.chi2MatchMCHMFT();
            auto mftmch_id = fwdtrack.matchMFTTrackId();
            auto fwd_muid = fwdtrack.mcParticle().globalIndex();
            auto fwd_pdg = fwdtrack.mcParticle().pdgCode();
            auto fwd_pt = fwdtrack.pt();
            histos.fill(HIST("MFTMCH_chi2"), mftmch_chi2);
            for(auto& mfttrack : mfttracks){
              if(mfttrack.has_collision() && mfttrack.has_mcParticle()){
                if(mfttrack.collisionId()==fwd_colId && mfttrack.globalIndex()==mftmch_id){
                  auto mft_muid = mfttrack.mcParticle().globalIndex();
                  if(fabs(fwd_pdg)==13 && fwd_muid==mft_muid){ //true matching
                    histos.fill(HIST("TrueMatch_pT"), fwd_pt);
                  }
                  if(fabs(fwd_pdg)==13 && fwd_muid!=mft_muid){ //fake matching
                    histos.fill(HIST("FakeMatch_pT"), fwd_pt);
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
    adaptAnalysisTask<CharmInfo>(cfgc),
    adaptAnalysisTask<MFTMacthingInfo>(cfgc)};
}