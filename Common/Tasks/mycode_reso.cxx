/// \author Koki Soeda
/// \since 13/09/2023

#include <iostream>
#include <fstream>
#include <array>
#include <cmath>
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "CommonConstants/GeomConstants.h"

#include "CCDB/BasicCCDBManager.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "Math/MatrixFunctions.h"
#include "MathUtils/Utils.h"
#include "TMath.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "MFTTracking/Tracker.h"
#include "Framework/ASoAHelpers.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DataFormatsParameters/GRPMagField.h"
#include <math.h>
#include <TLorentzVector.h>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::track;
using namespace evsel;
using o2::field::MagneticField;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

struct ResolutionTask{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int runNumber = -1;
  float Bz = 0;                                         // Magnetic field for MFT
  static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  o2::parameters::GRPMagField* grpmag = nullptr;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  void init(InitContext const&){
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    const AxisSpec axisD{300, 30, 60, "DistanceFromMFT(cm)"};
    const AxisSpec axisDeltaVtx{600, 0, 0.3, "#Delta Vertex (cm)"};
    const AxisSpec axisDeltaX{400, -1, 1, "#Delta(x_{mc}-x_{det}) (cm)"};
    const AxisSpec axisDeltaY{400, -1, 1, "#Delta(y_{mc}-y_{det}) (cm)"};
    const AxisSpec axisDistance{5000, 0, 5, "d(cm)"};
    const AxisSpec axisDeltaMFT{2500, 0, 5, "d(cm)"};
    const AxisSpec axisDeltaMFTAll{2500, 0, 5, "#Delta d(cm)"};
    const AxisSpec axisDeltaEta{501, -0.005, 5.005, "|#Delta #eta|"};
    const AxisSpec axisDeltaPhi{701, -0.005, 7.005, "|#Delta #phi|"};
    const AxisSpec axisP{500, 0, 50, "p(GeV/c)"};
    const AxisSpec axispT{500, 0, 50, "pT(GeV/c)"};
    const AxisSpec axisPID{5, -0.5, 4.5, ""};
    histos.add("Reso_atSecVertex", "Reso_atSecVertex", kTH2F, {axisD, axisDeltaVtx});
    histos.add("DeltaX_atSecVertexZ", "DeltaX_atSecVertexZ", kTH1F, {axisDeltaX});
    histos.add("DeltaY_atSecVertexZ", "DeltaY_atSecVertexZ", kTH1F, {axisDeltaY});
    histos.add("Dist_fromSecVertex", "Dist_fromSecVertex", kTH1F, {axisDistance});
    histos.add("Dist_fromSecVertex_Pair", "Dist_fromSecVertex_Pair", kTH1F, {axisDistance});
    histos.add("Dist_fromSecVertex_BG", "Dist_fromSecVertex_BG", kTH1F, {axisDistance});
    histos.add("AveDist_withOthers_atDecayPosz", "AveDist_withOthers_atDecayPosz", kTH1F, {axisDistance});
    histos.add("MuonDelta_atMFT_Linear", "MuonDelta_atMFT_Linear", kTH2F, {axisP, axisDeltaMFT});
    histos.add("MuonDelta_atMFT_Helix", "MuonDelta_atMFT_Helix", kTH2F, {axisP, axisDeltaMFT});
    histos.add("MFTTrack_p_pdg", "MFTTrack_p_pdg", kTH2F, {axisP, axisPID});
    histos.add("TrackDelta_atMFT_Linear", "TrackDelta_atMFT_Linear", kTH2F, {axisP, axisDeltaMFTAll});
    histos.add("TrackDelta_atMFT_Helix", "TrackDelta_atMFT_Helix", kTH2F, {axisP, axisDeltaMFTAll});
    histos.add("TrackDelta_atMFT_Linear_pT", "TrackDelta_atMFT_Linear_pT", kTH2F, {axispT, axisDeltaMFTAll});
    histos.add("TrackDelta_atMFT_Helix_pT", "TrackDelta_atMFT_Helix_pT", kTH2F, {axispT, axisDeltaMFTAll});
    histos.add("EtaDelta_atMFT", "EtaDelta_atMFT", kTH2F, {axisP, axisDeltaEta});
    histos.add("EtaDelta_atMFT_pT", "EtaDelta_atMFT_pT", kTH2F, {axisP, axisDeltaEta});
    histos.add("PhiDelta_atMFT", "PhiDelta_atMFT", kTH2F, {axisP, axisDeltaPhi});
    histos.add("PhiDelta_atMFT_pT", "PhiDelta_atMFT_pT", kTH2F, {axispT,axisDeltaPhi});

    auto hpdgcode = histos.get<TH2>(HIST("MFTTrack_p_pdg"));
    auto* x1 = hpdgcode->GetYaxis();
    x1->SetBinLabel(1, "Electron"); // 11
    x1->SetBinLabel(2, "Muon"); // 13
    x1->SetBinLabel(3, "Pion±"); // 211
    x1->SetBinLabel(4, "K±"); // 321
  }

  void initCCDB(ExtBCs::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current()
              << " A for run " << bc.runNumber()
              << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    runNumber = bc.runNumber();

    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    Bz = field->getBz(centerMFT);
    LOG(info) << "The field at the center of the MFT is Bz = " << Bz;
  }


  void process(soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions,
               soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA, aod::FwdTracksCov> const& fwdtracks, 
               soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
               aod::McParticles const& MCparticle,
               aod::AmbiguousMFTTracks const& amfttracks,
               aod::McCollisions const&,
               ExtBCs const& bcs)
  {
    initCCDB(bcs.begin());
    const double mftzpos = -46.0; // cm
    for(auto& fwdtrack : fwdtracks){
      if(!fwdtrack.has_mcParticle()) continue;
      if(fwdtrack.trackType()==1 || fwdtrack.trackType()==3 || fwdtrack.trackType()==4) continue;
      if(fwdtrack.eta()<=(-3.6) || fwdtrack.eta()>=(-2.5)) continue;
      //if(fwdtrack.pt()<0.5) continue;

      auto mcParticle_fwd = fwdtrack.mcParticle();
      int mcCol_id = -1;
      double Col_x, Col_y, Col_z;
      double mcCol_x, mcCol_y, mcCol_z;
      bool Has_MCCol = false;
      for(auto& t_col : collisions){
        if(!t_col.has_mcCollision()) continue;
        if(mcParticle_fwd.mcCollisionId()==t_col.mcCollisionId()){
          mcCol_id = t_col.mcCollisionId();
          Col_x = t_col.posX();
          Col_y = t_col.posY();
          Col_z = t_col.posZ();
          mcCol_x = t_col.mcCollision().posX();
          mcCol_y = t_col.mcCollision().posY();
          mcCol_z = t_col.mcCollision().posZ();
          Has_MCCol = true;
          break;
        }
      }
      if(fabs(mcCol_z)>10 || Has_MCCol==false) continue;
      if(fabs(mcParticle_fwd.pdgCode())!=13 || !mcParticle_fwd.has_mothers()) continue;

      int muid = mcParticle_fwd.globalIndex();
      double mu_p = mcParticle_fwd.p();
      auto mumom = mcParticle_fwd.mothers_first_as<aod::McParticles>();
      auto Daughters = mumom.daughters_as<aod::McParticles>();
      double convert_phi = 0;
      if(mcParticle_fwd.phi()>o2::math_utils::pi() && mcParticle_fwd.phi()<2*o2::math_utils::pi()){
        convert_phi = mcParticle_fwd.phi() - 2*(o2::math_utils::pi());
      }else{
        convert_phi = mcParticle_fwd.phi();
      }
      double mu_signedpt = mcParticle_fwd.pdgCode()/(fabs(mcParticle_fwd.pdgCode())*mcParticle_fwd.pt());

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
        PartType = 4; // LF Muon(BG)
      }

      double secverx, secvery, secverz, mueta;
      vector<int> vec_daughter;
      vec_daughter.clear();
      bool HasPair = false;
      for(auto& Daughter : Daughters){
        int d_pdg = fabs(Daughter.pdgCode());
        if(d_pdg==13 && Daughter.globalIndex()==muid){
          secverx = Daughter.vx();
          secvery = Daughter.vy();
          secverz = Daughter.vz();
          mueta = Daughter.eta();
        }
        vec_daughter.push_back(Daughter.globalIndex());
      }
      if(vec_daughter.size()>1) HasPair = true;
      if(!HasPair) continue;

      vector<double> muon_vec;
      SMatrix55 muoncovs(muon_vec.begin(), muon_vec.end());
      SMatrix5 muonpars(secverx, secvery, convert_phi, tan((o2::math_utils::pi()/2)-2*atan(exp(-mueta))), mu_signedpt);
      o2::track::TrackParCovFwd muonpropa(secverz, muonpars, muoncovs, fwdtrack.chi2());
      muonpropa.propagateToZlinear(fwdtrack.z());
      double mumftX_mc = muonpropa.getX();
      double mumftY_mc = muonpropa.getY();
      o2::track::TrackParCovFwd muonpropa_helix(secverz, muonpars, muoncovs, fwdtrack.chi2());
      muonpropa_helix.propagateToZhelix(fwdtrack.z(), Bz);
      double mumftX_mc_helix = muonpropa_helix.getX();
      double mumftY_mc_helix = muonpropa_helix.getY();
      //muonpropa.propagateToZlinear(mcCol_z);
      //double mudcaX_mc = muonpropa.getX();
      //double mudcaY_mc = muonpropa.getY();

      double n_x = (mumftX_mc-secverx)/sqrt(pow(mumftX_mc-secverx,2)+pow(mumftY_mc-secvery,2)+pow(mftzpos-secverz,2));
      double n_y = (mumftY_mc-secvery)/sqrt(pow(mumftX_mc-secverx,2)+pow(mumftY_mc-secvery,2)+pow(mftzpos-secverz,2));
      double n_z = (mftzpos-secverz)/sqrt(pow(mumftX_mc-secverx,2)+pow(mumftY_mc-secvery,2)+pow(mftzpos-secverz,2));
      double mc_d = sqrt(pow(mumftX_mc-secverx,2)+pow(mumftY_mc-secvery,2)+pow(mftzpos-secverz,2));

      vector<vector<double>> xy_atDecayz;
      xy_atDecayz.clear();
      for(auto& mfttrack : mfttracks){
        if(!mfttrack.has_mcParticle()) continue;
        if(mfttrack.mcParticle().mcCollisionId()!=mcCol_id) continue;
        auto mcParticle_mft = mfttrack.mcParticle();
        int mftid = mcParticle_mft.globalIndex();
        int pdg_abs = fabs(mcParticle_mft.pdgCode());
        if(pdg_abs==11){
          histos.fill(HIST("MFTTrack_p_pdg"), mcParticle_mft.p(), 0);
        }else if(pdg_abs==13){
          histos.fill(HIST("MFTTrack_p_pdg"), mcParticle_mft.p(), 1);
        }else if(pdg_abs==211){
          histos.fill(HIST("MFTTrack_p_pdg"), mcParticle_mft.p(), 2);
        }else if(pdg_abs==321){
          histos.fill(HIST("MFTTrack_p_pdg"), mcParticle_mft.p(), 3);
        }
        double convert_phi_mft = 0;
        if(mcParticle_mft.phi()>o2::math_utils::pi() && mcParticle_mft.phi()<2*o2::math_utils::pi()){
          convert_phi_mft = mcParticle_mft.phi() - 2*(o2::math_utils::pi());
        }else{
          convert_phi_mft = mcParticle_mft.phi();
        }
        double mft_signedpt = mcParticle_mft.pdgCode()/(fabs(mcParticle_mft.pdgCode())*mcParticle_mft.pt());

        vector<double> allcov;
        SMatrix55 allcovs(allcov.begin(), allcov.end());
        SMatrix5 allpars(mcParticle_mft.vx(), mcParticle_mft.vy(), convert_phi_mft, tan((o2::math_utils::pi()/2)-2*atan(exp(-mcParticle_mft.eta()))), mft_signedpt);
        o2::track::TrackParCovFwd mftpropa_linear(mcParticle_mft.vz(), allpars, allcovs, mfttrack.chi2());
        o2::track::TrackParCovFwd mftpropa_helix(mcParticle_mft.vz(), allpars, allcovs, mfttrack.chi2());
        mftpropa_linear.propagateToZlinear(mfttrack.z());
        mftpropa_helix.propagateToZhelix(mfttrack.z(), Bz);
        double all_delta_linear = sqrt(pow(mfttrack.x()-mftpropa_linear.getX(),2)+pow(mfttrack.y()-mftpropa_linear.getY(),2));
        double all_delta_helix = sqrt(pow(mfttrack.x()-mftpropa_helix.getX(),2)+pow(mfttrack.y()-mftpropa_helix.getY(),2));
        histos.fill(HIST("TrackDelta_atMFT_Linear"), mcParticle_mft.p(), all_delta_linear);
        histos.fill(HIST("TrackDelta_atMFT_Linear_pT"), mcParticle_mft.pt(), all_delta_linear);
        histos.fill(HIST("TrackDelta_atMFT_Helix"), mcParticle_mft.p(), all_delta_helix);
        histos.fill(HIST("TrackDelta_atMFT_Helix_pT"), mcParticle_mft.pt(), all_delta_helix);
        histos.fill(HIST("EtaDelta_atMFT"), mcParticle_mft.p(), fabs(mcParticle_mft.eta()-mfttrack.eta()));
        histos.fill(HIST("EtaDelta_atMFT_pT"), mcParticle_mft.pt(), fabs(mcParticle_mft.eta()-mfttrack.eta()));
        histos.fill(HIST("PhiDelta_atMFT"), mcParticle_mft.p(), fabs(convert_phi_mft-mfttrack.phi()));
        histos.fill(HIST("PhiDelta_atMFT_pT"), mcParticle_mft.pt(), fabs(convert_phi_mft-mfttrack.phi()));

        if(mftid==muid){
          double mftx = mfttrack.x();
          double mfty = mfttrack.y();
          double mftz = mfttrack.z();
          double mftphi = mfttrack.phi();
          double mfttgl = mfttrack.tgl();
          double mftsignedpt = mfttrack.signed1Pt();
          double mftchi2 = mfttrack.chi2();
          vector<double> covs;
          SMatrix55 mftcovs(covs.begin(), covs.end());
          SMatrix5 mftpars(mftx, mfty, mftphi, mfttgl, mftsignedpt);
          o2::track::TrackParCovFwd mftpropa(mftz, mftpars, mftcovs, mftchi2);
          mftpropa.propagateToZlinear(Col_z);
          double dcax = mftpropa.getX();
          double dcay = mftpropa.getY();
          double unit_x = (mftx-dcax)/sqrt(pow(mftx-dcax,2)+pow(mfty-dcay,2)+pow(mftz-Col_z,2));
          double unit_y = (mfty-dcay)/sqrt(pow(mftx-dcax,2)+pow(mfty-dcay,2)+pow(mftz-Col_z,2));
          double unit_z = (mftz-Col_z)/sqrt(pow(mftx-dcax,2)+pow(mfty-dcay,2)+pow(mftz-Col_z,2));
          double s = (n_x*(secverx-mftx) + n_y*(secvery-mfty) + n_z*(secverz-mftz))/(n_x*unit_x + n_y*unit_y + n_z*unit_z);
          double delta_vtx = sqrt(pow(secverx-(mftx+s*unit_x),2)+pow(secvery-(mfty+s*unit_y),2)+pow(secverz-(mftz+s*unit_z),2));
          histos.fill(HIST("Reso_atSecVertex"), mc_d, delta_vtx);

          double delta_mft_linear = sqrt(pow(mumftX_mc-mftx,2)+pow(mumftY_mc-mfty,2));
          histos.fill(HIST("MuonDelta_atMFT_Linear"), mu_p, delta_mft_linear);
          double delta_mft_helix = sqrt(pow(mumftX_mc_helix-mftx,2)+pow(mumftY_mc_helix-mfty,2));
          histos.fill(HIST("MuonDelta_atMFT_Helix"), mu_p, delta_mft_helix);

          double va = dcax - mftx;
          double vb = dcay - mfty;
          double vc = Col_z - mftz;
          double vx = secverx - mftx;
          double vy = secvery - mfty;
          double vz = secverz - mftz;
          double cross = fabs(vb*vz - vy*vc + vx*vb - va*vy + va*vz - vx*vc);
          double d = cross / sqrt(pow(va,2)+pow(vb,2)+pow(vc,2));
          histos.fill(HIST("Dist_fromSecVertex"), d);
          mftpropa.propagateToZlinear(secverz);
          histos.fill(HIST("DeltaX_atSecVertexZ"), secverx-mftpropa.getX());
          histos.fill(HIST("DeltaY_atSecVertexZ"), secvery-mftpropa.getY());
        }
        if(!vec_daughter.empty()){
          if(find(vec_daughter.begin(), vec_daughter.end(), mcParticle_mft.globalIndex())!=vec_daughter.end()){
            double pairmftx = mfttrack.x();
            double pairmfty = mfttrack.y();
            double pairmftz = mfttrack.z();
            double pairphi = mfttrack.phi();
            double pairtgl = mfttrack.tgl();
            double pairsignedpt = mfttrack.signed1Pt();
            double pairchi2 = mfttrack.chi2();
            vector<double> covpair;
            SMatrix55 mftpaircovs(covpair.begin(), covpair.end());
            SMatrix5 mftpairpars(pairmftx, pairmfty, pairphi, pairtgl, pairsignedpt);
            o2::track::TrackParCovFwd pairpropa(pairmftz, mftpairpars, mftpaircovs, pairchi2);
            pairpropa.propagateToZlinear(Col_z);
            double dcax = pairpropa.getX();
            double dcay = pairpropa.getY();
            double va = dcax - pairmftx;
            double vb = dcay - pairmfty;
            double vc = Col_z - pairmftz;
            double vx = secverx - pairmftx;
            double vy = secvery - pairmfty;
            double vz = secverz - pairmftz;
            double cross = fabs(vb*vz - vy*vc + vx*vb - va*vy + va*vz - vx*vc);
            double d = cross / sqrt(pow(va,2)+pow(vb,2)+pow(vc,2));
            histos.fill(HIST("Dist_fromSecVertex_Pair"), d);
            double pairphi_mc = 0;
            if(mcParticle_mft.phi()>o2::math_utils::pi() && mcParticle_mft.phi()<2*o2::math_utils::pi()){
              pairphi_mc = mcParticle_mft.phi() - 2*(o2::math_utils::pi());
            }else{
              pairphi_mc = mcParticle_mft.phi();
            }
            double pairsignedpt_mc = mcParticle_mft.pdgCode()/(fabs(mcParticle_mft.pdgCode())*mcParticle_mft.pt());
            vector<double> covpair_mc;
            SMatrix55 mftpaircovs_mc(covpair_mc.begin(), covpair_mc.end());
            SMatrix5 mftpairpars_mc(mcParticle_mft.vx(), mcParticle_mft.vy(), pairphi_mc, tan((o2::math_utils::pi()/2)-2*atan(exp(-(mcParticle_mft.eta())))), pairsignedpt_mc);
            o2::track::TrackParCovFwd pairpropa_mc(mcParticle_mft.vz(), mftpairpars_mc, mftpaircovs_mc, pairchi2);
            pairpropa_mc.propagateToZlinear(secverz);
            vector<double> pair_xy = {pairpropa_mc.getX(), pairpropa_mc.getY()};
            xy_atDecayz.push_back(pair_xy);
          }else{
            double bgmftx = mfttrack.x();
            double bgmfty = mfttrack.y();
            double bgmftz = mfttrack.z();
            double bgphi = mfttrack.phi();
            double bgtgl = mfttrack.tgl();
            double bgsignedpt = mfttrack.signed1Pt();
            double bgchi2 = mfttrack.chi2();
            vector<double> covbg;
            SMatrix55 mftbgcovs(covbg.begin(), covbg.end());
            SMatrix5 mftbgpars(bgmftx, bgmfty, bgphi, bgtgl, bgsignedpt);
            o2::track::TrackParCovFwd bgpropa(bgmftz, mftbgpars, mftbgcovs, bgchi2);
            bgpropa.propagateToZlinear(Col_z);
            double dcax = bgpropa.getX();
            double dcay = bgpropa.getY();
            double va = dcax - bgmftx;
            double vb = dcay - bgmfty;
            double vc = Col_z - bgmftz;
            double vx = secverx - bgmftx;
            double vy = secvery - bgmfty;
            double vz = secverz - bgmftz;
            double cross = fabs(vb*vz - vy*vc + vx*vb - va*vy + va*vz - vx*vc);
            double d = cross / sqrt(pow(va,2)+pow(vb,2)+pow(vc,2));
            histos.fill(HIST("Dist_fromSecVertex_BG"), d);

            double bgphi_mc = 0;
            if(mcParticle_mft.phi()>o2::math_utils::pi() && mcParticle_mft.phi()<2*o2::math_utils::pi()){
              bgphi_mc = mcParticle_mft.phi() - 2*(o2::math_utils::pi());
            }else{
              bgphi_mc = mcParticle_mft.phi();
            }
            double bgsignedpt_mc = mcParticle_mft.pdgCode()/(fabs(mcParticle_mft.pdgCode())*mcParticle_mft.pt());
            vector<double> covbg_mc;
            SMatrix55 mftbgcovs_mc(covbg_mc.begin(), covbg_mc.end());
            SMatrix5 mftbgpars_mc(mcParticle_mft.vx(), mcParticle_mft.vy(), bgphi_mc, tan((o2::math_utils::pi()/2)-2*atan(exp(-(mcParticle_mft.eta())))), bgsignedpt_mc);
            o2::track::TrackParCovFwd bgpropa_mc(mcParticle_mft.vz(), mftbgpars_mc, mftbgcovs_mc, bgchi2);
            bgpropa_mc.propagateToZlinear(secverz);
            vector<double> bg_xy = {bgpropa_mc.getX(), bgpropa_mc.getY()};
            xy_atDecayz.push_back(bg_xy);
          }
        }
      }
      if(xy_atDecayz.empty()) continue;
      for(int i=0; i<xy_atDecayz.size(); i++){
        double std_x = xy_atDecayz.at(i).at(0);
        double std_y = xy_atDecayz.at(i).at(1);
        vector<double> dist_all;
        for(int j=0; j<xy_atDecayz.size(); j++){
          if(j==i) continue;
          double diff_x = std_x - xy_atDecayz.at(j).at(0);
          double diff_y = std_y - xy_atDecayz.at(j).at(1);
          double dist_xy = sqrt(pow(diff_x,2)+pow(diff_y,2));
          dist_all.push_back(dist_xy);
        }
        double ave_dist = accumulate(dist_all.begin(), dist_all.end(), 0.0)/dist_all.size();
        histos.fill(HIST("AveDist_withOthers_atDecayPosz"), ave_dist);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ResolutionTask>(cfgc)
  };
}