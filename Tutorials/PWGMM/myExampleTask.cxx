/// \author Koki Soeda
/// \since 10/05/2023

#include <iostream>
#include <cmath>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "TDatabasePDG.h"
#include "MathUtils/Utils.h"

using namespace std;
using namespace o2;
using namespace o2::framework;

struct MyAnalysisTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsVer{"nBinsVer", 300, "N bins in Collisino Vertex histo"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisVertex{nBinsVer, -30, +30, "cm"};
    const AxisSpec axistrackType{6, -0.5, 5.5, "TrackType"};
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisChi2(505, -4.5, +500.5, "");
    const AxisSpec axisPDG(2001, -1000.5, +1000.5, "");
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}"};

    // create histograms
    histos.add("CollisionPos", "CollisionPos", kTH1F, {axisVertex});
    histos.add("TrackType", "TrackType", kTH1F, {axistrackType});
    histos.add("FwdParticleCounter", "FwdParticleCounter", kTH1F, {axisCounter});
    histos.add("FwdMuonCounter", "FwdMuonCounter", kTH1F, {axisCounter});
    histos.add("GlobalMuonTrackPDG", "GlobalMuonTrackPDG", kTH1F, {axisPDG});
    histos.add("Chi2GMT", "Chi2GMT", kTH1F, {axisChi2});
    histos.add("muonMomPDG", "muonMomPDG", kTH1F, {axisPDG});
    histos.add("MFTParticleCounter", "MFTParticleCounter", kTH1F, {axisCounter});
    histos.add("MFTMuonCounter", "MFTMuonCounter", kTH1F, {axisCounter});
    
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
      histos.fill(HIST("CollisionPos"), mcCol_z);

      //Information of FwdTracks
      for(auto& fwdtrack : fwdtracks){
        histos.fill(HIST("TrackType"), fwdtrack.trackType());

        if(!fwdtrack.has_mcParticle()) continue;
        auto mcParticle_fwd = fwdtrack.mcParticle();
        histos.fill(HIST("FwdParticleCounter"), 0.5);
        auto fwdTrackPDG = mcParticle_fwd.pdgCode();

        if(fabs(fwdTrackPDG)==13) histos.fill(HIST("FwdMuonCounter"), 0.5);

        if(fwdtrack.trackType()==0){//Required MFT-MCH-MID
          auto chi2GMT = fwdtrack.chi2MatchMCHMFT();
          histos.fill(HIST("GlobalMuonTrackPDG"), fwdTrackPDG);
          histos.fill(HIST("Chi2GMT"), chi2GMT);
          if(chi2GMT==-1) continue;
          if(fabs(fwdTrackPDG)!=13) continue;
          if(mcParticle_fwd.has_mothers()){
            auto mcMom = MCparticle.rawIteratorAt(mcParticle_fwd.mothersIds()[0]);
            auto mcMomPDG = mcMom.pdgCode();
            histos.fill(HIST("muonMomPDG"), mcMomPDG);
          }
        }
      }
      //Information of MFTTracks
      for(auto& mfttrack : mfttracks){
        if(!mfttrack.has_mcParticle()) continue;
        auto mcParticle_mft = mfttrack.mcParticle();
        histos.fill(HIST("MFTParticleCounter"), 0.5);
        if(fabs(mcParticle_mft.pdgCode())==13) histos.fill(HIST("MFTMuonCounter"), 0.5);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyAnalysisTask>(cfgc)};
}
