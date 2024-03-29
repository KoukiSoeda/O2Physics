// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief Femtodream Tutorial 1
/// \author Luca Barioglio, Anton Riedel

// O2 includes
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 1
// Example task illustrating how to acess information from different tables

namespace o2::aod
{

// TODO:
// join event selcection and multiplicity table to collision table

// TODO:
// join TPC information to fulltracks table

// TODO:
// create iterators of the joined tables

} // namespace o2::aod

struct CFTutorialTask1 {
  HistogramRegistry histos{
    "Histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    // Define your axes
    // Constant bin width axis
    AxisSpec vtxZAxis = {100, -20, 20};
    // Variable bin width axis
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                     1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0,
                                     2.2, 2.4, 2.8, 3.2, 3.6, 4.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    // Add histograms to histogram manager (as in the output object of in AliPhysics)
    histos.add("hZvtx_before_sel", ";Z (cm)", kTH1F, {vtxZAxis});
    histos.add("hZvtx_after_sel", ";Z (cm)", kTH1F, {vtxZAxis});
    histos.add("hP", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.5, 4.}});
    histos.add("hEta", ";#it{p} (GeV/#it{c})", kTH1F, {{100, -1.5, 1.5}});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hNsigmaTPC", ";#it{p} (GeV/#it{c}); n#sigma_{TPC}^{proton}", kTH2F,
               {{35, 0.5, 4.}, {100, -5., 5.}});
  }

  // Equivalent of the AliRoot task UserExec
  // TODO: use joint tables
  void process(aod::Collision const& coll, aod::Tracks const& inputTracks)
  {
    // Performing the event selection
    histos.fill(HIST("hZvtx_before_sel"), coll.posZ());
    if (fabs(coll.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());

    // Loop over tracks
    for (auto track : inputTracks) {
      if (fabs(track.eta()) > 0.8) {
        continue;
      }
      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
      // TODO
      // fill TPC information into histogram
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<CFTutorialTask1>(cfgc)};
  return workflow;
}
