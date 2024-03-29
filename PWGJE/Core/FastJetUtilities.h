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

/// \file FastJetUtilities.h
/// \brief Jet related utilities that require fastjet
///
/// \author Nima Zardoshti

#ifndef PWGJE_CORE_FASTJETUTILITIES_H_
#define PWGJE_CORE_FASTJETUTILITIES_H_

#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>

#include "PWGJE/Core/JetFinder.h"

enum class JetConstituentStatus {
  track = 0,
  cluster = 1,
  candidateHF = 2,
};

namespace FastJetUtilities
{

// Class defined to store additional info which is passed to the FastJet object
class fastjet_user_info : public fastjet::PseudoJet::UserInfoBase
{
  int status; // the status of each particle (Options are: TrueParticle (final state particles in generator event which arent special), HFParticle (heavy-flavour particle of interest in generator event), ThermalParticle (particles belonging to the thermal backgound), DecaySisterParticle (other particles poduced in the decay resulting in a non-prompt heavy-flavour particle of interest))
  int index;  // a number unique to each particle in the event

 public:
  fastjet_user_info()
  {
    status = -9;
    index = -9;
  }
  fastjet_user_info(int _status, int _index)
  {
    status = _status;
    index = _index;
  }
  ~fastjet_user_info() = default;
  void setStatus(int set) { status = set; }
  void setIndex(int set) { index = set; }
  int getStatus() const { return status; }
  int getIndex() const { return index; }
};

/**
 * Set the fastjet_user_info object when filling the jet constituents.
 *
 * @param constituents vector of constituents to be clustered.
 * @param index global index of constituent
 * @param status status of constituent type
 */

void setFastJetUserInfo(std::vector<fastjet::PseudoJet>& constituents, int index = -99999999, int status = static_cast<int>(JetConstituentStatus::track))
{
  fastjet_user_info* user_info = new fastjet_user_info(status, index); // FIXME: can setting this as a pointer be avoided?
  constituents.back().set_user_info(user_info);
  if (index != -99999999) { // FIXME: needed for constituent subtraction as user_info is not propagated, but need to be quite careful to make sure indices dont overlap between tracks, clusters and HF candidates. Current solution might not be optimal
    int i = index;
    if (status == static_cast<int>(JetConstituentStatus::track)) {
      i = i + 1;
    }
    if (status == static_cast<int>(JetConstituentStatus::cluster)) {
      i = -1 * (i + 1);
    }
    if (status == static_cast<int>(JetConstituentStatus::candidateHF)) {
      i = 0;
    }
    constituents.back().set_user_index(i); // FIXME: needed for constituent subtraction, but need to be quite careful to make sure indices dont overlap between tracks, clusters and HF candidates. Current solution might not be optimal
  }
}

/**
 * Add track as a pseudojet object to the fastjet vector
 *
 * @param constituent constituent to be added
 * @param constituents vector of constituents
 * @param index global index of constituent
 * @param status status of constituent type
 * @param status mass hypothesis for constituent
 */

template <typename T>
void fillTracks(const T& constituent, std::vector<fastjet::PseudoJet>& constituents, int index = -99999999, int status = static_cast<int>(JetConstituentStatus::track), double mass = JetFinder::mPion)
{
  if (status == static_cast<int>(JetConstituentStatus::track) || status == static_cast<int>(JetConstituentStatus::candidateHF)) {
    auto p = std::sqrt((constituent.px() * constituent.px()) + (constituent.py() * constituent.py()) + (constituent.pz() * constituent.pz()));
    auto energy = std::sqrt((p * p) + (mass * mass));
    constituents.emplace_back(constituent.px(), constituent.py(), constituent.pz(), energy);
  }
  setFastJetUserInfo(constituents, index, status);
}

/**
 * Add cluster as a pseudojet object to the fastjet vector
 *
 * @param constituent constituent to be added
 * @param constituents vector of constituents
 * @param index global index of constituent
 * @param status status of constituent type
 */

template <typename T>
void fillClusters(const T& constituent, std::vector<fastjet::PseudoJet>& constituents, int index = -99999999, int status = static_cast<int>(JetConstituentStatus::cluster))
{
  if (status == static_cast<int>(JetConstituentStatus::cluster)) {
    double clusterpt = constituent.energy() / std::cosh(constituent.eta());
    constituents.emplace_back(clusterpt * std::cos(constituent.phi()), clusterpt * std::sin(constituent.phi()), clusterpt * std::sinh(constituent.eta()), constituent.energy());
  }
  setFastJetUserInfo(constituents, index, status);
}

}; // namespace FastJetUtilities

#endif // PWGJE_CORE_FASTJETUTILITIES_H_
