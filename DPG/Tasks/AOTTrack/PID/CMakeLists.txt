# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

# PID QA

## TOF
o2physics_add_dpl_workflow(pid-tof-qa
                           SOURCES qaPIDTOF.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-tof-qa-mc
                           SOURCES qaPIDTOFMC.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

## TPC
o2physics_add_dpl_workflow(pid-tpc-qa
                           SOURCES qaPIDTPC.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-tpc-qa-mc
                           SOURCES qaPIDTPCMC.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)