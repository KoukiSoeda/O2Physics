/// \author Koki Soeda
/// \since 30.01.2023

#include <iostream>
#include <cmath>
#include <vector>
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "Common/Core/RecoDecay.h"
#include "TDatabasePDG.h"
#include "MathUtils/Utils.h"


using namespace std;
using namespace o2; 
using namespace o2::framework;
using namespace o2::soa;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

struct DCAandPCA {
   Service<TDatabasePDG> pdg;
   Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
   Configurable<float> zMax{"zMax", 50., "value for Zvtx cut"};
   
   HistogramRegistry registry{
      "registry",
      {
         {"mcCollision_Vertex", " ; z(cm); ", {HistType::kTH1F, {{10000, -25, 25}}}},
         {"Distance_D0", " ; z(cm); ", {HistType::kTH1F, {{1000, -0.5, 5}}}},         
         {"muTrack_dcax", " ; x(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"muTrack_dcay", " ; y(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"muTrack_mftx", " ; x(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"muTrack_mfty", " ; y(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"D0decayx", " ; x(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"D0decayy", " ; y(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"D0decayz", " ; z(cm); ", {HistType::kTH1F, {{5000, -50, 50}}}},
         {"pcax", " ; x(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"pcay", " ; y(cm); ", {HistType::kTH1F, {{10000, -100, 100}}}},
         {"pcaz", " ; z(cm); ", {HistType::kTH1F, {{5000, -50, 50}}}},
         {"PIDpurity", "", {HistType::kTH1F, {{102, -0.5, 100.5}}}},
         {"EstimatePDGcode", "", {HistType::kTH1F, {{2002, -1000.5, 1000.5}}}},
      }
   };

   void processMCDCA(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision,
                     MFTTracksLabeled const& tracks,
                     aod::McParticles const& particleMC,
                     aod::McCollisions const&)
   {
      if(collision.has_mcCollision()){
         if((collision.mcCollision().posZ() < zMax) && (collision.mcCollision().posZ() > -zMax)){
            if(!useEvSel || (useEvSel && collision.sel8())){
               float dcaXY;
               for(auto& track0 : tracks){
                  if(!track0.has_mcParticle()) continue;
                  auto particle0 = track0.mcParticle();
                  auto t_collision = track0.collision();
                  auto mc_collision = t_collision.globalIndex();
                  auto mc_col_verz = collision.mcCollision().posZ();
                  registry.fill(HIST("mcCollision_Vertex"), mc_col_verz);
                  //cout << mc_collision << "     :   Track0" << endl;

                  if(fabs(particle0.pdgCode())!=13) continue;
                  int estpdg;
                  int64_t truthK_ID,muonID,estID;
                  float vx_mu,vy_mu,vz_mu,px_mu,py_mu,pz_mu,t_mu,mft_mu_x,mft_mu_y,s_mu,dca_mu_x,dca_mu_y,mft_mu_z,dca_mu_z;
                  //float mft_det_x, mft_det_y, mft_det_z;
                  float vx_can,vy_can,vz_can,px_can,py_can,pz_can,t_can,mft_can_x,mft_can_y,mft_can_z,s_can,dca_can_x,dca_can_y,dca_can_z;
                  float s,t,a,b,c,d,e,f,g,h;
                  float mu_x,mu_y,mu_z,can_x,can_y,can_z,pca_x,pca_y,pca_z;
                  if(particle0.mcCollisionId() == collision.mcCollision().globalIndex()){
                     if(particle0.has_mothers()){
                        int signMUON = 0;
                        float closest_track = 10000;
                        auto mcMom = particleMC.rawIteratorAt(particle0.mothersIds()[0]);
                        auto Daughters = mcMom.daughters_as<aod::McParticles>();

                        if(fabs(mcMom.pdgCode())==421){
                           registry.fill(HIST("D0decayx"), particle0.vx());
                           registry.fill(HIST("D0decayy"), particle0.vy());
                           registry.fill(HIST("D0decayz"), particle0.vz());
                           registry.fill(HIST("Distance_D0"), fabs(mc_col_verz - particle0.vz()));
                        }

                        //Actual particles information
                        if(fabs(mcMom.pdgCode())==421){
                           signMUON++;
                           for(auto Daughter : Daughters){
                              if(fabs(Daughter.pdgCode())==321){
                                 truthK_ID = Daughter.globalIndex();
                              }
                              if(fabs(Daughter.pdgCode())==13){
                                 vx_mu = particle0.vx();
                                 vy_mu = particle0.vy();
                                 vz_mu = particle0.vz();
                                 px_mu = particle0.px();
                                 py_mu = particle0.py();
                                 pz_mu = particle0.pz();
                                 auto propaZ = track0.z();

                                 t_mu = (propaZ-vz_mu)/pz_mu;
                                 mft_mu_x = vx_mu + t_mu*px_mu;
                                 mft_mu_y = vy_mu + t_mu*py_mu;
                                 mft_mu_z = vz_mu + t_mu*pz_mu;
                                 s_mu = (mc_col_verz-vz_mu)/pz_mu;
                                 dca_mu_x = vx_mu + s_mu*px_mu;
                                 dca_mu_y = vy_mu + s_mu*py_mu;
                                 dca_mu_z = vz_mu + s_mu*pz_mu;
                                 
                                 //cout << mc_col_verz << ",  " << dca_mu_z << endl;

                                 if(particle0.globalIndex()==Daughter.globalIndex()) muonID = Daughter.globalIndex(); 
                                 //cout << "mftX: " << mft_det_x << "  mftY: " << mft_det_y << " mftZ: " << mft_det_z << endl;
                                 a = dca_mu_x;
                                 b = mft_mu_x;
                                 c = dca_mu_y;
                                 d = mft_mu_y;
                                 cout << b << "," << d << ",-46," << vx_mu << "," << vy_mu << "," << vz_mu << endl;
                                 //cout << "MFT_MC_x: " << b << ", MFT_x: " << track0.x() << ", MFT_MC_y: " << d << ", MFT_y: " << track0.y() << ", MFT_z: " << track0.z() <<endl;
                                 //cout << "mft_MC_x: " << b << "  mft_MC_y: " << d << endl;
                                 registry.fill(HIST("muTrack_dcax"), a);
                                 registry.fill(HIST("muTrack_dcay"), c);
                                 registry.fill(HIST("muTrack_mftx"), b);
                                 registry.fill(HIST("muTrack_mfty"), d);
                                 //cout << dca_mu_x << "," << dca_mu_y << "," << dca_mu_z << "," << mft_mu_x << "," << mft_mu_y << "," << mft_mu_z << endl;
                              }
                           }
                        }
                        if(signMUON==0) continue;
                        for(auto& track1 : tracks){
                           if(!track1.has_mcParticle()) continue;
                           auto t1_collision = track1.collision();
                           auto mc1_collision = t1_collision.globalIndex();
                           //cout << mc1_collision << endl;
                           
                           auto particle1 = track1.mcParticle();
                           if(particle1.globalIndex()==muonID) continue;
                           vx_can = particle1.vx();
                           vy_can = particle1.vy();
                           vz_can = particle1.vz();
                           px_can = particle1.px();
                           py_can = particle1.py();
                           pz_can = particle1.pz();

                           t_can = (-46-vz_can)/pz_can;
                           mft_can_x = vx_can + t_can*px_can;
                           mft_can_y = vy_can + t_can*py_can;
                           mft_can_z = vz_can + t_can*pz_can;
                           s_can = (mc_col_verz-vz_can)/pz_can;
                           dca_can_x = vx_can + s_can*px_can;
                           dca_can_y = vy_can + s_can*py_can;
                           dca_can_z = vz_can + s_can*pz_can;
                           e = dca_can_x;
                           f = mft_can_x;
                           g = dca_can_y;
                           h = mft_can_y;

                           float p11,p12,p13,p21,p22,p23,u11,u12,u13,u21,u22,u23;
                           u11 = b-a;
                           u12 = d-c;
                           u13 = mft_mu_z-dca_mu_z;
                           u21 = f-e;
                           u22 = h-g;
                           u23 = mft_can_z-dca_can_z;
                           p11 = a;
                           p12 = c;
                           p13 = dca_mu_z;
                           p21 = e;
                           p22 = g;
                           p23 = dca_can_z;

                           float norm = sqrt(pow((u12*u23-u13*u22),2)+pow((u13*u21-u11*u23),2)+pow((u11*u22-u12*u21),2));
                           float dot = sqrt(pow((u12*u23-u13*u22)*(p11-p21)+(u13*u21-u11*u23)*(p12-p22)+(u11*u22-u12*u21)*(p13-p23),2));
                           float pro = u11*u21+u12*u22+u13*u23;
                           if(dot!=0){
                              float dist = dot/norm;
                              if(dist<=closest_track){
                                 closest_track = dist;
                                 pca_x = (((p11+pro*u21)*(p21-p11))/(1-pro*pro))*u11;
                                 pca_y = (((p12+pro*u22)*(p22-p12))/(1-pro*pro))*u12;
                                 pca_z = (((p13+pro*u23)*(p23-p13))/(1-pro*pro))*u13;
                              }
                           }
                        /*
                           s = ((-e*f+a*f-g*h+c*h)*(b*f+d*h-pow(46,2)) + (e*b-a*b+d*g-c*d)*(pow(f,2)+pow(h,2)+pow(46,2))) / ((pow(b,2)+pow(d,2)+pow(46,2))*(pow(f,2)+pow(h,2)+pow(46,2)) - pow(b*f+d*h-pow(46,2),2));
                           t = ((b*f+d*h-pow(46,2))*s - e*f + a*f - g*h + c*h) / (pow(f,2)+pow(h,2)+pow(46,2));

                           mu_x = dca_mu_x+s*mft_mu_x;
                           mu_y = dca_mu_y+s*mft_mu_y;
                           mu_z = -46*s;
                           can_x = dca_can_x+t*mft_can_x;
                           can_y = dca_can_y+t*mft_can_y;
                           can_z = -46*t;

                           //float distance_track = pow(mu_x-can_x,2)+pow(mu_y-can_y,2)+pow(mu_z-can_z,2);
                           float distance_track = fabs((46*((a-e)*(h-d)+(c-g)*(b-f))) / sqrt(pow(46*(h-d),2)+pow(46*(b-f),2)+pow(b*h-d*f,2)));
                           
                           if(distance_track<=closest_track){
                              closest_track = distance_track;
                              pca_x = (mu_x+can_x)/2;
                              pca_y = (mu_y+can_y)/2;
                              pca_z = (mu_z+can_z)/2;
                              estpdg = particle1.pdgCode();
                              estID = particle1.globalIndex();
                           }
                        */
                        }
                        //cout << "Closest Distance: " << closest_track << ", x: " << pca_x << ", y: " << pca_y << ", z: " << pca_z << endl;
                        registry.fill(HIST("pcax"), pca_x);
                        registry.fill(HIST("pcay"), pca_y);
                        registry.fill(HIST("pcaz"), pca_z);
                        registry.fill(HIST("PIDpurity"), fabs(estID-truthK_ID));
                        registry.fill(HIST("EstimatePDGcode"), estpdg);
                        //cout << "Distance of each track: " << closest_track << endl;
                     }
                  }
               }
            }
         }
      }      
   }
   PROCESS_SWITCH(DCAandPCA, processMCDCA, "Get the DCAxy of MC information for MFT", true);
};

WorkflowSpec
   defineDataProcessing(ConfigContext const& cfgc)
{
   return WorkflowSpec{adaptAnalysisTask<DCAandPCA>(cfgc)};
}