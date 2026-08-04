//############################################################################
//### Name:        ShowerTrackHitDirection                                 ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using the         ###
//###              initial track the average direction of the initial hits ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"

namespace ShowerRecoTools {

  class ShowerTrackHitDirection : IShowerTool {

  public:
    ShowerTrackHitDirection(const fhicl::ParameterSet& pset);

    //Calculate the shower direction from the initial track hits.
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    //fcl
    int fVerbose;
    bool fUsePandoraVertex; //Direction from point defined as (Position of Hit - Vertex)
    //rather than (Position of Hit - Track Start Point)
    art::InputTag fHitModuleLabel;
    art::InputTag fPFParticleLabel;

    std::string fInitialTrackHitsInputLabel;
    std::string fShowerStartPositionInputLabel;
    std::string fInitialTrackInputLabel;
    std::string fShowerDirectionOutputLabel;
  };

  ShowerTrackHitDirection::ShowerTrackHitDirection(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fVerbose(pset.get<int>("Verbose"))
    , fUsePandoraVertex(pset.get<bool>("UsePandoraVertex"))
    , fHitModuleLabel(pset.get<art::InputTag>("HitModuleLabel"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fInitialTrackHitsInputLabel(pset.get<std::string>("InitialTrackHitsInputLabel"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fInitialTrackInputLabel(pset.get<std::string>("InitialTrackInputLabel"))
    , fShowerDirectionOutputLabel(pset.get<std::string>("ShowerDirectionOutputLabel"))
  {}

  int ShowerTrackHitDirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                                art::Event& Event,
                                                reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    //Check the Track Hits has been defined
    if (!ShowerEleHolder.CheckElement(fInitialTrackHitsInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerTrackHitDirection") << "Initial track hits not set" << std::endl;
      return 0;
    }

    //Check the start position is set.
    if (fUsePandoraVertex && !ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerTrackHitDirection")
          << "Start position not set, returning " << std::endl;
      return 0;
    }

    //Get the start poistion
    geo::Point_t StartPosition = {-999, -999, -99};
    if (fUsePandoraVertex) {
      ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, StartPosition);
    }
    else {
      //Check the Tracks has been defined
      if (!ShowerEleHolder.CheckElement(fInitialTrackInputLabel)) {
        if (fVerbose)
          mf::LogError("ShowerTrackHitDirection") << "Initial track not set" << std::endl;
        return 0;
      }
      recob::Track InitialTrack;
      ShowerEleHolder.GetElement(fInitialTrackInputLabel, InitialTrack);
      geo::Point_t Start_point = InitialTrack.Start();
      StartPosition = {Start_point.X(), Start_point.Y(), Start_point.Z()};
    }

    //Get the spacepoints associated to hits
    auto const hitHandle = Event.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel);

    //Get the spacepoint handle. We need to do this in 3D.
    const art::FindManyP<recob::SpacePoint>& fmsp =
      ShowerEleHolder.GetFindManyP<recob::SpacePoint>(hitHandle, Event, fPFParticleLabel);

    //Get the initial track hits.
    std::vector<art::Ptr<recob::Hit>> InitialTrackHits;
    ShowerEleHolder.GetElement(fInitialTrackHitsInputLabel, InitialTrackHits);

    //Calculate the mean direction
    geo::Vector_t meanDir = {0., 0., 0.};

    //Get the spacepoints associated to the track hit
    std::vector<art::Ptr<recob::SpacePoint>> intitaltrack_sp;
    int N = 0;
    for (auto const hit : InitialTrackHits) {
      std::vector<art::Ptr<recob::SpacePoint>> sps = fmsp.at(hit.key());
      for (auto const sp : sps) {
        intitaltrack_sp.push_back(sp);

        //Get the direction relative to the start positon
        auto const dir = (sp->position() - StartPosition).Unit();
        if (dir.R() == 0) { continue; }
        meanDir += dir;
        ++N;
      }
    }

    if (N > 0) {
      meanDir = meanDir.Unit();
      //NOTE need to include a direction otherwise this will not override any previous shower direciton tha did include an error :(
      geo::Vector_t meanDirErr = {-999, -999, -999}; //Not implementing an error
      ShowerEleHolder.SetElement(meanDir, meanDirErr, fShowerDirectionOutputLabel);
    }
    else {
      if (fVerbose)
        mf::LogError("ShowerTrackHitDirection")
          << "No points found to calculate a direction" << std::endl;
      return 1;
    }

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackHitDirection)
