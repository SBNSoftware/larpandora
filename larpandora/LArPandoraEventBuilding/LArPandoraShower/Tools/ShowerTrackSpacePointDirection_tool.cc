//############################################################################
//### Name:        ShowerTrackSpacePointDirection                          ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using the         ###
//###              the average direction of theinitial track spacepoints.  ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"

namespace ShowerRecoTools {

  class ShowerTrackSpacePointDirection : IShowerTool {

  public:
    ShowerTrackSpacePointDirection(const fhicl::ParameterSet& pset);

    //Calculate the direction using the initial track spacepoints
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    //fcl
    int fVerbose;
    bool fUsePandoraVertex; //Direction from point defined as
    //(Position of SP - Vertex) rather than
    //(Position of SP - Track Start Point).

    std::string fInitialTrackSpacePointsInputLabel;
    std::string fShowerStartPositionInputLabel;
    std::string fInitialTrackInputLabel;
    std::string fShowerDirectionOutputLabel;
  };

  ShowerTrackSpacePointDirection::ShowerTrackSpacePointDirection(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fVerbose(pset.get<int>("Verbose"))
    , fUsePandoraVertex(pset.get<bool>("UsePandoraVertex"))
    , fInitialTrackSpacePointsInputLabel(pset.get<std::string>("InitialTrackSpacePointsInputLabel"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fInitialTrackInputLabel(pset.get<std::string>("InitialTrackInputLabel"))
    , fShowerDirectionOutputLabel(pset.get<std::string>("ShowerDirectionOutputLabel"))
  {}

  int ShowerTrackSpacePointDirection::CalculateElement(
    const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    //Check the Track Hits has been defined
    if (!ShowerEleHolder.CheckElement(fInitialTrackSpacePointsInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerTrackSpacePointDirection")
          << "Initial track spacepoints not set" << std::endl;
      return 0;
    }

    //Check the start position is set.
    if (fUsePandoraVertex && !ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerTrackSpacePointDirection")
          << "Start position not set, returning " << std::endl;
      return 0;
    }

    //Get the start poistion
    geo::Point_t StartPosition = {-999, -999, -999};
    if (fUsePandoraVertex) {
      ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, StartPosition);
    }
    else {
      //Check the Tracks has been defined
      if (!ShowerEleHolder.CheckElement(fInitialTrackInputLabel)) {
        if (fVerbose)
          mf::LogError("ShowerTrackSpacePointDirection") << "Initial track not set" << std::endl;
        return 0;
      }
      recob::Track InitialTrack;
      ShowerEleHolder.GetElement(fInitialTrackInputLabel, InitialTrack);
      geo::Point_t Start_point = InitialTrack.Start();
      StartPosition = {Start_point.X(), Start_point.Y(), Start_point.Z()};
    }

    //Get the initial track hits.
    std::vector<art::Ptr<recob::SpacePoint>> intitaltrack_sp;
    ShowerEleHolder.GetElement(fInitialTrackSpacePointsInputLabel, intitaltrack_sp);

    //Calculate the mean direction
    geo::Vector_t meanDir = {0., 0., 0.};

    //Get the spacepoints associated to the track hit
    int N = 0;
    for (auto const& sp : intitaltrack_sp) {

      //Get the direction relative to the start positon
      auto const dir = (sp->position() - StartPosition).Unit();
      if (dir.R() == 0) { continue; }
      meanDir += dir;
      ++N;
    }

    if (N > 0) {
      meanDir = meanDir.Unit();
      //NOTE need to include a direction otherwise this will not override any previous shower direciton tha did include an error :(
      geo::Vector_t meanDirErr = {-999, -999, -999}; //Not implementing an error
      ShowerEleHolder.SetElement(meanDir, meanDirErr, fShowerDirectionOutputLabel);
    }
    else {
      if (fVerbose)
        mf::LogError("ShowerTrackSpacePointDirection")
          << "No points found to calculate a direction" << std::endl;
      return 1;
    }

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackSpacePointDirection)
