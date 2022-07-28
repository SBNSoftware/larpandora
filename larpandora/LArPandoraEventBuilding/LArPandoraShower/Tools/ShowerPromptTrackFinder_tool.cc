//############################################################################
//### Name:        ShowerPromptTrackFinder                                 ###
//### Author:      Bruno Zamorano                                          ###
//### Date:        19.07.22                                                ###
//### Description: Tool to find the prompt, 'track-like' bit of a shower.  ###
//###              Based on work made by Marina Bravo on nu-e scattering.  ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"

//Root Includes
#include "TPrincipal.h"

namespace ShowerRecoTools {

  class ShowerPromptTrackFinder : public IShowerTool {

  public:
    ShowerPromptTrackFinder(const fhicl::ParameterSet& pset);

    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    std::vector<art::Ptr<recob::SpacePoint>> RunIncrementalSpacePointFinder(
      const art::Event& Event,
      std::vector<art::Ptr<recob::SpacePoint>> const& sps,
      const art::FindManyP<recob::Hit>& fmh);

    void AddSpacePointsToSegment(std::vector<art::Ptr<recob::SpacePoint>>& segment,
                                 std::vector<art::Ptr<recob::SpacePoint>>& sps_pool,
                                 size_t num_sps_to_take);


    //Function to calculate the track direction using a charge-weighted 3D PCA calculation.
    void TrackPCA(std::vector<art::Ptr<recob::SpacePoint>>& sps);

    void TrackPCA(const detinfo::DetectorClocksData& clockData,
                             const detinfo::DetectorPropertiesData& detProp,
                             const std::vector<art::Ptr<recob::SpacePoint>>& sps,
                             const art::FindManyP<recob::Hit>& fmh);

    std::vector<TVector3> fPCAdirs;
    std::vector<float>   fPCAeigenValues;
    int fBestIndex;
    float fBestEigenValue;

    //Services
    art::InputTag fPFParticleLabel;
    int fVerbose;
    bool fUseShowerDirection;
    long unsigned int fStartFitSize;
    bool fChargeWeighted;
    bool fForwardHitsOnly;

    std::string fShowerStartPositionInputLabel;
    std::string fShowerDirectionInputLabel;
    std::string fInitialTrackHitsOutputLabel;
    std::string fInitialTrackSpacePointsOutputLabel;
  };

  ShowerPromptTrackFinder::ShowerPromptTrackFinder(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fUseShowerDirection(pset.get<bool>("UseShowerDirection"))
    , fStartFitSize(pset.get<int>("StartFitSize"))
    , fChargeWeighted(pset.get<bool>("ChargeWeighted"))
    , fForwardHitsOnly(pset.get<bool>("ForwardHitsOnly"))
    , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
    , fInitialTrackHitsOutputLabel(pset.get<std::string>("InitialTrackHitsOutputLabel"))
    , fInitialTrackSpacePointsOutputLabel(
        pset.get<std::string>("InitialTrackSpacePointsOutputLabel"))
  {
    if (fStartFitSize == 0) {
      throw cet::exception("ShowerPromptTrackFinder")
        << "We cannot make a track if you don't give us at leats one hit. Change fStartFitSize "
           "please to something sensible";
    }
  }

  int
  ShowerPromptTrackFinder::CalculateElement(
    const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerPromptTrackFinder")
          << "Start position not set, returning " << std::endl;
      return 1;
    }

    // Get the assocated pfParicle Handle
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    // Get the spacepoint - PFParticle assn
    const art::FindManyP<recob::SpacePoint>& fmspp =
      ShowerEleHolder.GetFindManyP<recob::SpacePoint>(pfpHandle, Event, fPFParticleLabel);

    // Get the spacepoints
    auto const spHandle = Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel);

    // Get the hits associated with the space points
    const art::FindManyP<recob::Hit>& fmh =
      ShowerEleHolder.GetFindManyP<recob::Hit>(spHandle, Event, fPFParticleLabel);

    // Get the SpacePoints
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints = fmspp.at(pfparticle.key());

    //We cannot progress with no spacepoints.
    if (spacePoints.empty()) {
      if (fVerbose)
        mf::LogError("ShowerPromptTrackFinder")
          << "No space points, returning " << std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, ShowerStartPosition);

    //Decide if the you want to use the direction of the shower or make one.
    if (fUseShowerDirection) {

      if (!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)) {
        if (fVerbose)
          mf::LogError("ShowerPromptTrackFinder")
            << "Direction not set, returning " << std::endl;
        return 1;
      }

      TVector3 ShowerDirection = {-999, -999, -999};
      ShowerEleHolder.GetElement(fShowerDirectionInputLabel, ShowerDirection);

      //Order the spacepoints
      IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePoints(
        spacePoints, ShowerStartPosition, ShowerDirection);
      //Remove the back hits if requird.
      if (fForwardHitsOnly) {
        int back_sps = 0;
        for (auto spacePoint : spacePoints) {
          double proj = IShowerTool::GetLArPandoraShowerAlg().SpacePointProjection(
            spacePoint, ShowerStartPosition, ShowerDirection);
          if (proj < 0) { ++back_sps; }
          if (proj > 0) { break; }
        }
        spacePoints.erase(spacePoints.begin(), spacePoints.begin() + back_sps);
      }
    }
    else {
      //Order the spacepoint using the magnitude away from the vertex
      IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePoints(spacePoints,
                                                                   ShowerStartPosition);
    }

    if (spacePoints.size() < 3) {
      if (fVerbose)
        mf::LogError("ShowerPromptTrackFinder")
          << "Not enough spacepoints. Bailing" << std::endl;
      return 1;
    }

    //Actually run the algorithm.
    std::vector<art::Ptr<recob::SpacePoint>> track_sps =
      RunIncrementalSpacePointFinder(Event, spacePoints, fmh);

    // Get the hits associated to the space points and seperate them by planes
    std::vector<art::Ptr<recob::Hit>> trackHits;
    for (auto const& spacePoint : track_sps) {
      std::vector<art::Ptr<recob::Hit>> hits = fmh.at(spacePoint.key());
      for (auto const& hit : hits) {
        trackHits.push_back(hit);
      }
    }

    //Add to the holder
    ShowerEleHolder.SetElement(trackHits, fInitialTrackHitsOutputLabel);
    ShowerEleHolder.SetElement(track_sps, fInitialTrackSpacePointsOutputLabel);

    return 0;
  }

  void ShowerPromptTrackFinder::TrackPCA(std::vector<art::Ptr<recob::SpacePoint>>& sps)
  {

    //Initialise the the PCA.
    TPrincipal* pca = new TPrincipal(3, "");

    //Normalise the spacepoints, charge-weight and add to the PCA.
    for (auto& sp : sps) {

      TVector3 sp_position = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp);

      double sp_coord[3];
      sp_coord[0] = sp_position.X();
      sp_coord[1] = sp_position.Y();
      sp_coord[2] = sp_position.Z();

      //Add to the PCA
      pca->AddRow(sp_coord);
    }

    //Evaluate the PCA
    pca->MakePrincipals();

    //Get the Eigenvectors and eigenvalues.
    const TMatrixD* Eigenvectors = pca->GetEigenVectors();
    const TVectorD* Eigenvalues  = pca->GetEigenValues();

    TVector3 Eigenvector = {(*Eigenvectors)[0][0], (*Eigenvectors)[1][0], (*Eigenvectors)[2][0]};
    float Eigenvalue = (*Eigenvalues)[0];

    delete pca;

    fPCAdirs.push_back(Eigenvector);
    fPCAeigenValues.push_back(Eigenvalue);
  }

  //Function to calculate the shower direction using a charge-weighted 3D PCA calculation.
  void ShowerPromptTrackFinder::TrackPCA(
    const detinfo::DetectorClocksData& clockData,
    const detinfo::DetectorPropertiesData& detProp,
    const std::vector<art::Ptr<recob::SpacePoint>>& sps,
    const art::FindManyP<recob::Hit>& fmh)
  {

    //Initialise the the PCA.
    TPrincipal* pca = new TPrincipal(3, "");

    float TotalCharge = 0;

    //Normalise the spacepoints, charge-weight and add to the PCA.
    for (auto& sp : sps) {

      TVector3 sp_position = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp);

      float wht = 1;

      if (fChargeWeighted) {

        //Get the charge.
        float Charge = IShowerTool::GetLArPandoraShowerAlg().SpacePointCharge(sp, fmh);
        //        std::cout << "Charge: " << Charge << std::endl;

        //Get the time of the spacepoint
        float Time = IShowerTool::GetLArPandoraShowerAlg().SpacePointTime(sp, fmh);

        //Correct for the lifetime at the moment.
        Charge *= std::exp((sampling_rate(clockData) * Time) / (detProp.ElectronLifetime() * 1e3));
        //        std::cout << "Charge: "<< Charge << std::endl;

        //Charge Weight
        wht *= std::sqrt(Charge / TotalCharge);
      }

      double sp_coord[3];
      sp_coord[0] = sp_position.X() * wht;
      sp_coord[1] = sp_position.Y() * wht;
      sp_coord[2] = sp_position.Z() * wht;

      //Add to the PCA
      pca->AddRow(sp_coord);
    }

    //Evaluate the PCA
    pca->MakePrincipals();

    //Get the Eigenvectors.
    const TMatrixD* Eigenvectors = pca->GetEigenVectors();
    const TVectorD* Eigenvalues  = pca->GetEigenValues();

    TVector3 Eigenvector = {(*Eigenvectors)[0][0], (*Eigenvectors)[1][0], (*Eigenvectors)[2][0]};
    float Eigenvalue = (*Eigenvalues)[0];

    delete pca;

    fPCAdirs.push_back(Eigenvector);
    fPCAeigenValues.push_back(Eigenvalue);
  }

  std::vector<art::Ptr<recob::SpacePoint>>
  ShowerPromptTrackFinder::RunIncrementalSpacePointFinder(
    const art::Event& Event,
    std::vector<art::Ptr<recob::SpacePoint>> const& sps,
    const art::FindManyP<recob::Hit>& fmh)
  {

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    //Create space point pool (yes we are copying the input vector because we're going to twiddle with it
    std::vector<art::Ptr<recob::SpacePoint>> sps_pool = sps;

    while (sps_pool.size() >= fStartFitSize) {

      std::vector<art::Ptr<recob::SpacePoint>> track_segment;
      AddSpacePointsToSegment(track_segment, sps_pool, (size_t)(fStartFitSize));

      // Now do the fit
      if (fChargeWeighted) TrackPCA(clockData, detProp, track_segment, fmh);
      else TrackPCA(track_segment);

      if(fPCAeigenValues.back() > fBestEigenValue){
        fBestEigenValue = fPCAeigenValues.back();
        fBestIndex = fPCAeigenValues.size()-1;
      }
    }

    std::vector<art::Ptr<recob::SpacePoint>> sps_pool_copy = sps;
    std::vector<art::Ptr<recob::SpacePoint>> best_track;

    for(int i = 0; i <= fBestIndex; i++){
      AddSpacePointsToSegment(best_track, sps_pool_copy, (size_t)(fStartFitSize));
    }

    return best_track;
  }

  void
  ShowerPromptTrackFinder::AddSpacePointsToSegment(
    std::vector<art::Ptr<recob::SpacePoint>>& segment,
    std::vector<art::Ptr<recob::SpacePoint>>& sps_pool,
    size_t num_sps_to_take)
  {
    size_t new_segment_size = segment.size() + num_sps_to_take;
    while (segment.size() < new_segment_size && sps_pool.size() > 0) {
      segment.push_back(sps_pool[0]);
      sps_pool.erase(sps_pool.begin());
    }
    return;
  }

}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPromptTrackFinder)
