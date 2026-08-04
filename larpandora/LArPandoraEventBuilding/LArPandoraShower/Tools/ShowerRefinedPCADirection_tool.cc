//############################################################################
//### Name:        ShowerRefinedPCADirection                               ###
//### Author:      ALex Wilkinson                                          ###
//### Date:        23.06.26                                                ###
//### Description: Tool for finding the shower direction using an          ###
//###              iterative weighted PCA on shower space points.          ###
//###              Aim is to find an axis weighted towards the early       ###
//###              and core parts of the shower.                           ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

//C++ Includes
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace ShowerRecoTools {

  class ShowerRefinedPCADirection : public IShowerTool {

  public:
    ShowerRefinedPCADirection(const fhicl::ParameterSet& pset);

    //Calculate the direction of the shower.
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    void InitialiseProducers() override;

    //Function to add the associations
    int AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr,
                        art::Event& Event,
                        reco::shower::ShowerElementHolder& ShowerEleHolder) override;

    struct PCAResult {
      recob::PCAxis pca;
      geo::Point_t centre;
      geo::Vector_t direction;
      int nPoints;
    };

    PCAResult CalculateWeightedPCA(const std::vector<art::Ptr<recob::SpacePoint>>& spacePoints,
                                   const std::vector<double>& weights) const;

    std::vector<double> CalculateWeights(
      const std::vector<art::Ptr<recob::SpacePoint>>& spacePoints,
      geo::Point_t const& startPosition,
      geo::Vector_t const& direction,
      detinfo::DetectorClocksData const* clockData,
      detinfo::DetectorPropertiesData const* detProp,
      art::FindManyP<recob::Hit> const* fmh) const;

    void OrientDirection(geo::Point_t const& startPosition,
                         geo::Point_t const& centre,
                         geo::Vector_t& direction) const;

    //fcl
    art::InputTag fPFParticleLabel;
    int fVerbose;
    bool fWritePCAxis;
    bool fChargeWeighted;
    bool fQuadraticLongitudinalWeight;
    bool fQuadraticTransverseWeight;
    unsigned int fMinPCAPoints;
    unsigned int fMaxIterations;
    double fConvergenceAngle;
    double fTransverseScale;
    double fLongitudinalScale;
    double fDistanceFloor;

    std::string fShowerStartPositionInputLabel;
    std::string fShowerDirectionInputLabel; //Requires an initial seed direction
    std::string fShowerDirectionOutputLabel;
    std::string fShowerCentreOutputLabel;
    std::string fShowerPCAOutputLabel;
  };

  ShowerRefinedPCADirection::ShowerRefinedPCADirection(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fWritePCAxis(pset.get<bool>("WritePCAxis"))
    , fChargeWeighted(pset.get<bool>("ChargeWeighted"))
    , fQuadraticLongitudinalWeight(pset.get<bool>("QuadraticLongitudinalWeight"))
    , fQuadraticTransverseWeight(pset.get<bool>("QuadraticTransverseWeight"))
    , fMinPCAPoints(pset.get<unsigned int>("MinPCAPoints"))
    , fMaxIterations(pset.get<unsigned int>("MaxIterations"))
    , fConvergenceAngle(pset.get<double>("ConvergenceAngle"))
    , fTransverseScale(pset.get<double>("TransverseScale"))
    , fLongitudinalScale(pset.get<double>("LongitudinalScale"))
    , fDistanceFloor(pset.get<double>("DistanceFloor"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
    , fShowerDirectionOutputLabel(pset.get<std::string>("ShowerDirectionOutputLabel"))
    , fShowerCentreOutputLabel(pset.get<std::string>("ShowerCentreOutputLabel"))
    , fShowerPCAOutputLabel(pset.get<std::string>("ShowerPCAOutputLabel"))
  {}

  void ShowerRefinedPCADirection::InitialiseProducers()
  {
    if (!fWritePCAxis) return;
    InitialiseProduct<std::vector<recob::PCAxis>>(fShowerPCAOutputLabel);
    InitialiseProduct<art::Assns<recob::Shower, recob::PCAxis>>("ShowerPCAxisAssn");
    InitialiseProduct<art::Assns<recob::PFParticle, recob::PCAxis>>("PFParticlePCAxisAssn");
  }

  int ShowerRefinedPCADirection::CalculateElement(
    const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerEleHolder)
  {
    if (std::abs(fLongitudinalScale) <= std::numeric_limits<double>::epsilon() ||
        std::abs(fTransverseScale) <= std::numeric_limits<double>::epsilon()) {
      if (fVerbose)
        mf::LogError("ShowerReginedPCADirection")
          << "Longitudinal or transverse scale parameter is too small (" << fLongitudinalScale
          << ", " << fTransverseScale << ")" << std::endl;
      return 1;
    }

    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerReginedPCADirection")
          << "Start position not set. Stopping." << std::endl;
      return 1;
    }

    if (!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerRefinedPCADirection")
          << "Seed direction not set. Stopping." << std::endl;
      return 1;
    }

    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
    const art::FindManyP<recob::SpacePoint>& fmspp =
      ShowerEleHolder.GetFindManyP<recob::SpacePoint>(pfpHandle, Event, fPFParticleLabel);

    std::vector<art::Ptr<recob::SpacePoint>> spacePoints = fmspp.at(pfparticle.key());
    if (spacePoints.size() < fMinPCAPoints) {
      if (fVerbose)
        mf::LogError("ShowerRefinedPCADirection")
          << "Not enough spacepoints for PCA (" << spacePoints.size() << ")" << std::endl;
      return 0;
    }

    geo::Point_t startPosition = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, startPosition);

    geo::Vector_t currentDirection = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel, currentDirection);

    if (currentDirection.R() == 0.) {
      if (fVerbose)
        mf::LogError("ShowerRefinedPCADirection")
          << "Seed direction has zero magnitude" << std::endl;
      return 1;
    }

    auto const spHandle = Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel);
    art::FindManyP<recob::Hit> const* fmh = nullptr;
    detinfo::DetectorClocksData const* clockData = nullptr;
    detinfo::DetectorPropertiesData const* detProp = nullptr;
    if (fChargeWeighted) {
      auto const clockDataForEvent =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
      auto const detPropForEvent =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event,
                                                                                clockDataForEvent);
      clockData = &clockDataForEvent;
      detProp = &detPropForEvent;
      fmh = &(ShowerEleHolder.GetFindManyP<recob::Hit>(spHandle, Event, fPFParticleLabel));
    }

    PCAResult currentResult;
    bool foundWeightedSolution = false;

    for (unsigned int iteration = 0; iteration < fMaxIterations; ++iteration) {
      std::vector<double> weights =
        CalculateWeights(spacePoints, startPosition, currentDirection, clockData, detProp, fmh);
      PCAResult nextResult = CalculateWeightedPCA(spacePoints, weights);

      if (nextResult.nPoints < static_cast<int>(fMinPCAPoints) || nextResult.direction.R() == 0.) {
        if (!foundWeightedSolution) {
          if (fVerbose)
            mf::LogWarning("ShowerRefinedPCADirection")
              << "Unable to calculate a weighted PCA from the seed direction, keeping the seed PCA"
              << std::endl;
          return nextResult.direction.R() == 0. ? 1 : 0;
        }

        if (fVerbose)
          mf::LogWarning("ShowerRefinedPCADirection")
            << "Weighted PCA dropped below the minimum number of usable spacepoints ("
            << fMinPCAPoints << "), "
            << "using the previous successful iteration" << std::endl;

        break;
      }

      OrientDirection(startPosition, nextResult.centre, nextResult.direction);

      const double dot =
        std::clamp(currentDirection.Unit().Dot(nextResult.direction.Unit()), -1., 1.);
      currentResult = nextResult;
      currentDirection = nextResult.direction;
      foundWeightedSolution = true;

      if (std::acos(dot) < fConvergenceAngle || iteration + 1 == fMaxIterations) {
        if (fVerbose)
          mf::LogInfo("ShowerRefinedPCADirection")
            << "Finished PCA refinement at iteration " << iteration + 1 << " / " << fMaxIterations
            << std::endl;
        break;
      }
    }

    //Not implementing errors
    geo::Point_t showerCentreErr = {-999, -999, -999};
    geo::Vector_t showerDirectionErr = {-999, -999, -999};

    ShowerEleHolder.SetElement(currentResult.centre, showerCentreErr, fShowerCentreOutputLabel);
    ShowerEleHolder.SetElement(currentResult.pca, fShowerPCAOutputLabel);
    ShowerEleHolder.SetElement(currentDirection, showerDirectionErr, fShowerDirectionOutputLabel);

    return 0;
  }

  ShowerRefinedPCADirection::PCAResult ShowerRefinedPCADirection::CalculateWeightedPCA(
    const std::vector<art::Ptr<recob::SpacePoint>>& spacePoints,
    const std::vector<double>& weights) const
  {
    if (spacePoints.size() != weights.size())
      throw cet::exception("ShowerRefinedPCADirection")
        << "Spacepoint and weight vectors are different sizes" << std::endl;

    double sumWeights = 0.;
    double x = 0.;
    double y = 0.;
    double z = 0.;
    int nPoints = 0;

    for (std::size_t i = 0; i < spacePoints.size(); ++i) {
      const double weight = weights[i];
      if (weight <= 0.) continue;

      auto const pos = spacePoints[i]->position();
      x += weight * pos.X();
      y += weight * pos.Y();
      z += weight * pos.Z();
      sumWeights += weight;
      ++nPoints;
    }

    if (sumWeights <= 0. || nPoints < 2) return {recob::PCAxis(), {}, {}, nPoints};

    geo::Point_t centre{x / sumWeights, y / sumWeights, z / sumWeights};

    double xx = 0.;
    double yy = 0.;
    double zz = 0.;
    double xy = 0.;
    double xz = 0.;
    double yz = 0.;

    for (std::size_t i = 0; i < spacePoints.size(); ++i) {
      const double weight = weights[i];
      if (weight <= 0.) continue;

      auto const delta = spacePoints[i]->position() - centre;
      xx += weight * delta.X() * delta.X();
      yy += weight * delta.Y() * delta.Y();
      zz += weight * delta.Z() * delta.Z();
      xy += weight * delta.X() * delta.Y();
      xz += weight * delta.X() * delta.Z();
      yz += weight * delta.Y() * delta.Z();
    }

    Eigen::Matrix3d covariance;
    covariance << xx, xy, xz, xy, yy, yz, xz, yz, zz;
    covariance /= sumWeights;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenMatrix(covariance);
    Eigen::Vector3d eigenValuesVector = eigenMatrix.eigenvalues();
    Eigen::Matrix3d eigenVectorsMatrix = eigenMatrix.eigenvectors();

    const double eigenValues[3] = {
      eigenValuesVector(2), eigenValuesVector(1), eigenValuesVector(0)};
    std::vector<std::vector<double>> eigenVectors = {
      {eigenVectorsMatrix(0, 2), eigenVectorsMatrix(1, 2), eigenVectorsMatrix(2, 2)},
      {eigenVectorsMatrix(0, 1), eigenVectorsMatrix(1, 1), eigenVectorsMatrix(2, 1)},
      {eigenVectorsMatrix(0, 0), eigenVectorsMatrix(1, 0), eigenVectorsMatrix(2, 0)}};
    const double avePos[3] = {centre.X(), centre.Y(), centre.Z()};

    recob::PCAxis pca(true, nPoints, eigenValues, eigenVectors, avePos);
    geo::Vector_t direction{
      eigenVectorsMatrix(0, 2), eigenVectorsMatrix(1, 2), eigenVectorsMatrix(2, 2)};

    return {pca, centre, direction, nPoints};
  }

  std::vector<double> ShowerRefinedPCADirection::CalculateWeights(
    const std::vector<art::Ptr<recob::SpacePoint>>& spacePoints,
    geo::Point_t const& startPosition,
    geo::Vector_t const& direction,
    detinfo::DetectorClocksData const* clockData,
    detinfo::DetectorPropertiesData const* detProp,
    art::FindManyP<recob::Hit> const* fmh) const
  {
    std::vector<double> weights(spacePoints.size(), 0.);
    const geo::Vector_t unitDirection = direction.Unit();
    double totalCharge = 0.;

    if (fmh != nullptr) //Means charge weighting is desired
    {
      for (std::size_t i = 0; i < spacePoints.size(); ++i) {
        const auto delta = spacePoints[i]->position() - startPosition;
        const double longitudinalDistance = delta.Dot(unitDirection);
        if (longitudinalDistance < 0.) // Ignore hits "behind" the start position
          continue;

        double charge =
          IShowerTool::GetLArPandoraShowerAlg().SpacePointCharge(spacePoints[i], *fmh);
        double time = IShowerTool::GetLArPandoraShowerAlg().SpacePointTime(spacePoints[i], *fmh);
        charge *=
          std::exp((sampling_rate(*clockData) * time) / (detProp->ElectronLifetime() * 1e3));
        totalCharge += std::max(charge, 0.);
      }
    }

    for (std::size_t i = 0; i < spacePoints.size(); ++i) {
      const auto delta = spacePoints[i]->position() - startPosition;
      double longDist = delta.Dot(unitDirection);
      if (longDist < 0.) // Ignore hits "behind" the start position
        continue;

      double transDist = (delta - longDist * unitDirection).R();

      longDist = std::max(longDist, fDistanceFloor);
      transDist = std::max(transDist, fDistanceFloor);
      const double longWeight = fQuadraticLongitudinalWeight ?
                                  1. / (1. + std::pow(longDist / fLongitudinalScale, 2)) :
                                  1. / (1. + longDist / fLongitudinalScale);
      const double transWeight = fQuadraticTransverseWeight ?
                                   1. / (1. + std::pow(transDist / fTransverseScale, 2)) :
                                   1. / (1. + transDist / fTransverseScale);

      double chargeWeight = 1.;
      if (fmh != nullptr && totalCharge > 0.) {
        double charge =
          IShowerTool::GetLArPandoraShowerAlg().SpacePointCharge(spacePoints[i], *fmh);
        double time = IShowerTool::GetLArPandoraShowerAlg().SpacePointTime(spacePoints[i], *fmh);
        charge *=
          std::exp((sampling_rate(*clockData) * time) / (detProp->ElectronLifetime() * 1e3));
        chargeWeight = std::sqrt(std::max(charge, 0.) / totalCharge);
      }

      weights[i] = transWeight * longWeight * chargeWeight;
    }

    return weights;
  }

  void ShowerRefinedPCADirection::OrientDirection(geo::Point_t const& startPosition,
                                                  geo::Point_t const& centre,
                                                  geo::Vector_t& direction) const
  {
    if ((centre - startPosition).R() == 0. || direction.R() == 0.) return;

    const geo::Vector_t generalDirection = (centre - startPosition).Unit();
    if (direction.Dot(generalDirection) < 0.) direction *= -1.;
  }

  int ShowerRefinedPCADirection::AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr,
                                                 art::Event& Event,
                                                 reco::shower::ShowerElementHolder& ShowerEleHolder)
  {
    if (!fWritePCAxis) return 0;

    if (!ShowerEleHolder.CheckElement(fShowerPCAOutputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerRefinedPCADirection: Add Assns") << "PCA not set." << std::endl;
      return 1;
    }

    int ptrSize = GetVectorPtrSize(fShowerPCAOutputLabel);
    const art::Ptr<recob::PCAxis> pcaPtr =
      GetProducedElementPtr<recob::PCAxis>(fShowerPCAOutputLabel, ShowerEleHolder, ptrSize - 1);
    const art::Ptr<recob::Shower> showerPtr =
      GetProducedElementPtr<recob::Shower>("shower", ShowerEleHolder);

    AddSingle<art::Assns<recob::Shower, recob::PCAxis>>(showerPtr, pcaPtr, "ShowerPCAxisAssn");
    AddSingle<art::Assns<recob::PFParticle, recob::PCAxis>>(pfpPtr, pcaPtr, "PFParticlePCAxisAssn");

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerRefinedPCADirection)
