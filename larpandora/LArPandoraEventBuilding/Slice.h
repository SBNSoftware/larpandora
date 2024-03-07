/**
 *  @file   larpandora/LArPandoraEventBuilding/Slice.h
 *
 *  @brief  header for the lar pandora slice class
 */

#ifndef LAR_PANDORA_SLICE_H
#define LAR_PANDORA_SLICE_H 1

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace lar_pandora {

  /**
 *  @brief Slice class
 */
  class Slice {
  public:
    /**
     *  @brief  Default constructor
     *
     *  @param  topologicalScore the topological score from Pandora
     *  @param  targetHypothesis the slice as reconstructed under the target hypothesis
     *  @param  crHypothesis the slice as reconstructed under the cosmic-ray hypothesis
     *  @param  isTarget if the slice has been identified as a target
     */
    Slice(const float topologicalScore,
          const PFParticleVector& targetHypothesis,
          const PFParticleVector& crHypothesis,
          const bool isTarget = false);

    /**
     *  @brief Get the topological score for the slice - closer to 1 means more likely to be the target slice
     */
    float GetTopologicalScore() const;

    /**
     *  @brief Get the slice as reconstructed under the target hypothesis
     */
    const PFParticleVector& GetTargetHypothesis() const;

    /**
     *  @brief Get the slice as reconstructed under the cosmic-ray hypothesis
     */
    const PFParticleVector& GetCosmicRayHypothesis() const;

    /**
     *  @brief Check if the slice has been identified as a target
     */
    bool IsTaggedAsTarget() const;

    /**
     *  @brief Tag the slice as a neutrino / test beam particle
     */
    void TagAsTarget();

    /**
     *  @brief Tag the slice as a cosmic
     */
    void TagAsCosmic();

  private:
    float m_topologicalScore; ///< The topological neutrino / beam particle score from Pandora
    PFParticleVector
      m_targetHypothesis; ///< The slice as reconstructed under the neutrino / beam particle hypothesis
    PFParticleVector m_crHypothesis; ///< The slice as reconstructed under the cosmic-ray hypothesis
    bool m_isTarget; ///< If the slice has been identified as a neutrino / beam particle
  };

  typedef std::vector<Slice> SliceVector;

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline Slice::Slice(const float topologicalScore,
                      const PFParticleVector& targetHypothesis,
                      const PFParticleVector& crHypothesis,
                      const bool isTarget)
    : m_topologicalScore(topologicalScore)
    , m_targetHypothesis(targetHypothesis)
    , m_crHypothesis(crHypothesis)
    , m_isTarget(isTarget)
  {}

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float Slice::GetTopologicalScore() const
  {
    return m_topologicalScore;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline const PFParticleVector& Slice::GetTargetHypothesis() const
  {
    return m_targetHypothesis;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline const PFParticleVector& Slice::GetCosmicRayHypothesis() const
  {
    return m_crHypothesis;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline bool Slice::IsTaggedAsTarget() const
  {
    return m_isTarget;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline void Slice::TagAsTarget()
  {
    m_isTarget = true;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline void Slice::TagAsCosmic()
  {
    m_isTarget = false;
  }

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_SLICE_H
