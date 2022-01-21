
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_MATTE_H
#define PBRT_MATERIALS_MATTE_H

// materials/matte.h*
#include "rtweekend.h"
#include "material.h"
#include "spectrum.h"

namespace pbrt {
    //class Spectrum;
// MatteMaterial Declarations
class MatteMaterial : public Material {
  public:
    // MatteMaterial Public Methods
    MatteMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,
                  const std::shared_ptr<Texture<Float>> &sigma,
                  const std::shared_ptr<Texture<Float>> &bumpMap)
        : Kd(Kd), sigma(sigma), bumpMap(bumpMap) {}
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // MatteMaterial Private Data
    std::shared_ptr<Texture<Spectrum>> Kd;
    std::shared_ptr<Texture<Float>> sigma, bumpMap;
};

//MatteMaterial *CreateMatteMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_MATTE_H
