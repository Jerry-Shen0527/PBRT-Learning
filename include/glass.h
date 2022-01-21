
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_GLASS_H
#define PBRT_MATERIALS_GLASS_H

// materials/glass.h*
//#include "pbrt.h"
#include "rtweekend.h"
#include "material.h"
#include "spectrum.h"

namespace pbrt {

// GlassMaterial Declarations
class GlassMaterial : public Material {
  public:
    // GlassMaterial Public Methods
    GlassMaterial(const std::shared_ptr<Texture<Spectrum>> &Kr,
                  const std::shared_ptr<Texture<Spectrum>> &Kt,
                  const std::shared_ptr<Texture<Float>> &uRoughness,
                  const std::shared_ptr<Texture<Float>> &vRoughness,
                  const std::shared_ptr<Texture<Float>> &index,
                  const std::shared_ptr<Texture<Float>> &bumpMap,
                  bool remapRoughness)
        : Kr(Kr),
          Kt(Kt),
          uRoughness(uRoughness),
          vRoughness(vRoughness),
          index(index),
          bumpMap(bumpMap),
          remapRoughness(remapRoughness) {}
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // GlassMaterial Private Data
    std::shared_ptr<Texture<Spectrum>> Kr, Kt;
    std::shared_ptr<Texture<Float>> uRoughness, vRoughness;
    std::shared_ptr<Texture<Float>> index;
    std::shared_ptr<Texture<Float>> bumpMap;
    bool remapRoughness;
};

//GlassMaterial *CreateGlassMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_GLASS_H
