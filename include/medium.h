#pragma once
#ifndef MEDIUM_H
#define MEDIUM_H

#include "spectrum.h"

// Medium Declarations
class Medium {
public:
    // Medium Interface
    virtual ~Medium() {}
    virtual Color Tr(const ray& ray, Sampler& sampler) const = 0;
    virtual Color Sample(const ray& ray, Sampler& sampler,
        MemoryArena& arena,
        MediumInteraction* mi) const = 0;
};

// Medium Declarations
class Medium {
public:
    // Medium Interface
    virtual ~Medium() {}
    virtual Color Tr(const ray& ray, Sampler& sampler) const = 0;
    virtual Color Sample(const ray& ray, Sampler& sampler,
        MemoryArena& arena,
        MediumInteraction* mi) const = 0;
};
#endif // !MEDIUM_H
