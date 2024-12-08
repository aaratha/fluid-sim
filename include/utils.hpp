#pragma once

#include <algorithm>
#include <cmath>

#include "raylib-cpp.hpp"

namespace rl = raylib;
using vec2 = rl::Vector2;

vec2 operator+(const vec2& a, const vec2& b);
vec2 operator-(const vec2& a, const vec2& b);

const int EXT_CIRCLE = 0;
const int INT_CIRCLE = 1;
const int EXT_RECTANGLE = 2;

struct Parameters {
    int screenWidth;
    int screenHeight;
    int particleCount;
    float particleRadius;
    float collisionDamping;
    float friction;
    float gravity;
    float smoothingMultiplier;
    int substeps;
    float targetDensity;
    float pressureMultiplier;
    float maxVelocity;
    float nearPressureMultiplier;
    float viscosity;
    float maxAcceleration;
    float mass;
    float mouseRadius;
    float mouseStrength;
};

// Precompute kernel coefficients during initialization
struct Kernels {
    float poly6;
    float spike;
    float spikePow3;
    float spikePow3Grad;
    float visc;
};

Kernels precomputeKernels(float radius);

vec2 lerp2D(vec2 a, vec2 b, float t);
float lerp1D(float a, float b, float t);

float poly6Kernel(Kernels& kernels, float radius, float dist);
float spikePow3Kernel(Kernels& kernels, float radius, float dist);
float spikePow3GradKernel(Kernels& kernels, float radius, float dist);
float spikeGradKernel(Kernels& kernels, float radius, float dist);
float viscKernel(Kernels& kernels, float radius, float dist);
float densityToPressure(float density, Parameters params);
float nearDensityToNearPressure(float nearDensity, Parameters params);
