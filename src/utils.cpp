#include "utils.hpp"
#include <raylib.h>
#include <cmath>

vec2 operator+(const vec2& a, const vec2& b) {
    return Vector2Add(a, b);
}
vec2 operator-(const vec2& a, const vec2& b) {
    return Vector2Subtract(a, b);
}

float lerp1D(float a, float b, float t) {
    return a + (b - a) * t;
}

Kernels precomputeKernels(float radius) {
    Kernels kernels;
    kernels.poly6 = 315 / (64 * PI * pow(radius, 9));
    kernels.spike = -45 / (PI * pow(radius, 6));
    kernels.spikePow3 = 45 / (PI * pow(radius, 6));
    kernels.spikePow3Grad = 135 / (PI * pow(radius, 6));
    kernels.visc = 45 / (PI * pow(radius, 5));
    return kernels;
}

float poly6Kernel(Kernels& kernels, float radius, float dist) {
    if (dist <= 0 || dist >= radius)
        return 0;
    float diff = radius * radius - dist * dist;
    return kernels.poly6 * diff * diff * diff;
}

float spikePow3Kernel(Kernels& kernels, float radius, float dist) {
    if (dist <= 0 || dist >= radius)
        return 0;
    float diff = radius - dist;
    return kernels.spikePow3 * diff * diff * diff;
}

float spikePow3GradKernel(Kernels& kernels, float radius, float dist) {
    if (dist <= 0 || dist >= radius)
        return 0;
    float diff = radius - dist;
    return -kernels.spikePow3Grad * diff * diff;
}

float spikeGradKernel(Kernels& kernels, float radius, float dist) {
    if (dist <= 0 || dist >= radius)
        return 0;
    float diff = radius - dist;
    return kernels.spike * diff * diff;
}

float viscKernel(Kernels& kernels, float radius, float dist) {
    if (dist <= 0 || dist >= radius)
        return 0;
    return kernels.visc * (radius - dist);
}

float densityToPressure(float density, Parameters params) {
    float densityError = density - params.targetDensity;
    return densityError * params.pressureMultiplier;
    // return std::max(0.0f, densityError * params.pressureMultiplier);
}

float nearDensityToNearPressure(float nearDensity, Parameters params) {
    // Linear response for short-range repulsion
    return nearDensity * params.nearPressureMultiplier;
}
