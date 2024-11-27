#include "utils.hpp"
#include <cmath>
#include <raylib.h>

vec2 operator+(const vec2 &a, const vec2 &b) { return Vector2Add(a, b); }
vec2 operator-(const vec2 &a, const vec2 &b) { return Vector2Subtract(a, b); }

Kernels precomputeKernels(float radius) {
  Kernels kernels;
  kernels.poly6 = 315 / (64 * PI * pow(radius, 9));
  kernels.spike = -45 / (PI * pow(radius, 6));
  kernels.visc = 45 / (PI * pow(radius, 5));
  return kernels;
}

float poly6Kernel(Kernels &kernels, float radius, float dist) {
  if (dist <= 0 || dist >= radius)
    return 0;
  float diff = radius * radius - dist * dist;
  return kernels.poly6 * diff * diff * diff * 100000;
}

float spikeGradKernel(Kernels &kernels, float radius, float dist) {
  if (dist <= 0 || dist >= radius)
    return 0;
  float diff = radius - dist;
  return kernels.spike * diff * diff * 100000;
}

float viscKernel(Kernels &kernels, float radius, float dist) {
  if (dist <= 0 || dist >= radius)
    return 0;
  return kernels.visc * (radius - dist) * 100000;
}

float nearSmoothingKernel(float radius, float dist) {
  if (0 >= dist || dist >= radius)
    return 0;

  // Spiky kernel with more balanced short-range behavior
  return -45 / (PI * pow(radius, 6));
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
