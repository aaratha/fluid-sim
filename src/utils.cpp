#include "utils.hpp"
#include <cmath>
#include <raylib.h>

Kernels precomputeKernels(float radius) {
  Kernels kernels;
  kernels.poly6 = 315 / (64 * PI * pow(radius, 9));
  kernels.spike = -45 / (PI * pow(radius, 6));
  kernels.visc = 15 / (2 * PI * radius * radius * radius);
  return kernels;
}

float poly6Kernel(Kernels &kernels, float radius, float dist) {
  if (dist <= 0 || dist >= radius)
    return 0;
  float diff = radius * radius - dist * dist;
  return kernels.poly6 * diff * diff * diff;
}

float spikeGradKernel(Kernels &kernels, float radius, float dist) {
  if (dist <= 0 || dist >= radius)
    return 0;
  return kernels.spike * (radius - dist);
}

float viscKernel(Kernels &kernels, float radius, float dist) {
  if (dist <= 0 || dist >= radius)
    return 0;
  float r2 = radius * radius;
  float r3 = r2 * radius;
  return kernels.visc * (-(dist * dist * dist) / (2 * r3) + (dist * dist) / r2 +
                         radius / (2 * dist) - 1);
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
}

float nearDensityToNearPressure(float nearDensity, Parameters params) {
  // Linear response for short-range repulsion
  return nearDensity * params.nearPressureMultiplier;
}
