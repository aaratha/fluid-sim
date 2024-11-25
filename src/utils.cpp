#include "utils.hpp"
#include <cmath>
#include <raylib.h>

float smoothingKernel(float radius, float dist) {
  if (0 >= dist || dist >= radius)
    return 0;

  float consts = 315 / (64 * PI * pow(radius, 9));
  return consts * pow((radius * radius - dist * dist), 3);
}

float nearSmoothingKernel(float radius, float dist) {
  if (0 >= dist || dist >= radius)
    return 0;

  // Spiky kernel with more balanced short-range behavior
  return -45 / (PI * pow(radius, 6));
}

float smoothingKernelGradient(float radius, float dist) {
  if (0 >= dist || dist >= radius)
    return 0;

  return -(45 / (PI * pow(radius, 6))) * (radius - dist);
}

float densityToPressure(float density, Parameters params) {
  float densityError = density - params.targetDensity;
  return densityError * params.pressureMultiplier;
}

float nearDensityToNearPressure(float nearDensity, Parameters params) {
  // Linear response for short-range repulsion
  return nearDensity * params.nearPressureMultiplier;
}
