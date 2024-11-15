#include "utils.hpp"
#include <cmath>
#include <raylib.h>

vec2 lerp2D(vec2 a, vec2 b, float t) {
  a.x += (b.x - a.x) * t;
  a.y += (b.y - a.y) * t;
  return a;
}

float lerp1D(float a, float b, float t) {
  a += (b - a) * t;
  return a;
}

float smoothingKernel(float radius, float dist) {
  float volume = PI * pow(radius, 8) / 4;
  float value = std::max(0.0f, radius * radius - dist * dist);
  return value * value * value / volume;
}

float smoothingKernelGradient(float radius, float dist) {
  if (dist >= radius)
    return 0;
  float f = radius * radius - dist * dist;
  float scale = -24 / (PI * pow(radius, 8));
  return scale * dist * f * f;
}

float densityToPressure(float density, Parameters params) {
  float densityError = density - params.targetDensity;
  float pressure = densityError * params.pressureMultiplier;
  return pressure;
}
