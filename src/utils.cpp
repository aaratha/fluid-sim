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
  if (dist >= radius)
    return 0;

  float volume = (PI * pow(radius, 4)) / 6;
  return (radius - dist) * (radius - dist) / volume;
}

float smoothingKernelGradient(float radius, float dist) {
  if (dist >= radius)
    return 0;
  float scale = 12 / (pow(radius, 4) * PI);
  return (dist - radius) * scale;
}

float densityToPressure(float density, Parameters params) {
  float densityError = density - params.targetDensity;
  float pressure = densityError * params.pressureMultiplier;
  return pressure;
}
