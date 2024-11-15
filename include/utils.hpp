#pragma once

#include <algorithm>
#include <cmath>

#include "raylib-cpp.hpp"

namespace rl = raylib;
using vec2 = rl::Vector2;

struct Parameters {
  int screenWidth;
  int screenHeight;
  int particleCount;
  float particleRadius;
  float collisionDamping;
  float friction;
  vec2 gravity;
  float smoothingRadius;
  int substeps;
  float targetDensity;
  float pressureMultiplier;
};

vec2 lerp2D(vec2 a, vec2 b, float t);
float lerp1D(float a, float b, float t);

float smoothingKernel(float radius, float dist);
float smoothingKernelGradient(float radius, float dist);
float densityToPressure(float density, Parameters params);
