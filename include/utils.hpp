#pragma once

#include "raylib-cpp.hpp"

namespace rl = raylib;
using vec2 = rl::Vector2;

vec2 lerp2D(vec2 a, vec2 b, float t);
float lerp1D(float a, float b, float t);
