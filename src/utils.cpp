#include "utils.hpp"

vec2 lerp2D(vec2 a, vec2 b, float t) {
  a.x += (b.x - a.x) * t;
  a.y += (b.y - a.y) * t;
  return a;
}

float lerp1D(float a, float b, float t) {
  a += (b - a) * t;
  return a;
}
