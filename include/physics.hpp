#pragma once

#include "utils.hpp"

#include <algorithm>
#include <iostream>
#include <math.h>
#include <omp.h>

class PhysObj {
private:
  vec2 pos;
  vec2 prev;
  vec2 targ;
  vec2 acc = vec2(0, 0);
  float friction = 0.98;
  float radius = 10;
  float density = 1;

public:
  PhysObj(vec2 pos, float radius);

  void updatePhysics(float dt);
  void accelerate(vec2 a);

  Color color = BLUE;

  vec2 getPos();
  vec2 getPrev();
  vec2 getTarg();
  float getRadius();
  float getDensity();
  vec2 getVelocity();
  void setPos(vec2 p);
  void setPrev(vec2 p);
  void setTarg(vec2 t);
  void setRadius(float r);
  void setDensity(float d);
};

struct Solver {
  void initializeCache(size_t particleCount);
  std::vector<PhysObj *> objects;
  vec2 g = vec2(0, 1000);
  std::vector<std::vector<float>> interactionCache; // Pairwise distance cache
  void update(float dt, Parameters params);
  void updateColor(int index);
  void applyCollisions(int i, Parameters params);
  void applyForces(int index, vec2 force, Parameters params);
  float calculateDensity(int i, Parameters params);
  vec2 calculatePressureForce(int i, Parameters params);
  void precomputeInteractions(Parameters params);
  ~Solver();
};
