#pragma once

#include "utils.hpp"

class PhysObj {
private:
  vec2 pos;
  vec2 prev;
  vec2 targ;
  vec2 acc;
  float friction = 0.98;

public:
  PhysObj(vec2 pos);

  void updatePhysics(float dt);
  void accelerate(vec2 a);

  vec2 getPos();
  vec2 getPrev();
  vec2 getTarg();
  void setPos(vec2 p);
  void setTarg(vec2 t);
};

struct Solver {
  std::vector<PhysObj *> objects;
  vec2 g = vec2(0, 1000);
  void update(float dt);
  void updatePositions(float dt);
  void applyCollisions();
  void applyForces();
};
