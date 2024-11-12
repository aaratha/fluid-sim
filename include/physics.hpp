#pragma once

#include "utils.hpp"

class PhysObj {
private:
  vec2 pos;
  vec2 prev;
  vec2 targ;
  vec2 acc;
  float friction = 0.98;
  float radius = 10;

public:
  PhysObj(vec2 pos);

  void updatePhysics(float dt);
  void accelerate(vec2 a);

  Color color = BLUE;

  vec2 getPos();
  vec2 getPrev();
  vec2 getTarg();
  float getRadius();
  vec2 getVelocity();
  void setPos(vec2 p);
  void setTarg(vec2 t);
};

struct Solver {
  std::vector<PhysObj *> objects;
  vec2 g = vec2(0, 1000);
  void update(float dt, int screenWidth, int screenHeight);
  void updateColor();
  void updatePositions(float dt);
  void applyCollisions(int screenWidth, int screenHeight);
  void applyForces();
};
