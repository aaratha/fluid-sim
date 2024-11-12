#include "physics.hpp"

PhysObj::PhysObj(vec2 pos) : pos(pos), prev(pos), targ(pos) {}

void PhysObj::updatePhysics(float dt) {
  const vec2 vel = pos - prev;
  prev = pos;
  pos += (vel * friction) + acc * dt * dt;
  acc = vec2(0, 0);
}

void PhysObj::accelerate(vec2 a) { acc += a; }

vec2 PhysObj::getPos() { return pos; }

vec2 PhysObj::getPrev() { return prev; }

vec2 PhysObj::getTarg() { return targ; }

void PhysObj::setPos(vec2 p) { pos = p; }

void PhysObj::setTarg(vec2 p) { targ = p; }

void Solver::update(float dt) {
  applyForces();
  updatePositions(dt);
  applyCollisions();
}

void Solver::updatePositions(float dt) {
  for (auto obj : objects) {
    obj->updatePhysics(dt);
  }
}

void Solver::applyCollisions() {}

void Solver::applyForces() {
  for (auto obj : objects) {
    // Apply gravity when the mouse button is not down
    obj->accelerate(g);
  }
}
