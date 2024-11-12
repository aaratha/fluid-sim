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

float PhysObj::getRadius() { return radius; }

vec2 PhysObj::getVelocity() { return pos - prev; }

void PhysObj::setPos(vec2 p) { pos = p; }

void PhysObj::setTarg(vec2 p) { targ = p; }

void Solver::update(float dt, int screenWidth, int screenHeight) {
  applyForces();
  updatePositions(dt);
  applyCollisions(screenWidth, screenHeight);
}

void Solver::updatePositions(float dt) {
  for (auto obj : objects) {
    obj->updatePhysics(dt);
  }
}

void Solver::applyCollisions(int screenWidth, int screenHeight) {
  for (int i = 0; i < objects.size(); i++) {
    // Object-object collisions
    for (int j = i + 1; j < objects.size(); j++) {
      vec2 diff = objects[i]->getPos() - objects[j]->getPos();
      float dist = Vector2Length(diff);
      float overlap = objects[i]->getRadius() + objects[j]->getRadius() - dist;
      if (overlap > 0) {
        vec2 normal = Vector2Normalize(diff);
        objects[i]->setPos(objects[i]->getPos() + normal * overlap / 2.0f);
        objects[j]->setPos(objects[j]->getPos() - normal * overlap / 2.0f);
      }
    }

    // Object-wall collisions
    float radius = objects[i]->getRadius();
    vec2 pos = objects[i]->getPos();
    if (pos.x - radius < 0) {
      objects[i]->setPos(vec2(radius, pos.y));
      objects[i]->accelerate(vec2(-objects[i]->getVelocity().x, 0));
    } else if (pos.x + radius > screenWidth) {
      objects[i]->setPos(vec2(screenWidth - radius, pos.y));
      objects[i]->accelerate(vec2(-objects[i]->getVelocity().x, 0));
    }
    if (pos.y - radius < 0) {
      objects[i]->setPos(vec2(pos.x, radius));
      objects[i]->accelerate(vec2(0, -objects[i]->getVelocity().y));
    } else if (pos.y + radius > screenHeight) {
      objects[i]->setPos(vec2(pos.x, screenHeight - radius));
      objects[i]->accelerate(vec2(0, -objects[i]->getVelocity().y));
    }
  }
}

void Solver::applyForces() {
  for (auto obj : objects) {
    // Apply gravity when the mouse button is not down
    obj->accelerate(g);
  }
}
