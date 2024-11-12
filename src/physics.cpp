#include "physics.hpp"

PhysObj::PhysObj(vec2 pos, float radius)
    : pos(pos), prev(pos + vec2(10, 0)), targ(pos), radius(radius) {}

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

void PhysObj::setPrev(vec2 p) { prev = p; }

void PhysObj::setTarg(vec2 p) { targ = p; }

void Solver::update(float dt, Parameters params) {
  for (int i = 0; i < objects.size(); i++) {
    objects[i]->updatePhysics(dt);
    applyForces(i);
    updateColor(i);
    applyCollisions(i, params);
  }
}

void Solver::updateColor(int index) {
  int rg = 50 * int(abs(Vector2Length(objects[index]->getVelocity())));
  objects[index]->color = (Color){static_cast<unsigned char>(rg),
                                  static_cast<unsigned char>(rg), 255, 255};
}

void Solver::applyForces(int index) {
  // Apply gravity when the mouse button is not down
  objects[index]->accelerate(g);
}

void Solver::applyCollisions(int i, Parameters params) {
  PhysObj *obj = objects[i];
  vec2 pos = obj->getPos();
  float radius = obj->getRadius();
  vec2 vel = obj->getVelocity();

  // Wall collisions with bounce effect
  const float bounceFactor =
      params.collisionDamping; // Reduces velocity on bounce
  bool collided = false;

  // Left wall
  if (pos.x - radius < 0) {
    pos.x = radius;
    vel.x = -vel.x * bounceFactor;
    collided = true;
  }
  // Right wall
  else if (pos.x + radius > params.screenWidth) {
    pos.x = params.screenWidth - radius;
    vel.x = -vel.x * bounceFactor;
    collided = true;
  }

  // Top wall
  if (pos.y - radius < 0) {
    pos.y = radius;
    vel.y = -vel.y * bounceFactor;
    collided = true;
  }
  // Bottom wall
  else if (pos.y + radius > params.screenHeight) {
    pos.y = params.screenHeight - radius;
    vel.y = -vel.y * bounceFactor;
    collided = true;
  }

  // Update position and velocity if collision occurred
  if (collided) {
    obj->setPos(pos);
    obj->setPrev(pos - vel); // Update previous position to reflect new velocity
  }

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
}
