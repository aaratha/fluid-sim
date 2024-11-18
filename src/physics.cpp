#include "physics.hpp"
#include "utils.hpp"
#include <algorithm>

// for init velocity: prev(pos + vec2(10, 0))
PhysObj::PhysObj(vec2 pos, float radius)
    : pos(pos), prev(pos), targ(pos), radius(radius) {}

void PhysObj::updatePhysics(float dt) {
  const vec2 vel = pos - prev;
  prev = pos;
  pos += (vel * friction) + acc * dt * dt;
  acc = vec2(0, 0);
  if (pos.x != pos.x) {
    std::cout << "pos: " << pos.x << ", " << pos.y << std::endl;
    CloseWindow();
  }
}

void PhysObj::accelerate(vec2 a) { acc += a; }

vec2 PhysObj::getPos() { return pos; }

vec2 PhysObj::getPrev() { return prev; }

vec2 PhysObj::getTarg() { return targ; }

float PhysObj::getRadius() { return radius; }

vec2 PhysObj::getVelocity() { return pos - prev; }

float PhysObj::getDensity() { return density; }

void PhysObj::setPos(vec2 p) { pos = p; }

void PhysObj::setPrev(vec2 p) { prev = p; }

void PhysObj::setTarg(vec2 p) { targ = p; }

void PhysObj::setRadius(float r) { radius = r; }

void PhysObj::setDensity(float d) { density = d; }

void Solver::initializeCache(size_t particleCount) {
  interactionCache.resize(particleCount,
                          std::vector<float>(particleCount, 0.0f));
}

void Solver::update(float dt, Parameters params) {
  precomputeInteractions(params); // Fill interactionCache and spatialGrid

  // Precompute neighbors
  std::vector<std::vector<int>> neighborLists(objects.size());
  for (int i = 0; i < objects.size(); ++i) {
    neighborLists[i] = getNeighbors(i);
  }

#pragma omp parallel for
  for (int i = 0; i < objects.size(); ++i) {
    vec2 velocity = objects[i]->getPos() - objects[i]->getPrev();
    objects[i]->setTarg(objects[i]->getPos() + velocity * 1 / 120.0f);

    cellSize = params.smoothingRadius;
    objects[i]->setDensity(calculateDensity(i, params, neighborLists[i]));
    objects[i]->setRadius(params.particleRadius);
    vec2 force = calculatePressureForce(i, params, neighborLists[i]);
    applyForces(i, force, params);
    objects[i]->updatePhysics(dt);
    applyCollisions(i, params);
    updateColor(i);
  }
}

void Solver::updateColor(int index) {
  int rg = 50 * int(abs(Vector2Length(objects[index]->getVelocity())));
  objects[index]->color = (Color){static_cast<unsigned char>(rg),
                                  static_cast<unsigned char>(rg), 255, 255};
}

void Solver::applyForces(int index, vec2 force, Parameters params) {
  // Apply gravity when the mouse button is not down
  objects[index]->accelerate(force / objects[index]->getDensity());
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
    vel.y *= -bounceFactor;
    collided = true;
  }
  // Bottom wall
  else if (pos.y + radius > params.screenHeight) {
    pos.y = params.screenHeight - radius;
    vel.y *= -bounceFactor;
    collided = true;
  }

  // Update position and velocity if collision occurred
  if (collided) {
    obj->setPos(pos);
    obj->setPrev(pos - vel); // Update previous position to reflect new velocity
  }
}

float Solver::calculateDensity(int i, Parameters params,
                               const std::vector<int> &neighbors) {
  float density = 0;

#pragma omp parallel for reduction(+ : density)
  for (int j : neighbors) {
    if (i == j)
      continue; // Can be removed if neighbors are pre-filtered

    float dist = Vector2Length(objects[j]->getPos() - objects[i]->getTarg());

    if (dist > params.smoothingRadius)
      continue;

    float influence = smoothingKernel(params.smoothingRadius, dist);
    density += influence;
  }

  return std::clamp(density, 1.0e-6f, 100.0f);
}

vec2 Solver::calculatePressureForce(int i, Parameters params,
                                    const std::vector<int> &neighbors) {
  vec2 pressureForce = vec2(0, 0);

#pragma omp parallel for reduction(- : pressureForce)
  for (int j : neighbors) {
    if (i == j)
      continue;

    vec2 offset = (objects[j]->getPos() - objects[i]->getTarg());
    float dist = Vector2Length(offset);

    if (dist > params.smoothingRadius)
      continue;

    vec2 dir = dist == 0 ? vec2(GetRandomValue(-2, 2), GetRandomValue(-2, 2))
                         : offset / dist;
    float slope = smoothingKernelGradient(params.smoothingRadius, dist);
    float density = objects[j]->getDensity();
    float sharedPressure =
        calculateSharedPressure(density, objects[i]->getDensity(), params);

    pressureForce.x -= dir.x * sharedPressure * slope / density;
    pressureForce.y -= dir.y * sharedPressure * slope / density;
  }
  return pressureForce;
}

float Solver::calculateSharedPressure(float densityA, float densityB,
                                      Parameters params) {
  float pressureA = densityToPressure(densityA, params);
  float pressureB = densityToPressure(densityB, params);
  return (pressureA + pressureB) / 2;
}

void Solver::precomputeInteractions(Parameters params) { buildSpatialGrid(); }

void Solver::buildSpatialGrid() {
  spatialGrid.clear();
  for (int i = 0; i < objects.size(); ++i) {
    vec2 pos = objects[i]->getPos();
    int cellX = static_cast<int>(pos.x / cellSize);
    int cellY = static_cast<int>(pos.y / cellSize);
    GridCell cell{cellX, cellY};
    spatialGrid[cell].push_back(i);
  }
}

// Get neighboring particles
std::vector<int> Solver::getNeighbors(int index) {
  std::vector<int> neighbors;
  vec2 pos = objects[index]->getPos();
  int cellX = static_cast<int>(pos.x / cellSize);
  int cellY = static_cast<int>(pos.y / cellSize);

  // Check the cell and surrounding cells
  const int cellDist = 2;
  for (int dx = -cellDist; dx <= cellDist; ++dx) {
    for (int dy = -cellDist; dy <= cellDist; ++dy) {
      GridCell neighborCell{cellX + dx, cellY + dy};
      if (spatialGrid.find(neighborCell) != spatialGrid.end()) {
        const auto &cellParticles = spatialGrid[neighborCell];
        neighbors.insert(neighbors.end(), cellParticles.begin(),
                         cellParticles.end());
      }
    }
  }
  return neighbors;
}

Solver::~Solver() {
  for (auto *obj : objects) {
    delete obj;
  }
}
