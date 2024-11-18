#pragma once

#include "utils.hpp"

#include <algorithm>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <unordered_map>

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

struct GridCell {
  int x, y;
  bool operator==(const GridCell &other) const {
    return x == other.x && y == other.y;
  }
};

// Hash function for GridCell
namespace std {
template <> struct hash<GridCell> {
  size_t operator()(const GridCell &cell) const {
    return hash<int>()(cell.x) ^ (hash<int>()(cell.y) << 1);
  }
};
} // namespace std

class Solver {
private:
  std::vector<std::vector<float>> interactionCache; // Pairwise distance cache
  std::unordered_map<GridCell, std::vector<int>> spatialGrid;

public:
  float cellSize;
  std::vector<PhysObj *> objects;
  void buildSpatialGrid();
  std::vector<int> getNeighbors(int index);
  void initializeCache(size_t particleCount);
  vec2 g = vec2(0, 1000);
  void update(float dt, Parameters params);
  void updateColor(int index);
  void applyCollisions(int i, Parameters params);
  void applyForces(int index, vec2 force, Parameters params);
  float calculateDensity(int i, Parameters params,
                         const std::vector<int> &neighbors);
  vec2 calculatePressureForce(int i, Parameters params,
                              const std::vector<int> &neighbors);
  float calculateSharedPressure(float densityA, float densityB,
                                Parameters params);
  void precomputeInteractions(Parameters params);
  ~Solver();
};
