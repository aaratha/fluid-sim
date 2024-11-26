#pragma once

#include "utils.hpp"

#include <algorithm>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <unordered_map>
#include <vector>

// Structure for GridCell remains the same
struct GridCell {
  int x, y;
  bool operator==(const GridCell &other) const {
    return x == other.x && y == other.y;
  }
};

namespace std {
template <> struct hash<GridCell> {
  size_t operator()(const GridCell &cell) const {
    return hash<int>()(cell.x) ^ (hash<int>()(cell.y) << 1);
  }
};
} // namespace std

struct Forces {
  vec2 tension;
  vec2 pressure;
  vec2 viscosity;
};

class Solver {
private:
  std::vector<std::vector<float>> interactionCache; // Pairwise distance cache
  std::unordered_map<GridCell, std::vector<int>> spatialGrid;

public:
  float cellSize;
  float smoothingRadius;
  float radius = 10;
  void buildSpatialGrid();

  std::vector<vec2> positions;
  std::vector<vec2> predictedPositions;
  std::vector<vec2> velocities;
  std::vector<float> densities;
  std::vector<float> nearDensities;
  std::vector<Color> colors;

  Kernels kernels;

  std::vector<int> getNeighbors(size_t index);
  void initializeCache(size_t particleCount);
  vec2 g = vec2(0, 1000);
  void update(float dt, Parameters params);
  void updateColor(size_t index, Parameters params);
  void applyCollisions(size_t i, Parameters params);
  void applyForces(size_t index, vec2 force, Parameters params);
  float calculateDensity(size_t i, Parameters params,
                         const std::vector<int> &neighbors);
  float calculateNearDensity(size_t i, Parameters params,
                             const std::vector<int> &neighbors);
  Forces calculateForces(size_t i, Parameters params,
                         const std::vector<int> &neighbors);

  float calculateSharedPressure(float densityA, float densityB,
                                Parameters params);
  void precomputeInteractions(Parameters params);
  ~Solver();
};
