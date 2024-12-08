#pragma once

#include "utils.hpp"

#include <math.h>
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>

// Structure for GridCell remains the same
struct GridCell {
    int x, y;
    bool operator==(const GridCell& other) const {
        return x == other.x && y == other.y;
    }
};

namespace std {
template <>
struct hash<GridCell> {
    size_t operator()(const GridCell& cell) const {
        return hash<int>()(cell.x) ^ (hash<int>()(cell.y) << 1);
    }
};
}  // namespace std

struct Forces {
    vec2 tension;
    vec2 pressure;
    vec2 viscosity;
    vec2 mouse;
};

class Obstacle {
   public:
    vec2 position;
    int shape;
    float radius;
    Rectangle rectangle;

    Obstacle(vec2 position, int shape, float radius, Rectangle rectangle);
    void lerpRadius(float target, float transition_time, float dt);
    float getRadius();
    void setRadius(float r);
    vec2 getPos();
    void setPos(vec2 pos);
};

class Solver {
   private:
    std::vector<std::vector<float>>
        interactionCache;  // Pairwise distance cache
    std::unordered_map<GridCell, std::vector<int>> spatialGrid;
    std::vector<float> densities;
    std::vector<vec2> velocities;
    std::vector<float> nearDensities;
    Kernels kernels;

   public:
    float cellSize;
    float smoothingRadius;
    float radius = 10;

    std::vector<Obstacle> obstacles;
    std::vector<vec2> positions;
    std::vector<vec2> predictedPositions;
    std::vector<Color> colors;

    void buildSpatialGrid();

    std::vector<int> getNeighbors(size_t index);
    void initializeCache(size_t particleCount);
    void update(float dt, Parameters params);
    void updateColor(size_t index, Parameters params);
    void applyCollisions(size_t i, Parameters params);
    void applyForces(size_t index, vec2 force, Parameters params);
    std::pair<float, float> calculateDensity(size_t i,
                                             Parameters params,
                                             const std::vector<int>& neighbors);
    Forces calculateForces(size_t i,
                           Parameters params,
                           const std::vector<int>& neighbors);
    void precomputeInteractions(Parameters params);
    ~Solver();
};
