#include "physics.hpp"
#include "utils.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <raylib.h>

const float EPS = 1e-6f;

void Solver::initializeCache(size_t particleCount) {
  interactionCache.resize(particleCount,
                          std::vector<float>(particleCount, 0.0f));
}

void Solver::update(float dt, Parameters params) {
  // Cache frequently accessed values
  smoothingRadius = params.smoothingMultiplier * params.particleRadius;
  kernels = precomputeKernels(smoothingRadius);
  cellSize = smoothingRadius;

  // Update spatial partitioning
  buildSpatialGrid();

  // Pre-allocate neighbor lists
  std::vector<std::vector<int>> neighborLists(positions.size());

// Step 1: Update positions and find neighbors
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < positions.size(); ++i) {
    // Add gravity
    velocities[i].y += params.gravity * dt;

    // Get neighbors first as they're needed for multiple calculations
    neighborLists[i] = getNeighbors(i);

    // Predict next position for density calculations
    predictedPositions[i] = positions[i] + velocities[i] * dt;
  }

// Step 2: Calculate densities
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < positions.size(); ++i) {
    std::pair<float, float> density_pair =
        calculateDensity(i, params, neighborLists[i]);
    densities[i] = density_pair.first;
    nearDensities[i] = density_pair.second;
  }

// Step 3: Calculate and apply forces
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < positions.size(); ++i) {
    Forces forces = calculateForces(i, params, neighborLists[i]);

    // Apply pressure acceleration
    vec2 acc =
        (forces.pressure + forces.tension + forces.viscosity + forces.mouse) /
        densities[i]; // / std::max(densities[i], 1e-6f);
    float acc_limit = params.maxAcceleration;
    float acc_mag = Vector2Length(acc);
    if (acc_mag > acc_limit * acc_limit) {
      acc = acc / sqrt(acc_mag) * acc_limit;
    }
    velocities[i] += acc * dt;

    // Apply velocity constraints
    float currentSpeed = Vector2Length(velocities[i]);
    if (currentSpeed > params.maxVelocity) {
      velocities[i] =
          Vector2Scale(Vector2Normalize(velocities[i]), params.maxVelocity);
    }
  }

// Step 4: Final position update and boundary handling
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < positions.size(); ++i) {
    // Apply friction
    velocities[i] = Vector2Scale(velocities[i], params.friction);

    // Update color based on velocity
    updateColor(i, params);

    // Update position
    positions[i] += velocities[i] * dt;

    // Handle collisions with boundaries
    applyCollisions(i, params);
  }
}

void Solver::updateColor(size_t index, Parameters params) {
  float velocityMagnitude = Vector2Length(velocities[index]);
  float normalizedVelocity =
      std::min(velocityMagnitude / params.maxVelocity, 1.0f);
  int colorValue = static_cast<int>(55 + normalizedVelocity * 200);
  colorValue = std::clamp(colorValue, 0, 255);

  colors[index] = Color{static_cast<unsigned char>(colorValue),
                        static_cast<unsigned char>(colorValue), 255, 255};
}

void Solver::applyCollisions(size_t i, Parameters params) {
  vec2 &pos = positions[i];
  vec2 &vel = velocities[i];
  float radius = params.particleRadius;
  float bounceFactor = params.collisionDamping;

  // Handle horizontal boundaries
  if (pos.x - radius < 0) {
    pos.x = radius;
    vel.x *= -bounceFactor;
  } else if (pos.x + radius > params.screenWidth) {
    pos.x = params.screenWidth - radius;
    vel.x *= -bounceFactor;
  }

  // Handle vertical boundaries
  if (pos.y - radius < 0) {
    pos.y = radius;
    vel.y *= -bounceFactor;
  } else if (pos.y + radius > params.screenHeight) {
    pos.y = params.screenHeight - radius;
    vel.y *= -bounceFactor;
  }
}

std::pair<float, float>
Solver::calculateDensity(size_t i, Parameters params,
                         const std::vector<int> &neighbors) {
  float density = 0.0f;
  float nearDensity = 0.0f;
  const vec2 &targetPos = predictedPositions[i];

#pragma omp parallel for reduction(+ : density)
  for (size_t j = 0; j < neighbors.size(); ++j) {
    int neighborIdx = neighbors[j];
    if (i == neighborIdx)
      continue;

    float dist = Vector2Distance(predictedPositions[neighborIdx], targetPos);

    if (dist >= smoothingRadius)
      continue;

    density += params.mass * poly6Kernel(kernels, smoothingRadius, dist);
    nearDensity +=
        params.mass * spikePow3Kernel(kernels, smoothingRadius, dist);
  }

  // std::cout << density << std::endl;
  //     return fmax(density, params.targetDensity);
  return std::pair<float, float>(Clamp(density, 1.0e-15, 100.0),
                                 Clamp(nearDensity, 1.0e-15, 100.0));

  // return std::max(density, params.targetDensity);
}

Forces Solver::calculateForces(size_t i, Parameters params,
                               const std::vector<int> &neighbors) {
  vec2 tensionForce = {0.0f, 0.0f};
  vec2 pressureForce = {0.0f, 0.0f};
  vec2 nearPressureForce = {0.0f, 0.0f};
  vec2 viscosityForce = {0.0f, 0.0f};
  vec2 mouseForce = {0.0f, 0.0f};

  const vec2 &targetPos = predictedPositions[i];
  const float targetDensity = densities[i];

  for (int neighborIdx : neighbors) {
    if (i == neighborIdx)
      continue;

    vec2 offset = predictedPositions[neighborIdx] - targetPos;
    float dist = Vector2Length(offset);

    if (dist >= smoothingRadius)
      continue;

    // Handle particles at exactly the same position
    vec2 dir;
    if (dist == 0.0f) {
      dir = vec2{static_cast<float>(GetRandomValue(-1, 1)),
                 static_cast<float>(GetRandomValue(-1, 1))};
      dir = Vector2Normalize(dir);
    } else {
      dir = Vector2Scale(offset, (1.0f / dist));
    }

    float slope = spikeGradKernel(kernels, smoothingRadius, dist);
    float nearSlope = spikePow3GradKernel(kernels, smoothingRadius, dist);
    float neighborDensity = std::max(densities[neighborIdx], 1e-6f);

    float pressure = densityToPressure(targetDensity, params);
    float neighborPressure = densityToPressure(neighborDensity, params);
    float nearPressure = nearDensityToNearPressure(nearDensities[i], params);
    float nearNeighborPressure =
        nearDensityToNearPressure(nearDensities[neighborIdx], params);

    float sharedPressure = (pressure + neighborPressure) * 0.5;
    float nearSharedPressure = (nearPressure + nearNeighborPressure) * 0.5;

    pressureForce =
        pressureForce + Vector2Scale(dir, params.mass * sharedPressure * slope /
                                              neighborDensity);
    pressureForce =
        pressureForce +
        Vector2Scale(dir, params.mass * nearSharedPressure * nearSlope /
                              nearDensities[neighborIdx]);

    vec2 velocityDiff = velocities[neighborIdx] - velocities[i];
    float velocityAlongDir = Vector2DotProduct(velocityDiff, dir);
    viscosityForce =
        viscosityForce +
        Vector2Scale(velocityDiff,
                     params.mass * params.viscosity *
                         viscKernel(kernels, smoothingRadius, dist) /
                         neighborDensity);

    if (IsMouseButtonDown(MOUSE_LEFT_BUTTON) ||
        IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
      float mouseStrength = params.mouseStrength;
      if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
        mouseStrength = -params.mouseStrength;
      }
      vec2 mousePos = GetMousePosition();
      vec2 mouseDir = Vector2Subtract(mousePos, targetPos);
      float mouseDist = std::max(Vector2Length(mouseDir), 1e-3f);
      if (mouseDist < params.mouseRadius) {
        mouseForce =
            Vector2Add(mouseForce,
                       Vector2Scale(Vector2Normalize(mouseDir), mouseStrength));
      }
    }

    // Calculate tension force
  }
  // std::cout << pressureForce.x << pressureForce.y << std::endl;
  return Forces{tensionForce, pressureForce, viscosityForce, mouseForce};
}

void Solver::precomputeInteractions(Parameters params) { buildSpatialGrid(); }

void Solver::buildSpatialGrid() {
  spatialGrid.clear();

  // Pre-calculate grid size based on particle density
  size_t estimatedCells = (positions.size() / 4) + 1;
  spatialGrid.reserve(estimatedCells);

#pragma omp parallel
  {
    std::unordered_map<GridCell, std::vector<int>> localGrid;

#pragma omp for nowait
    for (size_t i = 0; i < positions.size(); ++i) {
      const vec2 &pos = positions[i];
      GridCell cell{static_cast<int>(pos.x / cellSize),
                    static_cast<int>(pos.y / cellSize)};
      localGrid[cell].push_back(i);
    }

#pragma omp critical
    {
      for (const auto &[cell, particles] : localGrid) {
        auto &targetCell = spatialGrid[cell];
        targetCell.insert(targetCell.end(), particles.begin(), particles.end());
      }
    }
  }
}

std::vector<int> Solver::getNeighbors(size_t index) {
  std::vector<int> neighbors;
  neighbors.reserve(40); // Typical number of neighbors

  const vec2 &pos = predictedPositions[index];
  int cellX = static_cast<int>(pos.x / cellSize);
  int cellY = static_cast<int>(pos.y / cellSize);

  // Search neighboring cells
  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      GridCell neighborCell{cellX + dx, cellY + dy};
      auto it = spatialGrid.find(neighborCell);
      if (it != spatialGrid.end()) {
        const auto &cellParticles = it->second;
        neighbors.insert(neighbors.end(), cellParticles.begin(),
                         cellParticles.end());
      }
    }
  }

  return neighbors;
}

Solver::~Solver() {
  // No need to delete global arrays as they are automatically managed by
  // std::vector
  spatialGrid.clear();
  interactionCache.clear();
  positions.clear();
  predictedPositions.clear();
  velocities.clear();
  densities.clear();
  colors.clear();
}
