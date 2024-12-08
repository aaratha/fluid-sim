#include "physics.hpp"
#include <raylib.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include "utils.hpp"

const float EPS = 1e-6f;

Obstacle::Obstacle(vec2 position, int shape, float radius, Rectangle rectangle)
    : position(position), shape(shape), radius(radius), rectangle(rectangle) {}

void Obstacle::lerpRadius(float target, float transition_time, float dt) {
    // Calculate rate in seconds to complete transition
    float increment = dt / transition_time;
    radius = lerp1D(radius, target, increment);
}

float Obstacle::getRadius() {
    return radius;
}

void Obstacle::setRadius(float r) {
    radius = r;
}

vec2 Obstacle::getPos() {
    return position;
};

void Obstacle::setPos(vec2 pos) {
    position = pos;
};

void Solver::initializeCache(size_t particleCount) {
    interactionCache.resize(particleCount,
                            std::vector<float>(particleCount, 0.0f));
    positions.resize(particleCount);
    predictedPositions.resize(particleCount);
    velocities.resize(particleCount, vec2{0.0f, 0.0f});
    densities.resize(particleCount, 0.0f);
    nearDensities.resize(particleCount, 0.0f);
    colors.resize(particleCount);
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
        vec2 acc = (forces.pressure + forces.tension + forces.viscosity +
                    forces.mouse) /
                   densities[i];  // / std::max(densities[i], 1e-6f);
        float acc_limit = params.maxAcceleration;
        float acc_mag = Vector2Length(acc);
        if (acc_mag > acc_limit * acc_limit) {
            acc = acc / sqrt(acc_mag) * acc_limit;
        }
        velocities[i] += acc * dt;

        // Apply velocity constraints
        float currentSpeed = Vector2Length(velocities[i]);
        if (currentSpeed > params.maxVelocity) {
            velocities[i] = Vector2Scale(Vector2Normalize(velocities[i]),
                                         params.maxVelocity);
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
    vec2& pos = positions[i];
    vec2& predictedPos = predictedPositions[i];
    vec2& vel = velocities[i];
    float radius = params.particleRadius;
    float bounceFactor = params.collisionDamping;
    float collisionStiffness = 0.5f;  // Adjustable stiffness parameter

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

    // Handle obstacle collisions
    for (auto obstacle : obstacles) {
        if (obstacle.shape == EXT_CIRCLE || obstacle.shape == INT_CIRCLE) {
            vec2 obstaclePos = obstacle.position;
            float obstacleRadius = obstacle.radius;

            vec2 dir = Vector2Subtract(pos, obstaclePos);
            float dist = Vector2Length(dir);
            float overlap = radius + obstacleRadius - dist;

            if (obstacle.shape == EXT_CIRCLE) {
                // External circle: overlap occurs when distance is less than
                // sum of radii
                if (overlap > 0) {
                    vec2 normal = Vector2Normalize(dir);
                    pos += Vector2Scale(normal, overlap);

                    // Reflect velocity
                    float dot = Vector2DotProduct(vel, normal);
                    if (dot < 0) {
                        vel -= Vector2Scale(normal, 2.0f * dot);
                    }
                }
            } else if (obstacle.shape == INT_CIRCLE) {
                if (overlap < 0) {
                    vec2 normal = Vector2Normalize(dir);
                    pos += Vector2Scale(normal, overlap);

                    // Reflect velocity
                    float dot = Vector2DotProduct(vel, normal);
                    if (dot < 0) {
                        vel -= Vector2Scale(normal, 2.0f * dot);
                    }
                }
            } else if (obstacle.shape == EXT_RECTANGLE) {
                Rectangle rect = obstacle.rectangle;

                vec2 closest = {
                    std::clamp(pos.x, rect.x, rect.x + rect.width),
                    std::clamp(pos.y, rect.y, rect.y + rect.height)};

                vec2 dir = Vector2Subtract(pos, closest);
                float dist = Vector2Length(dir);

                if (dist < radius) {
                    vec2 normal = Vector2Normalize(dir);
                    float overlap = radius - dist;

                    // Soft constraint resolution
                    vec2 correction = normal * overlap * collisionStiffness;
                    pos += correction;

                    // Velocity modification with dampening
                    vel += correction * (1.0f / GetFrameTime());
                    vel *= 0.9f;  // Soft velocity reduction
                }
            }
        }
    }
}

std::pair<float, float> Solver::calculateDensity(
    size_t i,
    Parameters params,
    const std::vector<int>& neighbors) {
    float density = 0.0f;
    float nearDensity = 0.0f;
    const vec2& targetPos = predictedPositions[i];

#pragma omp parallel for reduction(+ : density)
    for (size_t j = 0; j < neighbors.size(); ++j) {
        int neighborIdx = neighbors[j];
        if (i == neighborIdx)
            continue;

        float dist =
            Vector2Distance(predictedPositions[neighborIdx], targetPos);

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

Forces Solver::calculateForces(size_t i,
                               Parameters params,
                               const std::vector<int>& neighbors) {
    vec2 tensionForce = {0.0f, 0.0f};
    vec2 pressureForce = {0.0f, 0.0f};
    vec2 nearPressureForce = {0.0f, 0.0f};
    vec2 viscosityForce = {0.0f, 0.0f};
    vec2 mouseForce = {0.0f, 0.0f};

    const vec2& targetPos = predictedPositions[i];
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
        float nearPressure =
            nearDensityToNearPressure(nearDensities[i], params);
        float nearNeighborPressure =
            nearDensityToNearPressure(nearDensities[neighborIdx], params);

        float sharedPressure = (pressure + neighborPressure) * 0.5;
        float nearSharedPressure = (nearPressure + nearNeighborPressure) * 0.5;

        pressureForce =
            pressureForce + Vector2Scale(dir, params.mass * sharedPressure *
                                                  slope / neighborDensity);
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
                mouseForce = Vector2Add(
                    mouseForce,
                    Vector2Scale(Vector2Normalize(mouseDir), mouseStrength));
            }
        }

        // Calculate tension force
    }
    // std::cout << pressureForce.x << pressureForce.y << std::endl;
    return Forces{tensionForce, pressureForce, viscosityForce, mouseForce};
}

void Solver::precomputeInteractions(Parameters params) {
    buildSpatialGrid();
}

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
            const vec2& pos = positions[i];
            GridCell cell{static_cast<int>(pos.x / cellSize),
                          static_cast<int>(pos.y / cellSize)};
            localGrid[cell].push_back(i);
        }

#pragma omp critical
        {
            for (const auto& [cell, particles] : localGrid) {
                auto& targetCell = spatialGrid[cell];
                targetCell.insert(targetCell.end(), particles.begin(),
                                  particles.end());
            }
        }
    }
}

std::vector<int> Solver::getNeighbors(size_t index) {
    std::vector<int> neighbors;
    neighbors.reserve(40);  // Typical number of neighbors

    const vec2& pos = predictedPositions[index];
    int cellX = static_cast<int>(pos.x / cellSize);
    int cellY = static_cast<int>(pos.y / cellSize);

    // Search neighboring cells
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            GridCell neighborCell{cellX + dx, cellY + dy};
            auto it = spatialGrid.find(neighborCell);
            if (it != spatialGrid.end()) {
                const auto& cellParticles = it->second;
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
