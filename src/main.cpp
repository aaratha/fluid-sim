#include <iostream>
#include <string>
#include <vector>

#include "physics.hpp"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#include "raylib-cpp.hpp"
#include "utils.hpp"

Parameters params = {.screenWidth = 500,
                     .screenHeight = 500,
                     .particleCount = 2000,
                     .particleRadius = 3.0f,
                     .collisionDamping = 0.95f,
                     .friction = 1.0f,
                     .gravity = 0.0,
                     .smoothingMultiplier = 13.0f,
                     .substeps = 8,
                     .targetDensity = 2.0f,
                     .pressureMultiplier = 100000.0f,
                     .maxVelocity = 400.0f,
                     .nearPressureMultiplier = 300.0,
                     .viscosity = 50.0,
                     .maxAcceleration = 200.0f,
                     .mass = 100000.0f};

int main(void) {
  // Initialization
  //---------------------------------------------------------
  const int screenWidth = params.screenWidth;
  const int screenHeight = params.screenHeight;

  rl::Window window(screenWidth, screenHeight, "SPH Fluid Simulation");

  // Initialize solver and particles
  Solver solver;
  solver.positions.resize(params.particleCount);
  solver.predictedPositions.resize(params.particleCount);
  solver.velocities.resize(params.particleCount, vec2{0.0f, 0.0f});
  solver.densities.resize(params.particleCount, 0.0f);
  solver.nearDensities.resize(params.particleCount, 0.0f);
  solver.colors.resize(params.particleCount);

  // Initialize particle positions in a grid
  int gridCols = static_cast<int>(sqrt(params.particleCount));
  int gridRows = params.particleCount / gridCols +
                 (params.particleCount % gridCols > 0 ? 1 : 0);

  float gridSpacing =
      params.particleRadius * 2.0f; // Spacing based on particle radius
  float startX =
      (screenWidth - gridCols * gridSpacing) / 2.0f; // Center horizontally
  float startY =
      (screenHeight - gridRows * gridSpacing) / 2.0f; // Center vertically

  for (int i = 0; i < params.particleCount; ++i) {
    int row = i / gridCols;
    int col = i % gridCols;

    vec2 gridPos = vec2{startX + col * gridSpacing, startY + row * gridSpacing};
    solver.positions[i] = gridPos;
    solver.predictedPositions[i] = gridPos;
  }

  bool pause = false; // Movement pause

  SetTargetFPS(120); // Set our game to run at 120 frames-per-second
  //----------------------------------------------------------

  // Main game loop
  while (!window.ShouldClose()) { // Detect window close button or ESC key
    // Update
    //-----------------------------------------------------
    // Pause handling
    if (IsKeyPressed(KEY_SPACE)) {
      pause = !pause;
    }

    if (!pause) {
      float dt = GetFrameTime();
      solver.update(dt, params);
    }

    // Rendering
    //-----------------------------------------------------
    window.ClearBackground(BLACK);

    BeginDrawing();

    for (int i = 0; i < params.particleCount; i++) {
      DrawCircleV(solver.positions[i], params.particleRadius, solver.colors[i]);
    }

    DrawCircle(GetMouseX(), GetMouseY(),
               params.smoothingMultiplier * params.particleRadius,
               (Color){0, 255, 0, 100});

    DrawFPS(10, 10);

    std::string smoothingToString =
        "smoothing radius: " + std::to_string(params.smoothingMultiplier);
    const char *smoothing = smoothingToString.c_str();
    GuiSliderBar((Rectangle){10, 40, 120, 20}, NULL, smoothing,
                 &params.smoothingMultiplier, 1, 50);

    std::string multiplierToString =
        "pressure multiplier: " + std::to_string(params.pressureMultiplier);
    const char *multiplier = multiplierToString.c_str();
    GuiSliderBar((Rectangle){10, 70, 120, 20}, NULL, multiplier,
                 &params.pressureMultiplier, 0.1, 100000);

    std::string accelerationToString =
        "max acceleration" + std::to_string(params.maxAcceleration);
    const char *acceleration = accelerationToString.c_str();
    GuiSliderBar((Rectangle){10, 100, 120, 20}, NULL, acceleration,
                 &params.maxAcceleration, 0, 100);

    std::string densityToString =
        "target density: " + std::to_string(params.targetDensity);
    const char *density = densityToString.c_str();
    GuiSliderBar((Rectangle){10, 130, 120, 20}, NULL, density,
                 &params.targetDensity, 0.0, 10);

    std::string radiusToString =
        "radius: " + std::to_string(params.particleRadius);
    const char *radius = radiusToString.c_str();
    GuiSliderBar((Rectangle){10, 160, 120, 20}, NULL, radius,
                 &params.particleRadius, 0, 20);

    std::string viscosityToString =
        "viscosity: " + std::to_string(params.viscosity);
    const char *viscosity = viscosityToString.c_str();
    GuiSliderBar((Rectangle){10, 190, 120, 20}, NULL, viscosity,
                 &params.viscosity, 0, 300);

    std::string gravityToString = "gravity: " + std::to_string(params.gravity);
    const char *gravity = gravityToString.c_str();
    GuiSliderBar((Rectangle){10, 220, 120, 20}, NULL, gravity, &params.gravity,
                 -1000, 1000);

    EndDrawing();
  }

  return 0;
}
