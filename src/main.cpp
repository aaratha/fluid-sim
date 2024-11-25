#include <iostream>
#include <string>
#include <vector>

#include "physics.hpp"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#include "raylib-cpp.hpp"
#include "utils.hpp"

Parameters params = {.screenWidth = 1280,
                     .screenHeight = 720,
                     .particleCount = 2000,
                     .particleRadius = 3.0f,
                     .collisionDamping = 1.0f,
                     .friction = 1.0f,
                     .gravity = 0.0,
                     .smoothingMultiplier = 20.0f,
                     .substeps = 8,
                     .targetDensity = 2.5f,
                     .pressureMultiplier = 2.0f,
                     .maxVelocity = 100.0f,
                     .nearPressureMultiplier = 300.0};

int main(void) {
  // Initialization
  //---------------------------------------------------------
  const int screenWidth = params.screenWidth;
  const int screenHeight = params.screenHeight;

  rl::Window window(screenWidth, screenHeight,
                    "raylib [shapes] example - collision area");

  // Initialize solver and particles
  Solver solver;
  solver.positions.resize(params.particleCount);
  solver.predictedPositions.resize(params.particleCount);
  solver.velocities.resize(params.particleCount, vec2{0.0f, 0.0f});
  solver.densities.resize(params.particleCount, 0.0f);
  solver.nearDensities.resize(params.particleCount, 0.0f);
  solver.colors.resize(params.particleCount);

  // Initialize particle positions randomly
  for (int i = 0; i < params.particleCount; ++i) {
    vec2 randPos =
        vec2{float(GetRandomValue(screenWidth / 3, screenWidth * 2 / 3)),
             float(GetRandomValue(screenHeight / 3, screenHeight * 2 / 3))};
    solver.positions[i] = randPos;
    solver.predictedPositions[i] = randPos;
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
                 &params.pressureMultiplier, 0.1, 2000);

    std::string nearMultiplierToString =
        "near pressure multiplier: " +
        std::to_string(params.nearPressureMultiplier);
    const char *nearMultiplier = nearMultiplierToString.c_str();
    GuiSliderBar((Rectangle){10, 100, 120, 20}, NULL, nearMultiplier,
                 &params.nearPressureMultiplier, -19, 500);

    std::string densityToString =
        "target density: " + std::to_string(params.targetDensity);
    const char *density = densityToString.c_str();
    GuiSliderBar((Rectangle){10, 130, 120, 20}, NULL, density,
                 &params.targetDensity, 0.1, 1000);

    std::string radiusToString =
        "radius: " + std::to_string(params.particleRadius);
    const char *radius = radiusToString.c_str();
    GuiSliderBar((Rectangle){10, 160, 120, 20}, NULL, radius,
                 &params.particleRadius, 0, 20);

    std::string gravityToString = "gravity: " + std::to_string(params.gravity);
    const char *gravity = gravityToString.c_str();
    GuiSliderBar((Rectangle){10, 190, 120, 20}, NULL, gravity, &params.gravity,
                 -100, 100);

    EndDrawing();
  }

  return 0;
}
