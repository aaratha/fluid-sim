#include <iostream>
#include <string>
#include <vector>

#include "physics.hpp"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#include "raylib-cpp.hpp"
#include "utils.hpp"

Parameters params = {.screenWidth = 800,
                     .screenHeight = 600,
                     .particleCount = 1000,
                     .particleRadius = 5.0f,
                     .collisionDamping = 1.0f,
                     .friction = 0.9f,
                     .gravity = vec2(0, 1000),
                     .smoothingRadius = 120.0f,
                     .substeps = 8,
                     .targetDensity = 2.75f,
                     .pressureMultiplier = 20.0f};

int main(void) {
  // test
  //  Initialization
  //---------------------------------------------------------
  const int screenWidth = params.screenWidth;
  const int screenHeight = params.screenHeight;

  rl::Window window(screenWidth, screenHeight,
                    "raylib [shapes] example - collision area");

  // push cards to solver
  Solver solver;
  for (int i = 0; i < params.particleCount; i++) {
    solver.objects.push_back(
        new PhysObj(vec2{float(GetRandomValue(0, screenWidth)),
                         float(GetRandomValue(0, screenHeight))},
                    params.particleRadius));
  }

  bool pause = false; // Movement pause

  SetTargetFPS(120); // Set our game to run at 60 frames-per-second
  //----------------------------------------------------------

  // Main game loop
  while (!window.ShouldClose()) { // Detect window close button or ESC key
    // Update
    //-----------------------------------------------------
    // Move box if not paused

    if (IsKeyPressed(KEY_SPACE))
      pause = !pause;

    if (!pause) {
      float dt = GetFrameTime();
      solver.update(dt, params);
    }

    window.ClearBackground(BLACK);

    BeginDrawing();

    for (int i = 0; i < solver.objects.size(); i++) {
      DrawCircleV(solver.objects[i]->getPos(), params.particleRadius,
                  solver.objects[i]->color);
    }

    DrawFPS(10, 10);

    std::string smoothingToString =
        "smoothing radius" + std::to_string(params.smoothingRadius);
    const char *smoothing = smoothingToString.c_str();
    GuiSliderBar((Rectangle){10, 40, 120, 20}, NULL, smoothing,
                 &params.smoothingRadius, params.particleRadius, 1000);

    std::string multiplierToString =
        "pressure multiplier" + std::to_string(params.pressureMultiplier);
    const char *multiplier = multiplierToString.c_str();
    GuiSliderBar((Rectangle){10, 70, 120, 20}, NULL, multiplier,
                 &params.pressureMultiplier, 0.1, 1000);

    std::string densityToString =
        "target density" + std::to_string(params.targetDensity);
    const char *density = densityToString.c_str();
    GuiSliderBar((Rectangle){10, 100, 120, 20}, NULL, density,
                 &params.targetDensity, 0.1, 10);

    std::string radiusToString =
        "radius" + std::to_string(params.particleRadius);
    const char *radius = radiusToString.c_str();
    GuiSliderBar((Rectangle){10, 130, 120, 20}, NULL, radius,
                 &params.particleRadius, 0, 20);

    EndDrawing();
  }

  solver.~Solver();

  return 0;
}
