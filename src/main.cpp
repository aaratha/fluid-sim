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
                     .particleCount = 300,
                     .particleRadius = 5.0f,
                     .collisionDamping = 0.8f,
                     .friction = 0.98f,
                     .gravity = vec2(0, 1000)};

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

  SetTargetFPS(60); // Set our game to run at 60 frames-per-second
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
      DrawCircleV(solver.objects[i]->getPos(), solver.objects[i]->getRadius(),
                  solver.objects[i]->color);
    }

    DrawFPS(10, 10);

    std::string gravToString = "Grav: " + std::to_string(params.gravity.y);
    const char *grav = gravToString.c_str();
    GuiSliderBar((Rectangle){10, 40, 120, 20}, NULL, grav, &params.gravity.y, 1,
                 1000);

    std::string bounceToString =
        "bounce: " + std::to_string(params.collisionDamping);
    const char *bounce = bounceToString.c_str();
    GuiSliderBar((Rectangle){10, 70, 120, 20}, NULL, bounce,
                 &params.collisionDamping, 0, 1);

    std::string radiusToString =
        "radius" + std::to_string(params.particleRadius);
    const char *radius = radiusToString.c_str();
    GuiSliderBar((Rectangle){10, 100, 120, 20}, NULL, radius,
                 &params.particleRadius, 0, 20);

    EndDrawing();
  }

  return 0;
}
