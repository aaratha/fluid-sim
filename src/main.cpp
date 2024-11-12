#include <iostream>
#include <vector>

#include "physics.hpp"
#include "raylib-cpp.hpp"
#include "utils.hpp"

Parameters params = {.screenWidth = 800,
                     .screenHeight = 600,
                     .particleCount = 200,
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

    if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
      DrawText("Left mouse button is pressed", 100, 40, 10, RED);
    }

    EndDrawing();
  }

  return 0;
}
