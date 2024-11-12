#include <iostream>
#include <vector>

#include "physics.hpp"
#include "raylib-cpp.hpp"
#include "utils.hpp"

const int PARTICLE_COUNT = 100;

int main(void) {
  // test
  //  Initialization
  //---------------------------------------------------------
  const int screenWidth = 800;
  const int screenHeight = 600;

  rl::Window window(screenWidth, screenHeight,
                    "raylib [shapes] example - collision area");

  // push cards to solver
  Solver solver;
  for (int i = 0; i < PARTICLE_COUNT; i++) {
    solver.objects.push_back(
        new PhysObj(vec2{float(GetRandomValue(0, screenWidth)),
                         float(GetRandomValue(0, screenHeight))}));
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
      solver.update(dt, screenWidth, screenHeight);
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
