#include <iostream>
#include <vector>

#include "physics.hpp"
#include "raylib-cpp.hpp"
#include "utils.hpp"

class Particle : public PhysObj {
public:
  Particle(vec2 pos) : PhysObj(pos) {
    radius = 20;
    color = RED;
  }
  float radius;
  Color color;
};

int main(void) {
  // test
  //  Initialization
  //---------------------------------------------------------
  const int screenWidth = 800;
  const int screenHeight = 600;

  rl::Window window(screenWidth, screenHeight,
                    "raylib [shapes] example - collision area");

  std::vector<Particle *> particles;
  for (int i = 0; i < 10; i++) {
    particles.push_back(
        new Particle(vec2{float(GetRandomValue(0, screenWidth)),
                          float(GetRandomValue(0, screenHeight))}));
  }
  // push cards to solver
  Solver solver;
  for (int i = 0; i < particles.size(); i++) {
    solver.objects.push_back(particles[i]);
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
      solver.update(dt);
    }

    window.ClearBackground(RAYWHITE);

    BeginDrawing();

    for (int i = 0; i < particles.size(); i++) {
      DrawCircleV(particles[i]->getPos(), particles[i]->radius,
                  particles[i]->color);
    }

    DrawFPS(10, 10);

    if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
      DrawText("Left mouse button is pressed", 100, 40, 10, RED);
    }

    EndDrawing();
  }

  return 0;
}
