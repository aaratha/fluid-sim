#!/usr/bin/env sh

# Check for the OS type
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
    # Windows with MinGW
    cmake . -B build -G "MinGW Makefiles"
else
    # Unix-based system (macOS/Linux)
    cmake . -B build -G "Unix Makefiles"
fi
cmake --build build
build/art-game
