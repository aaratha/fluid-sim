#!/usr/bin/env sh

# Check for the OS type

if [ "$1" = "web" ]; then
    emcmake cmake . -B build-web -DPLATFORM=Web -DCMAKE_BUILD_TYPE=Release
    emmake make -C build-web
    python -m http.server
else
    if [ "$OSTYPE" = "msys" ] || [ "$OSTYPE" = "win32" ]; then
        # Windows with MinGW
        cmake . -B build -G "MinGW Makefiles"
    else
        # Unix-based system (macOS/Linux)
        cmake . -B build -G "Unix Makefiles"
    fi
    cmake --build build
    build/fluid-sim
fi
