#!/usr/bin/env sh

# Check for the OS type and platform
if [ "$1" = "web" ]; then
    echo "Building for Web (Emscripten)..."
    emcmake cmake . -B build-web -DPLATFORM=Web -DCMAKE_BUILD_TYPE=Release
    emmake make -C build-web
    echo "Starting a local HTTP server for the Web build..."
    python3 -m http.server --directory build-web
else
    echo "Building for Desktop..."
    if [ "$(uname)" = "Darwin" ]; then
        # macOS
        cmake . -B build -G "Unix Makefiles" -DPLATFORM=Desktop
    elif [ "$(uname)" = "Linux" ]; then
        # Linux
        cmake . -B build -G "Unix Makefiles" -DPLATFORM=Desktop
    elif echo "$OSTYPE" | grep -q "msys"; then
        # Windows with MinGW
        cmake . -B build -G "MinGW Makefiles" -DPLATFORM=Desktop
    else
        echo "Unsupported OS: $OSTYPE"
        exit 1
    fi

    cmake --build build
    echo "Running the executable..."
    ./build/fluid-sim
fi
