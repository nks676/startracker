# CMake toolchain file: cross-compile for aarch64-linux-gnu (Raspberry Pi 4/5).
#
# Usage:
#   cmake -B build-arm64 \
#         -DCMAKE_TOOLCHAIN_FILE=cmake/aarch64-linux-gnu.cmake \
#         -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-arm64 --parallel
#
# Requires the host to have the cross toolchain installed, e.g. on
# Debian/Ubuntu:  sudo apt-get install g++-aarch64-linux-gnu
#
# This toolchain is intentionally minimal: it targets the same Linux ABI
# as the Pi's stock Raspberry Pi OS / Ubuntu Server install. No SIMD or
# CPU-specific tuning is set here — that lives in the main CMakeLists.

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR aarch64)

set(CMAKE_C_COMPILER   aarch64-linux-gnu-gcc)
set(CMAKE_CXX_COMPILER aarch64-linux-gnu-g++)

# Optional: a sysroot can be supplied at configure time via
# -DCMAKE_SYSROOT=/path/to/sysroot. We don't hardcode one because the
# Debian cross-gcc package ships its own multiarch sysroot under
# /usr/aarch64-linux-gnu and CMake picks it up automatically.

# When searching for libraries / headers / packages, prefer the target
# environment (cross sysroot) over the host. Programs (compilers etc.)
# must still come from the host.
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
