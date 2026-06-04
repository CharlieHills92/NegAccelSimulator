#!/bin/bash
# Improved IBSimu environment configuration
# Source this file before compiling: source ./setup_environment.sh

echo "Setting up IBSimu environment..."

# Load required modules
echo "Loading modules..."
module load GSL 2>/dev/null || echo "Warning: GSL module not available"
module load GTK3 2>/dev/null || echo "Warning: GTK3 module not available"  
module load cairo 2>/dev/null || echo "Warning: cairo module not available"
module load matplotlib 2>/dev/null || echo "Warning: matplotlib module not available"

# IBSimu paths - use the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export IBSIMU_ROOT="$SCRIPT_DIR"
export IBSIMU_DIR="${IBSIMU_ROOT}/libibsimu_patched"

# Library paths
export LDFLAGS="${LDFLAGS} -L${IBSIMU_DIR}/src/.libs -Wl,-rpath,${IBSIMU_DIR}/src/.libs"
# export LD_LIBRARY_PATH="${IBSIMU_DIR}/src/.libs${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export LD_LIBRARY_PATH="$PWD/libibsimu_patched/src/.libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# PKG_CONFIG paths
export PKG_CONFIG_PATH="${IBSIMU_DIR}${PKG_CONFIG_PATH:+:${PKG_CONFIG_PATH}}"

# Add module pkg-config paths if available
if [ -n "$EBROOTGSL" ]; then
    export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${EBROOTGSL}/lib/pkgconfig"
    echo "Added GSL pkg-config path: ${EBROOTGSL}/lib/pkgconfig"
fi

if [ -n "$EBROOTGTK3" ]; then
    export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${EBROOTGTK3}/lib/pkgconfig"
    echo "Added GTK3 pkg-config path: ${EBROOTGTK3}/lib/pkgconfig"
fi

if [ -n "$EBROOTCAIRO" ]; then
    export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${EBROOTCAIRO}/lib/pkgconfig"
    echo "Added cairo pkg-config path: ${EBROOTCAIRO}/lib/pkgconfig"
fi

# SLURM queue format
export SQUEUE_FORMAT="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R"

# Test if pkg-config works
echo "Testing pkg-config setup..."
if pkg-config --exists ibsimu-1.0.6dev; then
    echo "✓ IBSimu pkg-config working"
    echo "  Includes: $(pkg-config --cflags ibsimu-1.0.6dev)"
else
    echo "⚠ IBSimu pkg-config not working, Makefile will use fallback"
fi

if pkg-config --exists gsl; then
    echo "✓ GSL pkg-config working"
elif [ -n "$EBROOTGSL" ]; then
    echo "⚠ GSL pkg-config not working, but module loaded at: $EBROOTGSL"
else
    echo "✗ GSL not available"
fi

if pkg-config --exists cairo; then
    echo "✓ cairo pkg-config working" 
elif [ -n "$EBROOTCAIRO" ]; then
    echo "⚠ cairo pkg-config not working, but module loaded at: $EBROOTCAIRO"
else
    echo "✗ cairo not available"
fi

if pkg-config --exists gtk+-3.0; then
    echo "✓ GTK3 pkg-config working"
elif [ -n "$EBROOTGTK3" ]; then
    echo "⚠ GTK3 pkg-config not working, but module loaded at: $EBROOTGTK3"
else
    echo "✗ GTK3 not available"
fi

echo "Environment setup complete!"
echo ""
echo "Usage:"
echo "  source ./setup_environment.sh  # Setup environment"
echo "  make                           # Build project"
echo "  make debug                     # Build with debug flags"
echo "  make clean                     # Clean build files"
