#!/bin/bash
# =============================================================================
# deploy.sh - Deploy ibsimu simulation to shared location
#
# Copies all files needed to build and run the simulation to:
#   /work/projects/nbsimulations/ibsimu/
#
# After deployment, any user in the sdcc-nb group can:
#   cd /work/projects/nbsimulations/ibsimu/
#   source setup_environment.sh
#   make
#   ./runtest_new_v2 <scan_tag>
# =============================================================================

set -euo pipefail

SRC="/home/ITER/poggic/Simulations/ibsimu/ibsimu_newsolver_20250825"
IBSIMU_SRC="/home/ITER/poggic/Simulations/ibsimu"
SHARED_ROOT="/work/projects/nbsimulations/ibsimu"
DEST="${SHARED_ROOT}/ibsimu_MLsolver"

echo "============================================="
echo " Deploying ibsimu simulation"
echo " Source:      $SRC"
echo " Destination: $DEST"
echo " Shared libs: $SHARED_ROOT"
echo "============================================="

# --- 1. Create directory structure ---
echo ""
echo "[1/7] Creating directory structure..."
mkdir -p "$DEST"/{build,geom,Bfield_MITICA,Bfield_MTF,Bfield_SPIDER,Bfields,densprofiles,EAMCCsecondaries,param}
mkdir -p "$SHARED_ROOT"/{lib,include}

# --- 2. Copy source code, Makefile, and scripts ---
echo "[2/7] Copying source code..."
cp -v "$SRC"/*.cpp "$DEST/"
cp -v "$SRC"/*.h "$DEST/"
cp -v "$SRC"/Makefile "$DEST/"
cp -v "$SRC"/setup_environment.sh "$DEST/"
cp -v "$SRC"/plot_trajectories_vtk.py "$DEST/" 2>/dev/null || true

# --- 3. Copy scan files ---
echo "[3/7] Copying scan files..."
cp -v "$SRC"/*.scn "$DEST/"

# --- 4. Copy data directories ---
echo "[4/7] Copying data directories..."
echo "  geom..."
cp -rv "$SRC"/geom/* "$DEST/geom/"
echo "  Bfield_MITICA..."
cp -rv "$SRC"/Bfield_MITICA/* "$DEST/Bfield_MITICA/"
echo "  Bfield_MTF..."
cp -rv "$SRC"/Bfield_MTF/* "$DEST/Bfield_MTF/"
echo "  Bfield_SPIDER..."
cp -rv "$SRC"/Bfield_SPIDER/* "$DEST/Bfield_SPIDER/" 2>/dev/null || echo "  (Bfield_SPIDER is empty)"
echo "  Bfields..."
cp -rv "$SRC"/Bfields/* "$DEST/Bfields/"
echo "  densprofiles..."
cp -rv "$SRC"/densprofiles/* "$DEST/densprofiles/"
echo "  EAMCCsecondaries..."
cp -rv "$SRC"/EAMCCsecondaries/* "$DEST/EAMCCsecondaries/"
echo "  param..."
cp -rv "$SRC"/param/* "$DEST/param/"

# --- 5. Copy ibsimu library and headers to shared location ---
echo "[5/7] Copying ibsimu library and headers to shared location..."
cp -rv "$IBSIMU_SRC"/lib/* "$SHARED_ROOT/lib/"
cp -rv "$IBSIMU_SRC"/include/* "$SHARED_ROOT/include/"

# --- 6. Patch all paths for new location ---
echo "[6/7] Patching paths for new location..."

# Patch the pkgconfig file to reference the shared location instead of personal home
sed -i "s|^prefix=.*|prefix=${SHARED_ROOT}|" "$SHARED_ROOT/lib/pkgconfig/ibsimu-1.0.6dev.pc"

# Patch the .la file if present
if [ -f "$SHARED_ROOT/lib/libibsimu-1.0.6dev.la" ]; then
    sed -i "s|/home/ITER/poggic/Simulations/ibsimu|${SHARED_ROOT}|g" "$SHARED_ROOT/lib/libibsimu-1.0.6dev.la"
fi

# Patch Makefile: update IBSIMU_ROOT to point to the shared location
sed -i "s|^IBSIMU_ROOT = .*|IBSIMU_ROOT = ${SHARED_ROOT}|" "$DEST/Makefile"

# Patch setup_environment.sh: update IBSIMU_ROOT
sed -i "s|^export IBSIMU_ROOT=.*|export IBSIMU_ROOT=\"${SHARED_ROOT}\"|" "$DEST/setup_environment.sh"

# --- 7. Set permissions: readable and executable by all, writable by group ---
echo "[7/7] Setting permissions..."
# Directories: rwxrwsr-x (2775) - setgid to preserve group
find "$DEST" -type d -exec chmod 2775 {} \;
find "$SHARED_ROOT/lib" -type d -exec chmod 2775 {} \;
find "$SHARED_ROOT/include" -type d -exec chmod 2775 {} \;
# Source/header files: rw-rw-r-- (664)
find "$DEST" -type f \( -name "*.cpp" -o -name "*.h" -o -name "*.scn" -o -name "*.py" -o -name "*.fld" \) -exec chmod 664 {} \;
# Makefile and scripts: rwxrwxr-x (775)
chmod 775 "$DEST/Makefile"
chmod 775 "$DEST/setup_environment.sh"
chmod 775 "$DEST/deploy.sh" 2>/dev/null || true
# Shared library files
find "$SHARED_ROOT/lib" -type f -exec chmod 664 {} \;
find "$SHARED_ROOT/lib" -name "*.so*" -exec chmod 775 {} \;
# Shared include files
find "$SHARED_ROOT/include" -type f -exec chmod 664 {} \;
# Data files: rw-rw-r-- (664)
find "$DEST"/{densprofiles,EAMCCsecondaries,param,geom,Bfield_MITICA,Bfield_MTF,Bfield_SPIDER,Bfields} -type f -exec chmod 664 {} \; 2>/dev/null || true
# Build directory: writable so users can compile
chmod 2775 "$DEST/build"

echo ""
echo "============================================="
echo " Deployment complete!"
echo ""
echo " Total size:"
du -sh "$DEST"
echo ""
echo " To build and run (any user):"
echo "   cd $DEST"
echo "   source setup_environment.sh"
echo "   make clean && make"
echo "   ./runtest_new_v2 <scan_tag>"
echo "============================================="
