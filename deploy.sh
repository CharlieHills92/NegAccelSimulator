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

SRC="$(cd "$(dirname "$0")" && pwd)"
SHARED_ROOT="/work/projects/nbsimulations/ibsimu"
DEST="${SHARED_ROOT}/ibsimu_MLsolver"

echo "============================================="
echo " Deploying ibsimu simulation"
echo " Source:      $SRC"
echo " Destination: $DEST"
echo "============================================="

# --- 1. Create directory structure ---
echo ""
echo "[1/6] Creating directory structure..."
mkdir -p "$DEST"/{build,geom,Bfield_MITICA,Bfield_MTF,Bfield_SPIDER,Bfields,densprofiles,EAMCCsecondaries,param,libibsimu_patched}

# --- 2. Copy source code, Makefile, and scripts ---
echo "[2/7] Copying source code..."
cp -v "$SRC"/*.cpp "$DEST/"
cp -v "$SRC"/*.h "$DEST/"
cp -v "$SRC"/Makefile "$DEST/"
cp -v "$SRC"/setup_environment.sh "$DEST/"
cp -v "$SRC"/.gitignore "$DEST/" 2>/dev/null || true
cp -v "$SRC"/deploy.sh "$DEST/" 2>/dev/null || true
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

# --- 5. Copy patched IBSimu source tree ---
echo "[5/6] Copying patched IBSimu source tree..."
cp -rv "$SRC"/libibsimu_patched/* "$DEST/libibsimu_patched/"

# --- 6. Set permissions: readable and executable by all, writable by group ---
echo "[6/6] Setting permissions..."
# Directories: rwxrwsr-x (2775) - setgid to preserve group
find "$DEST" -type d -exec chmod 2775 {} \;
# Source/header files: rw-rw-r-- (664)
find "$DEST" -type f \( -name "*.cpp" -o -name "*.h" -o -name "*.hpp" -o -name "*.scn" -o -name "*.py" -o -name "*.fld" -o -name "*.pc" -o -name "*.in" -o -name "*.am" -o -name "*.ac" \) -exec chmod 664 {} \;
# Makefile and scripts: rwxrwxr-x (775)
chmod 775 "$DEST/Makefile"
chmod 775 "$DEST/setup_environment.sh"
chmod 775 "$DEST/deploy.sh" 2>/dev/null || true
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
