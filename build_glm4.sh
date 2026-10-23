#!/bin/bash
# Build GLM4 + AED3 from the pinned upstream source, identically on the local
# machine and on VUB Hydra.
#
#   bash build_glm4.sh
#
# Source: upstream/glm-source from AquaticEcoDynamics/glm-aed, checked out at
# the commit recorded in upstream/ (git -C upstream rev-parse HEAD), with the
# submodule revisions that commit pins (git -C upstream submodule status).
#
# Toolchain, chosen so both machines match exactly:
#   local : Ubuntu 24.04, system gfortran/gcc 13.3.0, system netCDF 4.9.2
#   Hydra : module GCC/13.3.0 + netCDF/4.9.2-gompi-2024a
# Built with --no-gui: no plotting window/libgd, which a batch cluster doesn't
# need. The resulting binary is bin/glm-built; build details go to
# bin/glm-built.info.
#
# Why build at all: the prebuilt GLM 4.0.0 binary needs glibc 2.38 and Hydra
# (Rocky Linux 9) has 2.34, so it can't run there.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/upstream/glm-source"

if command -v module >/dev/null 2>&1 && [ -n "${VSC_INSTITUTE:-}" ]; then
  module purge
  module load GCC/13.3.0 netCDF/4.9.2-gompi-2024a
  WHERE="hydra"
else
  # Make sure the system netCDF is used, not a conda one that may come first
  # in PATH (its libraries wouldn't be on the runtime library path).
  export PATH="/usr/bin:$PATH"
  WHERE="local"
fi
# build_glm.sh overwrites LIBRARY_PATH (export LIBRARY_PATH=$LIB), which
# discards the path a netCDF environment module sets -- the link then fails
# with "cannot find -lnetcdf". GLM's Makefile adds -L$(NETCDFLIB) when set.
export NETCDFLIB="$(nc-config --libdir)"
echo "Building on: $WHERE"
echo "NETCDFLIB  : $NETCDFLIB"
echo "nc-config  : $(command -v nc-config) ($(nc-config --version))"
echo "gfortran   : $(gfortran --version | head -1)"

cd "$SRC"
# Local patches (patches/<submodule>_*.patch, applied to that submodule of glm-source):
#   libaed-api_private_extc.patch -- with bioshade_feedback = .false., AED computes its own
#     light with self-shading on a private copy of the extinction coefficients, so GLM's
#     heating is untouched (upstream: AED always writes its extinction into GLM's array).
#   GLM_ice_turbulent_wind.patch -- over ice, sensible (and, while melting, latent) heat used the
#     ~zero wind meant for under-ice mixing (WindSp = 0.0001), against the code's own comment; uses the
#     real wind. Without it high-Arctic lakes never melt out (global test, 2026-10-22).
#   GLM_aed_light_below_ice.patch -- AED's surface irradiance (I_0, used by AED's own light) was the raw
#     daily-mean forcing shortwave (no albedo, no daily cycle, no ice/snow); uses the shortwave
#     entering the water from GLM. Without it phytoplankton bloomed under ice in full sun.
#   GLM_daily_shortwave.patch -- daily forcing (subdaily = .false.) spreads the day's shortwave as a
#     half-sine over the day length. Three bugs fixed: (1) declination sign (-23.45 sin(...) gave
#     northern lakes winter day lengths in summer, polar night in June at 74N); (2) acos argument
#     clamped to [-1, 1] (polar day/night gave NaN -> no sunlight); (3) half-sine centred on
#     solar noon (it was centred on clock 0 = midnight: light 20:00-05:00 in hourly output).
# Each patched submodule is first reset to its pinned commit (tracked files only), then all its
# patches are applied in name order; any failure stops the build. So the binary is always
# "pinned source + exactly these patches", on a fresh checkout or a reused one.
for sub in $(ls "$HERE"/patches/*.patch 2>/dev/null | xargs -n1 basename | cut -d_ -f1 | sort -u); do
  git -C "$SRC/$sub" checkout -- .
  for p in "$HERE"/patches/"${sub}"_*.patch; do
    git -C "$SRC/$sub" apply "$p"
    echo "patch applied: $(basename "$p")"
  done
done
# Always build from a clean tree, so nothing compiled on another machine (or
# with another toolchain) gets linked in. Upstream's clean.sh skips
# libaed-api, so clean that one explicitly.
./clean.sh > /dev/null 2>&1
(cd libaed-api && make distclean > /dev/null 2>&1) || true
./build_glm.sh --no-gui --dont-package 2>&1 | tee "$HERE/build_$WHERE.log"

# install with an atomic rename: GLM runs already started keep the old file, new runs get the new one
cp "$SRC/GLM/glm" "$HERE/bin/glm-built.new" && mv -f "$HERE/bin/glm-built.new" "$HERE/bin/glm-built"
{
  echo "built on        : $WHERE ($(hostname)) $(date -u +%FT%TZ)"
  echo "glm-aed commit  : $(git -C "$HERE/upstream" rev-parse HEAD)"
  git -C "$HERE/upstream" submodule status --recursive | sed 's/^/  /'
  echo "local patches   :"
  (cd "$HERE/patches" && sha256sum *.patch 2>/dev/null | sed 's/^/  /')
  echo "gfortran        : $(gfortran --version | head -1)"
  echo "netCDF          : $(nc-config --version)"
  echo "glibc           : $(ldd --version | head -1)"
  echo "flags           : --no-gui --dont-package"
} > "$HERE/bin/glm-built.info"
cat "$HERE/bin/glm-built.info"
"$HERE/bin/glm-built" --version 2>&1 | grep -E "Version|API|AED version|build date" | sed 's/^ *//'
