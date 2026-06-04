nsimp_plasmacyl_SOURCES = nsimp_plasmacyl.cpp
nsimp_plasma2d_SOURCES = nsimp_plasma2d.cpp
nsimp_plasma3d_SOURCES = nsimp_plasma3d.cpp
particles_bfield2d_SOURCES = particles_bfield2d.cpp
particles_bfieldcyl_SOURCES = particles_bfieldcyl.cpp
particles_cyl2_SOURCES = particles_cyl2.cpp
particles_cyl_SOURCES = particles_cyl.cpp
particles_diagnostic_SOURCES = particles_diagnostic.cpp
particles_edip_SOURCES = particles_edip.cpp
particles_mirror_SOURCES = particles_mirror.cpp
particles_relativistic_SOURCES = particles_relativistic.cpp
particles_trajectory_SOURCES = particles_trajectory.cpp
plasma2d_SOURCES = plasma2d.cpp
plasma3d_SOURCES = plasma3d.cpp
plasmacyl_SOURCES = plasmacyl.cpp
solver1d_scharge_SOURCES = solver1d_scharge.cpp
solver2d_coax_SOURCES = solver2d_coax.cpp
solver3d_sphere_SOURCES = solver3d_sphere.cpp
solvercyl_coax_SOURCES = solvercyl_coax.cpp

check_PROGRAMS = \
nsimp_plasmacyl \
nsimp_plasma2d \
nsimp_plasma3d \
particles_bfield2d \
particles_bfieldcyl \
particles_cyl2 \
particles_cyl \
particles_diagnostic \
particles_edip \
particles_mirror \
particles_relativistic \
particles_trajectory \
plasma2d \
plasma3d \
plasmacyl \
solver1d_scharge \
solver2d_coax \
solver3d_sphere \
solvercyl_coax

TESTS = $(check_PROGRAMS)

EXTRA_LTLIBRARIES = libibsimutest.la
libibsimutest_la_SOURCES = \
ibsimutest.cpp \
ibsimutest.hpp

AM_CPPFLAGS = $(CPPFLAGS) -Wall -g -I$(top_srcdir)/src
AM_LDFLAGS = $(LIBS) -L$(top_srcdir)/src @rdynamic_flag@
LDADD = $(top_builddir)/src/libibsimu-1.0.6dev.la libibsimutest.la -lm

if CSG
AM_CPPFLAGS += `pkg-config --cflags csg`
AM_LDFLAGS += `pkg-config --libs csg`
endif

EXTRA_DIST = \
solver2d_coax.gpl

CLEANFILES = *~ *.png *.eps *.svg *.pdf *.ps *.dat *_vout.txt plasmacyl_emit1.txt beam3d.path libibsimutest.la

