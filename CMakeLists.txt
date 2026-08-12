# ACES_PROGRAM superbuild.
#
# Roadmap step 3, phase 2: replaces the thin bash orchestrator from phase 1
# (configure/build.sh/run_tests.sh, commit bb3194e4) with a real CMake
# superbuild of aces2 (ACES II), aces3 (ACES III), and pyaces.
#
# This does NOT reimplement any subsystem's own compilation/physics logic:
#   aces2  -- still builds via its own Makefiles/GNUmakefile (gmake)
#   aces3  -- still builds via its own autotools (configure + gmake)
#   pyaces -- still builds via its own build.sh (f2py + objcopy surgery)
# CMake's job here is only to drive each one with the shared toolchain
# (Intel compilers + Intel MKL -- this project's already-settled choice,
# see ../aces3/configure_line.ufhpc and pyaces/build.sh) and wire them into
# one dependency graph (aces2 must exist before pyaces links against it;
# aces3 is independent) and one `ctest` entry point for the fast suites.
#
# Modeled on ~/Develop/aces4/CMakeLists.txt for house style (option(),
# cache variables, a cmake/ module directory) but deliberately NOT on its
# architecture: aces4 was a from-scratch, single-codebase CMake build that
# compiled everything itself via add_library/add_executable. Doing that
# here would mean rewriting three separate legacy Fortran/C/C++ build
# systems (one of them, pyaces/build.sh, encoding over a dozen hard-won
# symbol-collision fixes) in CMake's own language -- exactly the kind of
# from-scratch rewrite this two-phase plan exists to avoid repeating (see
# NOTES.md). So there is no enable_language()/find_package(MPI)/
# find_package(LAPACK) here: none of these targets compile anything
# directly under CMake's own compiler detection, so there is nothing for
# those checks to usefully gate at `cmake` configure time. Each subsystem's
# own build system still does its own compiler/library discovery, exactly
# as before -- CMake's ExternalProject_Add targets below just invoke it.
cmake_minimum_required(VERSION 3.18)
project(ACES_PROGRAM NONE)

include(ExternalProject)

# --- shared toolchain -------------------------------------------------
# aces2/aces3/pyaces each `module load` this same set before building or
# testing (mirrors phase-1's ./configure --modules). Not validated at
# `cmake` configure time -- like phase 1, that's left to each build step's
# own `module load`, since modules aren't necessarily loaded in the shell
# that just runs `cmake -S . -B build`.
set(ACES_MODULES "intel/2025.1.0 openmpi/5.0.7 python/3.12" CACHE STRING
    "module load set shared by aces2/aces3/pyaces")
set(ACES_NP 4 CACHE STRING
    "default MPI rank / parallel job count for the fast test suites")

option(ACES_BUILD_ACES2  "Build ACES II"                            ON)
option(ACES_BUILD_ACES3  "Build ACES III"                           ON)
option(ACES_BUILD_PYACES "Build pyaces (requires ACES_BUILD_ACES2)" ON)

if(ACES_BUILD_PYACES AND NOT ACES_BUILD_ACES2)
  message(FATAL_ERROR
    "ACES_BUILD_PYACES requires ACES_BUILD_ACES2 -- pyaces links against "
    "aces2's own built libraries (pyaces/build.sh's $ACES=.../aces2).")
endif()

set(ACES2_SOURCE_DIR  ${CMAKE_SOURCE_DIR}/aces2)
set(ACES3_SOURCE_DIR  ${CMAKE_SOURCE_DIR}/aces3)
set(PYACES_SOURCE_DIR ${CMAKE_SOURCE_DIR}/pyaces)

# aces2 installs in-place into its own source tree's bin/,lib/ (the classic
# ACES II layout pyaces/build.sh already expects at ../aces2/{bin,lib}) --
# just always THIS checkout's copy, never the hardcoded main-tree default
# that used to be baked into aces2/Makefiles/GNUmakefile.src. See the
# 2026-08-12 comment there for the landmine this fixes.
set(ACES2_INSTALL_BIN ${ACES2_SOURCE_DIR}/bin)
set(ACES2_INSTALL_LIB ${ACES2_SOURCE_DIR}/lib)

#----------------------------------------------------------------------
# aces2 (ACES II) -- GNUmake
#----------------------------------------------------------------------
if(ACES_BUILD_ACES2)
  ExternalProject_Add(aces2
    SOURCE_DIR        ${ACES2_SOURCE_DIR}
    DOWNLOAD_COMMAND  ""
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE   1
    BUILD_COMMAND     bash ${CMAKE_SOURCE_DIR}/cmake/build_aces2.sh
                       "${ACES_MODULES}" "${ACES2_SOURCE_DIR}"
                       "${ACES2_INSTALL_BIN}" "${ACES2_INSTALL_LIB}"
    INSTALL_COMMAND   ""
    BUILD_ALWAYS      1
    USES_TERMINAL_BUILD 1
  )
endif()

#----------------------------------------------------------------------
# aces3 (ACES III) -- autotools
#----------------------------------------------------------------------
if(ACES_BUILD_ACES3)
  ExternalProject_Add(aces3
    SOURCE_DIR        ${ACES3_SOURCE_DIR}
    DOWNLOAD_COMMAND  ""
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE   1
    BUILD_COMMAND     bash ${CMAKE_SOURCE_DIR}/cmake/build_aces3.sh
                       "${ACES_MODULES}" "${ACES3_SOURCE_DIR}"
    INSTALL_COMMAND   ""
    BUILD_ALWAYS      1
    USES_TERMINAL_BUILD 1
  )
endif()

#----------------------------------------------------------------------
# pyaces -- f2py + build.sh's own symbol-collision surgery
#----------------------------------------------------------------------
if(ACES_BUILD_PYACES)
  ExternalProject_Add(pyaces
    SOURCE_DIR        ${PYACES_SOURCE_DIR}
    DOWNLOAD_COMMAND  ""
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE   1
    BUILD_COMMAND     bash ${CMAKE_SOURCE_DIR}/cmake/build_pyaces.sh
                       "${PYACES_SOURCE_DIR}"
    INSTALL_COMMAND   ""
    BUILD_ALWAYS      1
    USES_TERMINAL_BUILD 1
    DEPENDS           aces2
  )
endif()

#----------------------------------------------------------------------
# Fast test suites -- ctest wraps each subsystem's own runner, then gates
# the combined rollup against the zero-point-20260812 baseline. Assumes
# the relevant subsystem(s) are already built (mirrors phase-1's
# run_tests.sh, which never invoked build.sh itself either).
#----------------------------------------------------------------------
enable_testing()

set(ACES_TEST_LOG_DIR ${CMAKE_BINARY_DIR}/test_logs)
file(MAKE_DIRECTORY ${ACES_TEST_LOG_DIR})

set(ACES2_TEST_LOG  ${ACES_TEST_LOG_DIR}/aces2_fast_suite.log)
set(ACES3_TEST_LOG  ${ACES_TEST_LOG_DIR}/aces3_fast_suite.log)
set(PYACES_TEST_LOG ${ACES_TEST_LOG_DIR}/pyaces_suite.log)

set(ACES_ROLLUP_FIXTURES "")

if(ACES_BUILD_ACES2)
  add_test(NAME aces2_fast_suite
    COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/run_suite_aces2.sh
            "${ACES_MODULES}" "${ACES2_SOURCE_DIR}" "${ACES_NP}" "${ACES2_TEST_LOG}")
  set_tests_properties(aces2_fast_suite PROPERTIES FIXTURES_SETUP aces2_ran)
  list(APPEND ACES_ROLLUP_FIXTURES aces2_ran)
endif()

if(ACES_BUILD_ACES3)
  add_test(NAME aces3_fast_suite
    COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/run_suite_aces3.sh
            "${ACES_MODULES}" "${ACES3_SOURCE_DIR}" "${ACES_NP}" "${ACES3_TEST_LOG}")
  set_tests_properties(aces3_fast_suite PROPERTIES FIXTURES_SETUP aces3_ran)
  list(APPEND ACES_ROLLUP_FIXTURES aces3_ran)
endif()

if(ACES_BUILD_PYACES)
  add_test(NAME pyaces_suite
    COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/run_suite_pyaces.sh
            "${ACES_MODULES}" "${PYACES_SOURCE_DIR}" "${PYACES_TEST_LOG}")
  set_tests_properties(pyaces_suite PROPERTIES FIXTURES_SETUP pyaces_ran)
  list(APPEND ACES_ROLLUP_FIXTURES pyaces_ran)
endif()

if(ACES_ROLLUP_FIXTURES)
  add_test(NAME rollup
    COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/rollup.sh
            "${ACES2_TEST_LOG}" "${ACES3_TEST_LOG}" "${PYACES_TEST_LOG}")
  set_tests_properties(rollup PROPERTIES FIXTURES_REQUIRED "${ACES_ROLLUP_FIXTURES}")
endif()
