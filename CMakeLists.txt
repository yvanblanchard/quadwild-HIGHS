cmake_minimum_required(VERSION 3.18)
project(quadwild LANGUAGES C CXX)

# Command line QuadWild application (mirrors quadwild/quadwild.pro + libs/libs.pri + libs/quadretopology/*.pri)
# The quadrangulation ILP is solved with HiGHS (MIT) instead of Gurobi.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(QUADWILD_USE_COMISO "Use CoMISo for the field computation (requires a prebuilt CoMISo + BLAS)" OFF)
set(QUADWILD_HIGHS_VERSION "v1.11.0" CACHE STRING "HiGHS git tag downloaded when HiGHS is not installed")
set(BOOST_INCLUDE_DIR "" CACHE PATH "Boost headers directory (downloaded automatically if empty and not found)")

set(QW_ROOT             ${CMAKE_CURRENT_SOURCE_DIR})
set(LIBS_PATH           ${QW_ROOT}/libs)
set(QUADRETOPOLOGY_PATH ${LIBS_PATH}/quadretopology)
set(PATTERNS_PATH       ${QUADRETOPOLOGY_PATH}/patterns)
set(LPSOLVER_PATH       ${PATTERNS_PATH}/patterns/lp_solve_5.5)
set(OPENMESH_PATH       ${PATTERNS_PATH}/patterns/OpenMesh)
set(PATCHGEN_PATH       ${PATTERNS_PATH}/patterns/patchgen)
set(LIBIGL_PATH         ${LIBS_PATH}/libigl)
set(VCGLIB_PATH         ${LIBS_PATH}/vcglib)
set(EIGEN_PATH          ${LIBS_PATH}/eigen)
set(COMISO_PATH         ${LIBS_PATH}/CoMISo)
set(XFIELDTRACER_PATH   ${LIBS_PATH}/xfield_tracer)
set(MESHFIELD_PATH      ${QW_ROOT}/components/field_computation)
set(MESHTRACE_PATH      ${QW_ROOT}/components/field_tracing)
set(QUADRANGULATE_PATH  ${QW_ROOT}/components/quad_from_patches)

if(NOT EXISTS ${VCGLIB_PATH}/vcg OR NOT EXISTS ${XFIELDTRACER_PATH}/tracing)
    message(FATAL_ERROR "Submodules missing: run 'git submodule update --init --recursive'")
endif()

# ---------------------------------------------------------------- Boost (header only)
if(NOT BOOST_INCLUDE_DIR)
    if(POLICY CMP0167)
        cmake_policy(SET CMP0167 NEW)
    endif()
    find_package(Boost CONFIG QUIET)
    if(Boost_FOUND)
        set(BOOST_INCLUDE_DIR ${Boost_INCLUDE_DIRS})
    else()
        # Extract only the boost/ headers: the full archive has doc paths longer than
        # the Windows MAX_PATH limit, which breaks extraction/cleanup of the build directory
        set(_boost_root ${CMAKE_BINARY_DIR}/_deps/boost)
        set(BOOST_INCLUDE_DIR ${_boost_root}/boost_1_86_0)
        if(NOT EXISTS ${BOOST_INCLUDE_DIR}/boost/version.hpp)
            message(STATUS "Boost not found: downloading headers")
            file(MAKE_DIRECTORY ${_boost_root})
            file(DOWNLOAD https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
                 ${_boost_root}/boost.tar.gz STATUS _boost_dl)
            list(GET _boost_dl 0 _boost_dl_code)
            if(NOT _boost_dl_code EQUAL 0)
                message(FATAL_ERROR "Boost download failed: ${_boost_dl}. Set -DBOOST_INCLUDE_DIR=...")
            endif()
            execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf boost.tar.gz boost_1_86_0/boost
                            WORKING_DIRECTORY ${_boost_root} RESULT_VARIABLE _boost_tar)
            file(REMOVE ${_boost_root}/boost.tar.gz)
            if(NOT EXISTS ${BOOST_INCLUDE_DIR}/boost/version.hpp)
                message(FATAL_ERROR "Boost extraction failed (${_boost_tar})")
            endif()
        endif()
    endif()
endif()

# ---------------------------------------------------------------- HiGHS (MIT) - MILP solver
# Uses an installed HiGHS (find_package) or downloads and builds it as a static library.
find_package(highs CONFIG QUIET)
if(NOT highs_FOUND)
    include(FetchContent)
    message(STATUS "HiGHS not found: downloading ${QUADWILD_HIGHS_VERSION}")
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(BUILD_CXX_EXE OFF CACHE BOOL "" FORCE)
    set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(ZLIB OFF CACHE BOOL "" FORCE)
    set(FAST_BUILD ON CACHE BOOL "" FORCE)
    FetchContent_Declare(highs
        GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
        GIT_TAG ${QUADWILD_HIGHS_VERSION}
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(highs)
endif()

# ---------------------------------------------------------------- lp_solve (C)
set(LPSOLVE_SOURCES
    bfp/bfp_LUSOL/lp_LUSOL.c bfp/bfp_LUSOL/LUSOL/lusol.c colamd/colamd.c fortify.c ini.c
    lp_Hash.c lp_MDO.c lp_MPS.c lp_SOS.c lp_crash.c lp_lib.c lp_matrix.c lp_mipbb.c
    lp_params.c lp_presolve.c lp_price.c lp_pricePSE.c lp_report.c lp_rlp.c lp_scale.c
    lp_simplex.c lp_utils.c lp_wlp.c shared/commonlib.c shared/mmio.c shared/myblas.c yacc_read.c)
list(TRANSFORM LPSOLVE_SOURCES PREPEND ${LPSOLVER_PATH}/)

set(LPSOLVE_INCLUDES
    ${LPSOLVER_PATH} ${LPSOLVER_PATH}/bfp ${LPSOLVER_PATH}/bfp/bfp_LUSOL
    ${LPSOLVER_PATH}/bfp/bfp_LUSOL/LUSOL ${LPSOLVER_PATH}/colamd ${LPSOLVER_PATH}/shared)
set(LPSOLVE_DEFINES
    LoadInverseLib=FALSE LoadLanguageLib=FALSE LoadableBlasLib=0 INVERSE_ACTIVE=INVERSE_LUSOL
    RoleIsExternalInvEngine NOFORTIFY CHECK_SOLUTION YY_NEVER_INTERACTIVE PARSER_LP LP_MAXLINELEN=0)

add_library(qw_lpsolve STATIC ${LPSOLVE_SOURCES})
target_include_directories(qw_lpsolve PUBLIC ${LPSOLVE_INCLUDES})
target_compile_definitions(qw_lpsolve PUBLIC ${LPSOLVE_DEFINES})
if(MSVC)
    target_compile_definitions(qw_lpsolve PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(qw_lpsolve PRIVATE /w)
else()
    target_compile_options(qw_lpsolve PRIVATE -w)
endif()

# ---------------------------------------------------------------- quadwild executable
set(OPENMESH_SOURCES
    Core/Geometry/NormalConeT.cc Core/IO/BinaryHelper.cc Core/IO/IOManager.cc Core/IO/OMFormat.cc
    Core/IO/OMFormatT.cc Core/IO/reader/BaseReader.cc Core/IO/reader/OBJReader.cc
    Core/IO/reader/OFFReader.cc Core/IO/reader/OMReader.cc Core/IO/reader/PLYReader.cc
    Core/IO/reader/STLReader.cc Core/IO/writer/BaseWriter.cc Core/IO/writer/OBJWriter.cc
    Core/IO/writer/OFFWriter.cc Core/IO/writer/OMWriter.cc Core/IO/writer/PLYWriter.cc
    Core/IO/writer/STLWriter.cc Core/IO/writer/VTKWriter.cc Core/Mesh/ArrayKernel.cc
    Core/Mesh/ArrayKernelT.cc Core/Mesh/BaseKernel.cc Core/Mesh/PolyConnectivity.cc
    Core/Mesh/PolyMeshT.cc Core/Mesh/TriConnectivity.cc Core/Mesh/TriMeshT.cc
    Core/System/omstream.cc Core/Utils/BaseProperty.cc Core/Utils/Endian.cc
    Core/Utils/RandomNumberGenerator.cc Core/Utils/SingletonT.cc)
list(TRANSFORM OPENMESH_SOURCES PREPEND ${OPENMESH_PATH}/)
# qmake compiles .cc as C++; make sure CMake does too
set_source_files_properties(${OPENMESH_SOURCES} PROPERTIES LANGUAGE CXX)

set(PATCHGEN_SOURCES
    adjust_parameter.cpp get_boundary_geometry.cpp get_constraint_matrix.cpp get_constraint_rhs.cpp
    get_default_parameter.cpp get_param_str.cpp get_variable.cpp get_variable_indicators.cpp
    switch_pattern.cpp)
list(TRANSFORM PATCHGEN_SOURCES PREPEND ${PATCHGEN_PATH}/)

# Only the sources that are not #included by their own header: principal_curvature, quad_from_patches,
# quadretopology, functions and qr_{charts,convert,ilp,patterns,utils} are header-included (templates)
set(QUADRETOPOLOGY_SOURCES
    ${QUADRETOPOLOGY_PATH}/quadretopology/includes/qr_mapping.cpp
    ${PATTERNS_PATH}/patterns/Patch.cpp
    ${PATTERNS_PATH}/patterns/determine_geometry.cpp
    ${PATTERNS_PATH}/patterns/generate_patch.cpp)

add_executable(quadwild
    ${QW_ROOT}/quadwild/quadwild.cpp
    ${QUADRANGULATE_PATH}/load_save.cpp
    ${VCGLIB_PATH}/wrap/ply/plylib.cpp
    ${QUADRETOPOLOGY_SOURCES}
    ${PATCHGEN_SOURCES}
    ${OPENMESH_SOURCES})

target_include_directories(quadwild PRIVATE
    ${QW_ROOT}/quadwild
    ${QUADRETOPOLOGY_PATH}
    ${PATTERNS_PATH}
    ${PATTERNS_PATH}/patterns
    ${VCGLIB_PATH}
    ${EIGEN_PATH}
    ${MESHFIELD_PATH}
    ${MESHTRACE_PATH}
    ${QUADRANGULATE_PATH}
    ${XFIELDTRACER_PATH})
target_include_directories(quadwild SYSTEM PRIVATE ${LIBIGL_PATH}/include ${BOOST_INCLUDE_DIR})

target_compile_definitions(quadwild PRIVATE GLEW_STATIC INCLUDE_TEMPLATES $<$<NOT:$<CONFIG:Debug>>:NDEBUG>)
target_link_libraries(quadwild PRIVATE qw_lpsolve highs::highs)

if(QUADWILD_USE_COMISO)
    find_library(COMISO_LIB NAMES CoMISo PATHS ${COMISO_PATH}/build/Build/lib/CoMISo REQUIRED)
    target_compile_definitions(quadwild PRIVATE COMISO_FIELD)
    target_include_directories(quadwild PRIVATE ${COMISO_PATH}/.. ${COMISO_PATH}/gmm/include)
    target_sources(quadwild PRIVATE ${LIBIGL_PATH}/include/igl/copyleft/comiso/nrosy.cpp)
    target_link_libraries(quadwild PRIVATE ${COMISO_LIB})
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(quadwild PRIVATE OpenMP::OpenMP_CXX)
endif()

if(MSVC)
    # _HAS_AUTO_PTR_ETC: the bundled Eigen uses std::binder2nd, removed from the C++17 MSVC STL
    target_compile_definitions(quadwild PRIVATE NOMINMAX _USE_MATH_DEFINES _SCL_SECURE_NO_DEPRECATE
                                                _CRT_SECURE_NO_WARNINGS _HAS_AUTO_PTR_ETC=1)
    target_compile_options(quadwild PRIVATE /bigobj /permissive- /Zc:__cplusplus /W0 /EHsc)
else()
    target_compile_options(quadwild PRIVATE -w)
endif()

# basic_setup.txt is read from the working directory
add_custom_command(TARGET quadwild POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${QW_ROOT}/quadwild/basic_setup.txt
            ${QW_ROOT}/quadwild/basic_setup_mechanical.txt
            ${QW_ROOT}/quadwild/basic_setup_organic.txt
            $<TARGET_FILE_DIR:quadwild>)

install(TARGETS quadwild RUNTIME DESTINATION bin)
install(FILES ${QW_ROOT}/quadwild/basic_setup.txt ${QW_ROOT}/quadwild/basic_setup_mechanical.txt
              ${QW_ROOT}/quadwild/basic_setup_organic.txt DESTINATION bin)
