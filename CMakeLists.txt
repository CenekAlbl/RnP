cmake_minimum_required (VERSION 2.6)
project (RnP)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()


set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if (EIGEN_FOUND)
  include_directories(${EIGEN_INCLUDE_DIRS})
endif (EIGEN_FOUND)

if(NOT DEFINED Python_EXECUTABLE)
    find_package (Python COMPONENTS Interpreter Development)
    if(Python_FOUND)
        message(WARNING "Python ${Python_VERSION} found and will be used. If you want to use a different version, please specify Python_EXECUTABLE variable")
    else()
        message(WARNING "Python was not found automatically. If you want to compile the python bindings, please specify the Python_EXECUTABLE variable")
    endif(Python_FOUND)
else()
    message(WARNING "Python executable specified by user: ${Python_EXECUTABLE}")
    set(Python_FOUND True)
endif()

find_package(Matlab COMPONENTS MX_LIBRARY MEX_COMPILER)




include_directories(
    c++
)


add_library(rnp STATIC
  c++/r6pDoubleLin.cpp c++/r6pSingleLin.cpp c++/R7Pf.cpp c++/R7Pfr.cpp c++/iterative.cpp c++/rnp.cpp c++/sturm.cpp c++/utils.cpp 
)

target_include_directories(rnp PRIVATE .)
target_include_directories(rnp PUBLIC ${EIGEN_INCLUDE_DIRS})

set_target_properties(rnp PROPERTIES PUBLIC_HEADER c++/rnp.h)
set_target_properties(rnp PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_link_libraries (rnp Eigen3::Eigen)

install(TARGETS rnp
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION arch
    PUBLIC_HEADER DESTINATION include)

enable_testing()

ADD_EXECUTABLE(test_r6p_doublelin tests/test_r6p_doublelin.cpp)
ADD_EXECUTABLE(test_r6p_singlelin tests/test_r6p_doublelin.cpp)
ADD_EXECUTABLE(test_r7pf tests/test_r7pf.cpp)
ADD_EXECUTABLE(test_r7pfr tests/test_r7pfr.cpp)
ADD_EXECUTABLE(test_rs_projections tests/test_rs_projections.cpp)

ADD_TEST(test_r6p_doublelin test_r6p_doublelin)
ADD_TEST(test_r6p_singlelin test_r6p_doublelin)
ADD_TEST(test_r7pf test_r7pf)
ADD_TEST(test_r7pfr test_r7pfr)
ADD_TEST(test_rs_projections test_rs_projections)

target_include_directories(test_r6p_doublelin PRIVATE .)
target_include_directories(test_r6p_doublelin PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (test_r6p_doublelin Eigen3::Eigen rnp)

target_include_directories(test_r6p_singlelin PRIVATE .)
target_include_directories(test_r6p_singlelin PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (test_r6p_singlelin Eigen3::Eigen rnp)

target_include_directories(test_r7pf PRIVATE .)
target_include_directories(test_r7pf PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (test_r7pf Eigen3::Eigen rnp)

target_include_directories(test_r7pfr PRIVATE .)
target_include_directories(test_r7pfr PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (test_r7pfr Eigen3::Eigen rnp)

target_include_directories(test_rs_projections PRIVATE .)
target_include_directories(test_rs_projections PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (test_rs_projections Eigen3::Eigen rnp)


if(Matlab_FOUND)
    matlab_add_mex(
        NAME r6pSingleLin
        SRC "matlab/r6pSingleLinMex.cpp"
        LINK_TO rnp
    )
endif(Matlab_FOUND)

if(Python_FOUND)

    set(SETUP_PY_IN "${CMAKE_SOURCE_DIR}/cython/setup.py.in")
    set(SETUP_PY "${CMAKE_SOURCE_DIR}/cython/setup.py")

    #set(PY_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build/pytimestamp")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    add_custom_command(OUTPUT "rnppy" 
        COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cython/setup.py build_ext --inplace --include-dirs "${CMAKE_SOURCE_DIR}/c++/" --library-dirs "${CMAKE_CURRENT_BINARY_DIR}"
        DEPENDS rnp)

    target_link_libraries (rnp)

    add_custom_target(rnppyt ALL DEPENDS rnppy rnp)

    install(CODE "execute_process(COMMAND ${Python_EXECUTABLE} ${SETUP_PY} install)")

endif(Python_FOUND)
