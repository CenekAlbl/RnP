cmake_minimum_required (VERSION 2.6)
project (RnP)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if (EIGEN_FOUND)
  include_directories(${EIGEN_INCLUDE_DIRS})
endif (EIGEN_FOUND)

find_package (Python COMPONENTS Interpreter Development)



include_directories(
    c++
)


add_library(rnp SHARED
  c++/r6p.cpp c++/sturm.cpp
)

target_include_directories(rnp PRIVATE .)
target_include_directories(rnp PUBLIC ${EIGEN_INCLUDE_DIRS})

set_target_properties(rnp PROPERTIES PUBLIC_HEADER c++/rnp.h)

target_link_libraries (rnp Eigen3::Eigen)

install(TARGETS rnp
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include)


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