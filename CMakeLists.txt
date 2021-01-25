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


find_package(Matlab COMPONENTS MX_LIBRARY MEX_COMPILER)




include_directories(
    c++
)


add_library(rnp STATIC
  c++/r6p_2lin.cpp c++/r6p_1lin.cpp c++/r7pf.cpp c++/r7pfr.cpp c++/r6p_iter.cpp c++/iterative.cpp c++/rnp.cpp c++/sturm.cpp c++/utils.cpp
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

ADD_EXECUTABLE(test_r6p_2lin tests/test_r6p_2lin.cpp)
ADD_EXECUTABLE(test_r6p_1lin tests/test_r6p_1lin.cpp)
ADD_EXECUTABLE(test_r6p_iter tests/test_r6p_iter.cpp)
ADD_EXECUTABLE(test_r7pf tests/test_r7pf.cpp)
ADD_EXECUTABLE(test_r7pfr tests/test_r7pfr.cpp)
ADD_EXECUTABLE(test_rs_projections tests/test_rs_projections.cpp)

ADD_EXECUTABLE(benchmark c++/benchmark.cpp)

ADD_TEST(test_r6p_2lin test_r6p_2lin)
ADD_TEST(test_r6p_1lin test_r6p_1lin)
ADD_TEST(test_r6p_iter test_r6p_iter)
ADD_TEST(test_r7pf test_r7pf)
ADD_TEST(test_r7pfr test_r7pfr)
ADD_TEST(test_rs_projections test_rs_projections)

target_include_directories(test_r6p_2lin PRIVATE .)
target_include_directories(test_r6p_2lin PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (test_r6p_2lin Eigen3::Eigen rnp)

target_include_directories(test_r6p_1lin PRIVATE .)
target_include_directories(test_r6p_1lin PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (test_r6p_1lin Eigen3::Eigen rnp)

target_include_directories(test_r6p_iter PRIVATE .)
target_include_directories(test_r6p_iter PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (test_r6p_iter Eigen3::Eigen rnp)

target_include_directories(test_r7pf PRIVATE .)
target_include_directories(test_r7pf PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (test_r7pf Eigen3::Eigen rnp)

target_include_directories(test_r7pfr PRIVATE .)
target_include_directories(test_r7pfr PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (test_r7pfr Eigen3::Eigen rnp)

target_include_directories(test_rs_projections PRIVATE .)
target_include_directories(test_rs_projections PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (test_rs_projections Eigen3::Eigen rnp)

target_include_directories(benchmark PRIVATE .)
target_include_directories(benchmark PUBLIC ${EIGEN_INCLUDE_DIRS})
target_link_libraries (benchmark Eigen3::Eigen rnp)


if(Matlab_FOUND)
    matlab_add_mex(
        NAME r6p1lin
        SRC "matlab/r6p1lin.cpp"
        LINK_TO rnp
    )
    matlab_add_mex(
        NAME r6p2lin
        SRC "matlab/r6p2lin.cpp"
        LINK_TO rnp
    )
    matlab_add_mex(
	    NAME r6pIter
	SRC "matlab/r6pIter.cpp"
        LINK_TO rnp
    )
    matlab_add_mex(
        NAME r7pf
        SRC "matlab/r7pf.cpp"
        LINK_TO rnp
    )
    matlab_add_mex(
        NAME r7pfr
        SRC "matlab/r7pfr.cpp"
        LINK_TO rnp
    )
endif(Matlab_FOUND)

if(PYTHON_BINDINGS)

	if(NOT DEFINED Python_EXECUTABLE)
	    find_package (Python COMPONENTS Interpreter Development)
	    if(Python_FOUND)
		message(WARNING "Python ${Python_VERSION} found and will be used. If you want to use a different version, please specify Python_EXECUTABLE variable")
	    else()
		message(WARNING "Python was not found automatically. If you want to compile the python bindings, please specify the Python_EXECUTABLE variable")
	    endif(Python_FOUND)
	else()
	    message(WARNING "Python executable specified by user: ${Python_EXECUTABLE}")
	endif()

    	add_subdirectory(extern/pybind11)

	pybind11_add_module(pyrnp EXCLUDE_FROM_ALL c++/python_bindings.cpp)

	target_include_directories(pyrnp PRIVATE .)

	target_include_directories(pyrnp PUBLIC ${EIGEN_INCLUDE_DIRS})

	target_link_libraries(pyrnp PRIVATE Eigen3::Eigen rnp)

endif(PYTHON_BINDINGS)
