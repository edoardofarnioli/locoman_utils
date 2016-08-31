cmake_minimum_required(VERSION 2.8.12)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 OLD)
    cmake_policy(SET CMP0005 NEW)
    cmake_policy(SET CMP0017 NEW)
endif(COMMAND cmake_policy)
include(ExternalProject)
PROJECT(locoman_utils)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

FIND_PACKAGE(YARP REQUIRED)
FIND_PACKAGE(idynutils REQUIRED)
# to avoid header not found when using idynutils 
FIND_PACKAGE(Eigen3 REQUIRED)

INCLUDE_DIRECTORIES(include ${YARP_INCLUDE_DIRS} ${idynutils_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

# for every file in INCLUDES CMake already sets the property HEADER_FILE_ONLY
file(GLOB_RECURSE ${PROJECT_NAME}_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/include/locoman/utils" *.h*)

ADD_LIBRARY(${PROJECT_NAME} SHARED  src/declarations.cpp
                        	    ${${PROJECT_NAME}_INCLUDES})

TARGET_LINK_LIBRARIES(${PROJECT_NAME}   ${YARP_LIBRARIES})

##############################################################################
# use YCM to export LOCOMAN_UTILS so that it can be found using find_package #
##############################################################################

set(VARS_PREFIX "locoman_utils")

set(${PROJECT_NAME}_MAJOR_VERSION 0)
set(${PROJECT_NAME}_MINOR_VERSION 0)
set(${PROJECT_NAME}_PATCH_VERSION 1)
set(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${${PROJECT_NAME}_PATCH_VERSION})

find_package(YCM REQUIRED)
include(YCMDefaultDirs)
ycm_default_dirs(${VARS_PREFIX})

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
                                      		  "$<INSTALL_INTERFACE:${${VARS_PREFIX}_INSTALL_INCLUDEDIR}>")

set_target_properties(${PROJECT_NAME}   PROPERTIES VERSION ${${VARS_PREFIX}_VERSION}
                            SOVERSION ${${VARS_PREFIX}_VERSION})

install(DIRECTORY include/ 
        DESTINATION "${${VARS_PREFIX}_INSTALL_INCLUDEDIR}"
        FILES_MATCHING PATTERN "*.h*")

install(TARGETS ${PROJECT_NAME}  
        EXPORT ${PROJECT_NAME}
        ARCHIVE DESTINATION "${${VARS_PREFIX}_INSTALL_BINDIR}" COMPONENT lib
        RUNTIME DESTINATION "${${VARS_PREFIX}_INSTALL_BINDIR}" COMPONENT bin
        LIBRARY DESTINATION "${${VARS_PREFIX}_INSTALL_LIBDIR}" COMPONENT shlib)
        
#enabling it will add all GYM dependencies as dependencies for third party users
set_property(GLOBAL APPEND PROPERTY ${VARS_PREFIX}_TARGETS ${PROJECT_NAME})


include(InstallBasicPackageFiles)
install_basic_package_files(${PROJECT_NAME}     VARS_PREFIX ${VARS_PREFIX}
                                    		VERSION ${${VARS_PREFIX}_VERSION}
                                    		COMPATIBILITY SameMajorVersion
                                    		TARGETS_PROPERTY ${VARS_PREFIX}_TARGETS
                                    		NO_CHECK_REQUIRED_COMPONENTS_MACRO)

include(AddUninstallTarget)

