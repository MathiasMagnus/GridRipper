###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Example2-HiggsParticle example
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 2014/08/01
#
###############################################################################

# The project name and version, release numbers
project (STL-Example2-HiggsParticle)


# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
include_directories (${PROJECT_SOURCE_DIR}/inc/)
include_directories (${Gripper_ADDITIONAL_INCLUDE_DIRS})


# Adding source code files according to configuration
set (STL_Example2_HiggsParticle_BUILD)
set (STL_Example2_HiggsParticle_HDRS)
set (STL_Example2_HiggsParticle_SRCS)

set (STL_Example2_HiggsParticle_HDRS ${STL_Example2_HiggsParticle_HDRS} ${PROJECT_SOURCE_DIR}/inc/STL_Example2_HiggsParticle.hpp)

set (STL_Example2_HiggsParticle_SRCS ${STL_Example2_HiggsParticle_SRCS} ${PROJECT_SOURCE_DIR}/src/STL_Example2_HiggsParticle.cpp)

set (STL_Example2_HiggsParticle_BUILD ${STL_Example2_HiggsParticle_BUILD} ${STL_Example2_HiggsParticle_HDRS} ${STL_Example2_HiggsParticle_SRCS})


# Adding library target for build
add_executable (STL-Example2-HiggsParticle ${STL_Example2_HiggsParticle_BUILD})
target_link_libraries (STL-Example2-HiggsParticle Gripper ${Gripper_ADDITIONAL_LIBS})


# Create project groups for IDEs
set_target_properties (STL-Example2-HiggsParticle PROPERTIES FOLDER "Examples//STL")

# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Example2_HiggsParticle_HDRS})
source_group ("Sources" FILES ${STL_Example2_HiggsParticle_SRCS})


# Add platform specific compile options
if(MSVC_VERSION)
#set_target_properties (STL-Example2-HiggsParticle PROPERTIES
#                       COMPILE_FLAGS ${MSVC_COMPILER_FLAGS}
#                       LINK_FLAGS ${MSVC_LINKER_FLAGS})
endif(MSVC_VERSION)