# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = "/Users/loganmorrison/Library/Application Support/JetBrains/Toolbox/apps/CLion/ch-0/193.6911.21/CLion.app/Contents/bin/cmake/mac/bin/cmake"

# The command to remove a file.
RM = "/Users/loganmorrison/Library/Application Support/JetBrains/Toolbox/apps/CLion/ch-0/193.6911.21/CLion.app/Contents/bin/cmake/mac/bin/cmake" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/loganmorrison/Documents/git_hub/Haliax

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release

# Include any dependencies generated for this target.
include lanre/test/CMakeFiles/test_dm_models_darksun.dir/depend.make

# Include the progress variables for this target.
include lanre/test/CMakeFiles/test_dm_models_darksun.dir/progress.make

# Include the compile flags for this target's objects.
include lanre/test/CMakeFiles/test_dm_models_darksun.dir/flags.make

lanre/test/CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.o: lanre/test/CMakeFiles/test_dm_models_darksun.dir/flags.make
lanre/test/CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.o: ../lanre/test/test_dm_models_darksun.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lanre/test/CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.o"
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.o -c /Users/loganmorrison/Documents/git_hub/Haliax/lanre/test/test_dm_models_darksun.cpp

lanre/test/CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.i"
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/loganmorrison/Documents/git_hub/Haliax/lanre/test/test_dm_models_darksun.cpp > CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.i

lanre/test/CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.s"
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/loganmorrison/Documents/git_hub/Haliax/lanre/test/test_dm_models_darksun.cpp -o CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.s

# Object files for target test_dm_models_darksun
test_dm_models_darksun_OBJECTS = \
"CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.o"

# External object files for target test_dm_models_darksun
test_dm_models_darksun_EXTERNAL_OBJECTS =

../bin/test_dm_models_darksun: lanre/test/CMakeFiles/test_dm_models_darksun.dir/test_dm_models_darksun.cpp.o
../bin/test_dm_models_darksun: lanre/test/CMakeFiles/test_dm_models_darksun.dir/build.make
../bin/test_dm_models_darksun: lib/libgtest_main.a
../bin/test_dm_models_darksun: lib/libgtest.a
../bin/test_dm_models_darksun: lanre/test/CMakeFiles/test_dm_models_darksun.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../bin/test_dm_models_darksun"
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_dm_models_darksun.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lanre/test/CMakeFiles/test_dm_models_darksun.dir/build: ../bin/test_dm_models_darksun

.PHONY : lanre/test/CMakeFiles/test_dm_models_darksun.dir/build

lanre/test/CMakeFiles/test_dm_models_darksun.dir/clean:
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test && $(CMAKE_COMMAND) -P CMakeFiles/test_dm_models_darksun.dir/cmake_clean.cmake
.PHONY : lanre/test/CMakeFiles/test_dm_models_darksun.dir/clean

lanre/test/CMakeFiles/test_dm_models_darksun.dir/depend:
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/loganmorrison/Documents/git_hub/Haliax /Users/loganmorrison/Documents/git_hub/Haliax/lanre/test /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/CMakeFiles/test_dm_models_darksun.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lanre/test/CMakeFiles/test_dm_models_darksun.dir/depend

