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
CMAKE_SOURCE_DIR = /Users/loganmorrison/Desktop/Haliax

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/loganmorrison/Desktop/Haliax/cmake-build-release

# Include any dependencies generated for this target.
include lanre/test/CMakeFiles/test_gam.dir/depend.make

# Include the progress variables for this target.
include lanre/test/CMakeFiles/test_gam.dir/progress.make

# Include the compile flags for this target's objects.
include lanre/test/CMakeFiles/test_gam.dir/flags.make

lanre/test/CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.o: lanre/test/CMakeFiles/test_gam.dir/flags.make
lanre/test/CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.o: ../lanre/test/test_interpolate_univariate_spline.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/loganmorrison/Desktop/Haliax/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lanre/test/CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.o"
	cd /Users/loganmorrison/Desktop/Haliax/cmake-build-release/lanre/test && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.o -c /Users/loganmorrison/Desktop/Haliax/lanre/test/test_interpolate_univariate_spline.cpp

lanre/test/CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.i"
	cd /Users/loganmorrison/Desktop/Haliax/cmake-build-release/lanre/test && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/loganmorrison/Desktop/Haliax/lanre/test/test_interpolate_univariate_spline.cpp > CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.i

lanre/test/CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.s"
	cd /Users/loganmorrison/Desktop/Haliax/cmake-build-release/lanre/test && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/loganmorrison/Desktop/Haliax/lanre/test/test_interpolate_univariate_spline.cpp -o CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.s

# Object files for target test_gam
test_gam_OBJECTS = \
"CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.o"

# External object files for target test_gam
test_gam_EXTERNAL_OBJECTS =

bin/test_gam: lanre/test/CMakeFiles/test_gam.dir/test_interpolate_univariate_spline.cpp.o
bin/test_gam: lanre/test/CMakeFiles/test_gam.dir/build.make
bin/test_gam: lanre/test/CMakeFiles/test_gam.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/loganmorrison/Desktop/Haliax/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/test_gam"
	cd /Users/loganmorrison/Desktop/Haliax/cmake-build-release/lanre/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_gam.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lanre/test/CMakeFiles/test_gam.dir/build: bin/test_gam

.PHONY : lanre/test/CMakeFiles/test_gam.dir/build

lanre/test/CMakeFiles/test_gam.dir/clean:
	cd /Users/loganmorrison/Desktop/Haliax/cmake-build-release/lanre/test && $(CMAKE_COMMAND) -P CMakeFiles/test_gam.dir/cmake_clean.cmake
.PHONY : lanre/test/CMakeFiles/test_gam.dir/clean

lanre/test/CMakeFiles/test_gam.dir/depend:
	cd /Users/loganmorrison/Desktop/Haliax/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/loganmorrison/Desktop/Haliax /Users/loganmorrison/Desktop/Haliax/lanre/test /Users/loganmorrison/Desktop/Haliax/cmake-build-release /Users/loganmorrison/Desktop/Haliax/cmake-build-release/lanre/test /Users/loganmorrison/Desktop/Haliax/cmake-build-release/lanre/test/CMakeFiles/test_gam.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lanre/test/CMakeFiles/test_gam.dir/depend
