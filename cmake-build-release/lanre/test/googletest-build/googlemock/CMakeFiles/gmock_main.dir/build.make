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
include lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/depend.make

# Include the progress variables for this target.
include lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/progress.make

# Include the compile flags for this target's objects.
include lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/flags.make

lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o: lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/flags.make
lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o: lanre/test/googletest-src/googlemock/src/gmock_main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o"
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-build/googlemock && /usr/local/opt/llvm/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmock_main.dir/src/gmock_main.cc.o -c /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-src/googlemock/src/gmock_main.cc

lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock_main.dir/src/gmock_main.cc.i"
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-build/googlemock && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-src/googlemock/src/gmock_main.cc > CMakeFiles/gmock_main.dir/src/gmock_main.cc.i

lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock_main.dir/src/gmock_main.cc.s"
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-build/googlemock && /usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-src/googlemock/src/gmock_main.cc -o CMakeFiles/gmock_main.dir/src/gmock_main.cc.s

# Object files for target gmock_main
gmock_main_OBJECTS = \
"CMakeFiles/gmock_main.dir/src/gmock_main.cc.o"

# External object files for target gmock_main
gmock_main_EXTERNAL_OBJECTS =

lib/libgmock_main.a: lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/src/gmock_main.cc.o
lib/libgmock_main.a: lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/build.make
lib/libgmock_main.a: lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ../../../../lib/libgmock_main.a"
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-build/googlemock && $(CMAKE_COMMAND) -P CMakeFiles/gmock_main.dir/cmake_clean_target.cmake
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-build/googlemock && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gmock_main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/build: lib/libgmock_main.a

.PHONY : lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/build

lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/clean:
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-build/googlemock && $(CMAKE_COMMAND) -P CMakeFiles/gmock_main.dir/cmake_clean.cmake
.PHONY : lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/clean

lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/depend:
	cd /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/loganmorrison/Documents/git_hub/Haliax /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-src/googlemock /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-build/googlemock /Users/loganmorrison/Documents/git_hub/Haliax/cmake-build-release/lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lanre/test/googletest-build/googlemock/CMakeFiles/gmock_main.dir/depend

