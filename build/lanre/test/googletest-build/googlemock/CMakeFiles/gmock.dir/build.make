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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.15.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.15.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/loganmorrison/Desktop/Haliax

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/loganmorrison/Desktop/Haliax/build

# Include any dependencies generated for this target.
include lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/depend.make

# Include the progress variables for this target.
include lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/progress.make

# Include the compile flags for this target's objects.
include lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/flags.make

lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/src/gmock-all.cc.o: lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/flags.make
lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/src/gmock-all.cc.o: lanre/test/googletest-src/googlemock/src/gmock-all.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/loganmorrison/Desktop/Haliax/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/src/gmock-all.cc.o"
	cd /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-build/googlemock && /usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmock.dir/src/gmock-all.cc.o -c /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-src/googlemock/src/gmock-all.cc

lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/src/gmock-all.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock.dir/src/gmock-all.cc.i"
	cd /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-build/googlemock && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-src/googlemock/src/gmock-all.cc > CMakeFiles/gmock.dir/src/gmock-all.cc.i

lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/src/gmock-all.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock.dir/src/gmock-all.cc.s"
	cd /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-build/googlemock && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-src/googlemock/src/gmock-all.cc -o CMakeFiles/gmock.dir/src/gmock-all.cc.s

# Object files for target gmock
gmock_OBJECTS = \
"CMakeFiles/gmock.dir/src/gmock-all.cc.o"

# External object files for target gmock
gmock_EXTERNAL_OBJECTS =

lib/libgmockd.a: lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/src/gmock-all.cc.o
lib/libgmockd.a: lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/build.make
lib/libgmockd.a: lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/loganmorrison/Desktop/Haliax/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ../../../../lib/libgmockd.a"
	cd /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-build/googlemock && $(CMAKE_COMMAND) -P CMakeFiles/gmock.dir/cmake_clean_target.cmake
	cd /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-build/googlemock && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gmock.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/build: lib/libgmockd.a

.PHONY : lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/build

lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/clean:
	cd /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-build/googlemock && $(CMAKE_COMMAND) -P CMakeFiles/gmock.dir/cmake_clean.cmake
.PHONY : lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/clean

lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/depend:
	cd /Users/loganmorrison/Desktop/Haliax/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/loganmorrison/Desktop/Haliax /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-src/googlemock /Users/loganmorrison/Desktop/Haliax/build /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-build/googlemock /Users/loganmorrison/Desktop/Haliax/build/lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lanre/test/googletest-build/googlemock/CMakeFiles/gmock.dir/depend
