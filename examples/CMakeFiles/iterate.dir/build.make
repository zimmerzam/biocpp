# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/zimmer/Desktop/biocpp/examples

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zimmer/Desktop/biocpp/examples

# Include any dependencies generated for this target.
include CMakeFiles/iterate.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/iterate.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/iterate.dir/flags.make

CMakeFiles/iterate.dir/iterate.cpp.o: CMakeFiles/iterate.dir/flags.make
CMakeFiles/iterate.dir/iterate.cpp.o: iterate.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/zimmer/Desktop/biocpp/examples/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/iterate.dir/iterate.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) -DBIOCPP_INCLUDE_PDB $(CXX_FLAGS) -o CMakeFiles/iterate.dir/iterate.cpp.o -c /home/zimmer/Desktop/biocpp/examples/iterate.cpp

CMakeFiles/iterate.dir/iterate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/iterate.dir/iterate.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) -DBIOCPP_INCLUDE_PDB $(CXX_FLAGS) -E /home/zimmer/Desktop/biocpp/examples/iterate.cpp > CMakeFiles/iterate.dir/iterate.cpp.i

CMakeFiles/iterate.dir/iterate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/iterate.dir/iterate.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) -DBIOCPP_INCLUDE_PDB $(CXX_FLAGS) -S /home/zimmer/Desktop/biocpp/examples/iterate.cpp -o CMakeFiles/iterate.dir/iterate.cpp.s

CMakeFiles/iterate.dir/iterate.cpp.o.requires:
.PHONY : CMakeFiles/iterate.dir/iterate.cpp.o.requires

CMakeFiles/iterate.dir/iterate.cpp.o.provides: CMakeFiles/iterate.dir/iterate.cpp.o.requires
	$(MAKE) -f CMakeFiles/iterate.dir/build.make CMakeFiles/iterate.dir/iterate.cpp.o.provides.build
.PHONY : CMakeFiles/iterate.dir/iterate.cpp.o.provides

CMakeFiles/iterate.dir/iterate.cpp.o.provides.build: CMakeFiles/iterate.dir/iterate.cpp.o

# Object files for target iterate
iterate_OBJECTS = \
"CMakeFiles/iterate.dir/iterate.cpp.o"

# External object files for target iterate
iterate_EXTERNAL_OBJECTS =

iterate: CMakeFiles/iterate.dir/iterate.cpp.o
iterate: CMakeFiles/iterate.dir/build.make
iterate: CMakeFiles/iterate.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable iterate"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/iterate.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/iterate.dir/build: iterate
.PHONY : CMakeFiles/iterate.dir/build

CMakeFiles/iterate.dir/requires: CMakeFiles/iterate.dir/iterate.cpp.o.requires
.PHONY : CMakeFiles/iterate.dir/requires

CMakeFiles/iterate.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/iterate.dir/cmake_clean.cmake
.PHONY : CMakeFiles/iterate.dir/clean

CMakeFiles/iterate.dir/depend:
	cd /home/zimmer/Desktop/biocpp/examples && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zimmer/Desktop/biocpp/examples /home/zimmer/Desktop/biocpp/examples /home/zimmer/Desktop/biocpp/examples /home/zimmer/Desktop/biocpp/examples /home/zimmer/Desktop/biocpp/examples/CMakeFiles/iterate.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/iterate.dir/depend

