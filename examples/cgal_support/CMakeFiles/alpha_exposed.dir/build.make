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
CMAKE_SOURCE_DIR = /home/zimmer/Desktop/biocpp/examples/cgal_support

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zimmer/Desktop/biocpp/examples/cgal_support

# Include any dependencies generated for this target.
include CMakeFiles/alpha_exposed.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/alpha_exposed.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/alpha_exposed.dir/flags.make

CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o: CMakeFiles/alpha_exposed.dir/flags.make
CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o: alpha_exposed.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/zimmer/Desktop/biocpp/examples/cgal_support/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) -DBIOCPP_INCLUDE_PDB -DBIOCPP_INCLUDE_STANDARD -DBIOCPP_INCLUDE_CGAL $(CXX_FLAGS) -o CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o -c /home/zimmer/Desktop/biocpp/examples/cgal_support/alpha_exposed.cpp

CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) -DBIOCPP_INCLUDE_PDB -DBIOCPP_INCLUDE_STANDARD -DBIOCPP_INCLUDE_CGAL $(CXX_FLAGS) -E /home/zimmer/Desktop/biocpp/examples/cgal_support/alpha_exposed.cpp > CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.i

CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) -DBIOCPP_INCLUDE_PDB -DBIOCPP_INCLUDE_STANDARD -DBIOCPP_INCLUDE_CGAL $(CXX_FLAGS) -S /home/zimmer/Desktop/biocpp/examples/cgal_support/alpha_exposed.cpp -o CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.s

CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o.requires:
.PHONY : CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o.requires

CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o.provides: CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o.requires
	$(MAKE) -f CMakeFiles/alpha_exposed.dir/build.make CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o.provides.build
.PHONY : CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o.provides

CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o.provides.build: CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o

# Object files for target alpha_exposed
alpha_exposed_OBJECTS = \
"CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o"

# External object files for target alpha_exposed
alpha_exposed_EXTERNAL_OBJECTS =

alpha_exposed: CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o
alpha_exposed: /usr/lib/libCGAL.so
alpha_exposed: /usr/lib/libgmpxx.so
alpha_exposed: /usr/lib/libmpfr.so
alpha_exposed: /usr/lib/libgmp.so
alpha_exposed: /usr/lib/libboost_thread-mt.so
alpha_exposed: /usr/lib/libCGAL.so
alpha_exposed: /usr/lib/libgmpxx.so
alpha_exposed: /usr/lib/libmpfr.so
alpha_exposed: /usr/lib/libgmp.so
alpha_exposed: /usr/lib/libboost_thread-mt.so
alpha_exposed: CMakeFiles/alpha_exposed.dir/build.make
alpha_exposed: CMakeFiles/alpha_exposed.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable alpha_exposed"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/alpha_exposed.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/alpha_exposed.dir/build: alpha_exposed
.PHONY : CMakeFiles/alpha_exposed.dir/build

CMakeFiles/alpha_exposed.dir/requires: CMakeFiles/alpha_exposed.dir/alpha_exposed.cpp.o.requires
.PHONY : CMakeFiles/alpha_exposed.dir/requires

CMakeFiles/alpha_exposed.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/alpha_exposed.dir/cmake_clean.cmake
.PHONY : CMakeFiles/alpha_exposed.dir/clean

CMakeFiles/alpha_exposed.dir/depend:
	cd /home/zimmer/Desktop/biocpp/examples/cgal_support && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zimmer/Desktop/biocpp/examples/cgal_support /home/zimmer/Desktop/biocpp/examples/cgal_support /home/zimmer/Desktop/biocpp/examples/cgal_support /home/zimmer/Desktop/biocpp/examples/cgal_support /home/zimmer/Desktop/biocpp/examples/cgal_support/CMakeFiles/alpha_exposed.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/alpha_exposed.dir/depend

