# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
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

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/kt-fitz/human-ergodicObj

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kt-fitz/human-ergodicObj

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/kt-fitz/human-ergodicObj/CMakeFiles /home/kt-fitz/human-ergodicObj/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/kt-fitz/human-ergodicObj/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named mda_cartpend

# Build rule for target.
mda_cartpend: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 mda_cartpend
.PHONY : mda_cartpend

# fast build rule for target.
mda_cartpend/fast:
	$(MAKE) -f CMakeFiles/mda_cartpend.dir/build.make CMakeFiles/mda_cartpend.dir/build
.PHONY : mda_cartpend/fast

#=============================================================================
# Target rules for targets named mig_mda_cartpend

# Build rule for target.
mig_mda_cartpend: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 mig_mda_cartpend
.PHONY : mig_mda_cartpend

# fast build rule for target.
mig_mda_cartpend/fast:
	$(MAKE) -f CMakeFiles/mig_mda_cartpend.dir/build.make CMakeFiles/mig_mda_cartpend.dir/build
.PHONY : mig_mda_cartpend/fast

#=============================================================================
# Target rules for targets named sac_cartpend

# Build rule for target.
sac_cartpend: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 sac_cartpend
.PHONY : sac_cartpend

# fast build rule for target.
sac_cartpend/fast:
	$(MAKE) -f CMakeFiles/sac_cartpend.dir/build.make CMakeFiles/sac_cartpend.dir/build
.PHONY : sac_cartpend/fast

mda_cartpend.o: mda_cartpend.cpp.o

.PHONY : mda_cartpend.o

# target to build an object file
mda_cartpend.cpp.o:
	$(MAKE) -f CMakeFiles/mda_cartpend.dir/build.make CMakeFiles/mda_cartpend.dir/mda_cartpend.cpp.o
.PHONY : mda_cartpend.cpp.o

mda_cartpend.i: mda_cartpend.cpp.i

.PHONY : mda_cartpend.i

# target to preprocess a source file
mda_cartpend.cpp.i:
	$(MAKE) -f CMakeFiles/mda_cartpend.dir/build.make CMakeFiles/mda_cartpend.dir/mda_cartpend.cpp.i
.PHONY : mda_cartpend.cpp.i

mda_cartpend.s: mda_cartpend.cpp.s

.PHONY : mda_cartpend.s

# target to generate assembly for a file
mda_cartpend.cpp.s:
	$(MAKE) -f CMakeFiles/mda_cartpend.dir/build.make CMakeFiles/mda_cartpend.dir/mda_cartpend.cpp.s
.PHONY : mda_cartpend.cpp.s

mig_mda_cartpend.o: mig_mda_cartpend.cpp.o

.PHONY : mig_mda_cartpend.o

# target to build an object file
mig_mda_cartpend.cpp.o:
	$(MAKE) -f CMakeFiles/mig_mda_cartpend.dir/build.make CMakeFiles/mig_mda_cartpend.dir/mig_mda_cartpend.cpp.o
.PHONY : mig_mda_cartpend.cpp.o

mig_mda_cartpend.i: mig_mda_cartpend.cpp.i

.PHONY : mig_mda_cartpend.i

# target to preprocess a source file
mig_mda_cartpend.cpp.i:
	$(MAKE) -f CMakeFiles/mig_mda_cartpend.dir/build.make CMakeFiles/mig_mda_cartpend.dir/mig_mda_cartpend.cpp.i
.PHONY : mig_mda_cartpend.cpp.i

mig_mda_cartpend.s: mig_mda_cartpend.cpp.s

.PHONY : mig_mda_cartpend.s

# target to generate assembly for a file
mig_mda_cartpend.cpp.s:
	$(MAKE) -f CMakeFiles/mig_mda_cartpend.dir/build.make CMakeFiles/mig_mda_cartpend.dir/mig_mda_cartpend.cpp.s
.PHONY : mig_mda_cartpend.cpp.s

sac_cartpend.o: sac_cartpend.cpp.o

.PHONY : sac_cartpend.o

# target to build an object file
sac_cartpend.cpp.o:
	$(MAKE) -f CMakeFiles/sac_cartpend.dir/build.make CMakeFiles/sac_cartpend.dir/sac_cartpend.cpp.o
.PHONY : sac_cartpend.cpp.o

sac_cartpend.i: sac_cartpend.cpp.i

.PHONY : sac_cartpend.i

# target to preprocess a source file
sac_cartpend.cpp.i:
	$(MAKE) -f CMakeFiles/sac_cartpend.dir/build.make CMakeFiles/sac_cartpend.dir/sac_cartpend.cpp.i
.PHONY : sac_cartpend.cpp.i

sac_cartpend.s: sac_cartpend.cpp.s

.PHONY : sac_cartpend.s

# target to generate assembly for a file
sac_cartpend.cpp.s:
	$(MAKE) -f CMakeFiles/sac_cartpend.dir/build.make CMakeFiles/sac_cartpend.dir/sac_cartpend.cpp.s
.PHONY : sac_cartpend.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... mda_cartpend"
	@echo "... mig_mda_cartpend"
	@echo "... sac_cartpend"
	@echo "... rebuild_cache"
	@echo "... mda_cartpend.o"
	@echo "... mda_cartpend.i"
	@echo "... mda_cartpend.s"
	@echo "... mig_mda_cartpend.o"
	@echo "... mig_mda_cartpend.i"
	@echo "... mig_mda_cartpend.s"
	@echo "... sac_cartpend.o"
	@echo "... sac_cartpend.i"
	@echo "... sac_cartpend.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

