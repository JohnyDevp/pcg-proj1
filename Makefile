# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/johnny/Code/MIT/pcg-proj1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/johnny/Code/MIT/pcg-proj1

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/johnny/Code/MIT/pcg-proj1/CMakeFiles /home/johnny/Code/MIT/pcg-proj1//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/johnny/Code/MIT/pcg-proj1/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named gen

# Build rule for target.
gen: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 gen
.PHONY : gen

# fast build rule for target.
gen/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/gen.dir/build.make CMakeFiles/gen.dir/build
.PHONY : gen/fast

#=============================================================================
# Target rules for targets named h5Helper

# Build rule for target.
h5Helper: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 h5Helper
.PHONY : h5Helper

# fast build rule for target.
h5Helper/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/h5Helper.dir/build.make CMakeFiles/h5Helper.dir/build
.PHONY : h5Helper/fast

#=============================================================================
# Target rules for targets named nbodyCpu

# Build rule for target.
nbodyCpu: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 nbodyCpu
.PHONY : nbodyCpu

# fast build rule for target.
nbodyCpu/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbodyCpu.dir/build.make CMakeFiles/nbodyCpu.dir/build
.PHONY : nbodyCpu/fast

#=============================================================================
# Target rules for targets named nbody0

# Build rule for target.
nbody0: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 nbody0
.PHONY : nbody0

# fast build rule for target.
nbody0/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody0.dir/build.make CMakeFiles/nbody0.dir/build
.PHONY : nbody0/fast

#=============================================================================
# Target rules for targets named nbody1

# Build rule for target.
nbody1: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 nbody1
.PHONY : nbody1

# fast build rule for target.
nbody1/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody1.dir/build.make CMakeFiles/nbody1.dir/build
.PHONY : nbody1/fast

#=============================================================================
# Target rules for targets named nbody2

# Build rule for target.
nbody2: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 nbody2
.PHONY : nbody2

# fast build rule for target.
nbody2/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody2.dir/build.make CMakeFiles/nbody2.dir/build
.PHONY : nbody2/fast

#=============================================================================
# Target rules for targets named nbody3

# Build rule for target.
nbody3: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 nbody3
.PHONY : nbody3

# fast build rule for target.
nbody3/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody3.dir/build.make CMakeFiles/nbody3.dir/build
.PHONY : nbody3/fast

#=============================================================================
# Target rules for targets named nbody4

# Build rule for target.
nbody4: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 nbody4
.PHONY : nbody4

# fast build rule for target.
nbody4/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody4.dir/build.make CMakeFiles/nbody4.dir/build
.PHONY : nbody4/fast

Commons/gen.o: Commons/gen.cpp.o
.PHONY : Commons/gen.o

# target to build an object file
Commons/gen.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/gen.dir/build.make CMakeFiles/gen.dir/Commons/gen.cpp.o
.PHONY : Commons/gen.cpp.o

Commons/gen.i: Commons/gen.cpp.i
.PHONY : Commons/gen.i

# target to preprocess a source file
Commons/gen.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/gen.dir/build.make CMakeFiles/gen.dir/Commons/gen.cpp.i
.PHONY : Commons/gen.cpp.i

Commons/gen.s: Commons/gen.cpp.s
.PHONY : Commons/gen.s

# target to generate assembly for a file
Commons/gen.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/gen.dir/build.make CMakeFiles/gen.dir/Commons/gen.cpp.s
.PHONY : Commons/gen.cpp.s

Commons/h5Helper.o: Commons/h5Helper.cpp.o
.PHONY : Commons/h5Helper.o

# target to build an object file
Commons/h5Helper.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/h5Helper.dir/build.make CMakeFiles/h5Helper.dir/Commons/h5Helper.cpp.o
.PHONY : Commons/h5Helper.cpp.o

Commons/h5Helper.i: Commons/h5Helper.cpp.i
.PHONY : Commons/h5Helper.i

# target to preprocess a source file
Commons/h5Helper.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/h5Helper.dir/build.make CMakeFiles/h5Helper.dir/Commons/h5Helper.cpp.i
.PHONY : Commons/h5Helper.cpp.i

Commons/h5Helper.s: Commons/h5Helper.cpp.s
.PHONY : Commons/h5Helper.s

# target to generate assembly for a file
Commons/h5Helper.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/h5Helper.dir/build.make CMakeFiles/h5Helper.dir/Commons/h5Helper.cpp.s
.PHONY : Commons/h5Helper.cpp.s

Cpu/main.o: Cpu/main.cpp.o
.PHONY : Cpu/main.o

# target to build an object file
Cpu/main.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbodyCpu.dir/build.make CMakeFiles/nbodyCpu.dir/Cpu/main.cpp.o
.PHONY : Cpu/main.cpp.o

Cpu/main.i: Cpu/main.cpp.i
.PHONY : Cpu/main.i

# target to preprocess a source file
Cpu/main.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbodyCpu.dir/build.make CMakeFiles/nbodyCpu.dir/Cpu/main.cpp.i
.PHONY : Cpu/main.cpp.i

Cpu/main.s: Cpu/main.cpp.s
.PHONY : Cpu/main.s

# target to generate assembly for a file
Cpu/main.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbodyCpu.dir/build.make CMakeFiles/nbodyCpu.dir/Cpu/main.cpp.s
.PHONY : Cpu/main.cpp.s

Cpu/nbody.o: Cpu/nbody.cpp.o
.PHONY : Cpu/nbody.o

# target to build an object file
Cpu/nbody.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbodyCpu.dir/build.make CMakeFiles/nbodyCpu.dir/Cpu/nbody.cpp.o
.PHONY : Cpu/nbody.cpp.o

Cpu/nbody.i: Cpu/nbody.cpp.i
.PHONY : Cpu/nbody.i

# target to preprocess a source file
Cpu/nbody.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbodyCpu.dir/build.make CMakeFiles/nbodyCpu.dir/Cpu/nbody.cpp.i
.PHONY : Cpu/nbody.cpp.i

Cpu/nbody.s: Cpu/nbody.cpp.s
.PHONY : Cpu/nbody.s

# target to generate assembly for a file
Cpu/nbody.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbodyCpu.dir/build.make CMakeFiles/nbodyCpu.dir/Cpu/nbody.cpp.s
.PHONY : Cpu/nbody.cpp.s

Step0/main.o: Step0/main.cu.o
.PHONY : Step0/main.o

# target to build an object file
Step0/main.cu.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody0.dir/build.make CMakeFiles/nbody0.dir/Step0/main.cu.o
.PHONY : Step0/main.cu.o

Step0/main.i: Step0/main.cu.i
.PHONY : Step0/main.i

# target to preprocess a source file
Step0/main.cu.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody0.dir/build.make CMakeFiles/nbody0.dir/Step0/main.cu.i
.PHONY : Step0/main.cu.i

Step0/main.s: Step0/main.cu.s
.PHONY : Step0/main.s

# target to generate assembly for a file
Step0/main.cu.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody0.dir/build.make CMakeFiles/nbody0.dir/Step0/main.cu.s
.PHONY : Step0/main.cu.s

Step0/nbody.o: Step0/nbody.cu.o
.PHONY : Step0/nbody.o

# target to build an object file
Step0/nbody.cu.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody0.dir/build.make CMakeFiles/nbody0.dir/Step0/nbody.cu.o
.PHONY : Step0/nbody.cu.o

Step0/nbody.i: Step0/nbody.cu.i
.PHONY : Step0/nbody.i

# target to preprocess a source file
Step0/nbody.cu.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody0.dir/build.make CMakeFiles/nbody0.dir/Step0/nbody.cu.i
.PHONY : Step0/nbody.cu.i

Step0/nbody.s: Step0/nbody.cu.s
.PHONY : Step0/nbody.s

# target to generate assembly for a file
Step0/nbody.cu.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody0.dir/build.make CMakeFiles/nbody0.dir/Step0/nbody.cu.s
.PHONY : Step0/nbody.cu.s

Step1/main.o: Step1/main.cu.o
.PHONY : Step1/main.o

# target to build an object file
Step1/main.cu.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody1.dir/build.make CMakeFiles/nbody1.dir/Step1/main.cu.o
.PHONY : Step1/main.cu.o

Step1/main.i: Step1/main.cu.i
.PHONY : Step1/main.i

# target to preprocess a source file
Step1/main.cu.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody1.dir/build.make CMakeFiles/nbody1.dir/Step1/main.cu.i
.PHONY : Step1/main.cu.i

Step1/main.s: Step1/main.cu.s
.PHONY : Step1/main.s

# target to generate assembly for a file
Step1/main.cu.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody1.dir/build.make CMakeFiles/nbody1.dir/Step1/main.cu.s
.PHONY : Step1/main.cu.s

Step1/nbody.o: Step1/nbody.cu.o
.PHONY : Step1/nbody.o

# target to build an object file
Step1/nbody.cu.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody1.dir/build.make CMakeFiles/nbody1.dir/Step1/nbody.cu.o
.PHONY : Step1/nbody.cu.o

Step1/nbody.i: Step1/nbody.cu.i
.PHONY : Step1/nbody.i

# target to preprocess a source file
Step1/nbody.cu.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody1.dir/build.make CMakeFiles/nbody1.dir/Step1/nbody.cu.i
.PHONY : Step1/nbody.cu.i

Step1/nbody.s: Step1/nbody.cu.s
.PHONY : Step1/nbody.s

# target to generate assembly for a file
Step1/nbody.cu.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody1.dir/build.make CMakeFiles/nbody1.dir/Step1/nbody.cu.s
.PHONY : Step1/nbody.cu.s

Step2/main.o: Step2/main.cu.o
.PHONY : Step2/main.o

# target to build an object file
Step2/main.cu.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody2.dir/build.make CMakeFiles/nbody2.dir/Step2/main.cu.o
.PHONY : Step2/main.cu.o

Step2/main.i: Step2/main.cu.i
.PHONY : Step2/main.i

# target to preprocess a source file
Step2/main.cu.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody2.dir/build.make CMakeFiles/nbody2.dir/Step2/main.cu.i
.PHONY : Step2/main.cu.i

Step2/main.s: Step2/main.cu.s
.PHONY : Step2/main.s

# target to generate assembly for a file
Step2/main.cu.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody2.dir/build.make CMakeFiles/nbody2.dir/Step2/main.cu.s
.PHONY : Step2/main.cu.s

Step2/nbody.o: Step2/nbody.cu.o
.PHONY : Step2/nbody.o

# target to build an object file
Step2/nbody.cu.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody2.dir/build.make CMakeFiles/nbody2.dir/Step2/nbody.cu.o
.PHONY : Step2/nbody.cu.o

Step2/nbody.i: Step2/nbody.cu.i
.PHONY : Step2/nbody.i

# target to preprocess a source file
Step2/nbody.cu.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody2.dir/build.make CMakeFiles/nbody2.dir/Step2/nbody.cu.i
.PHONY : Step2/nbody.cu.i

Step2/nbody.s: Step2/nbody.cu.s
.PHONY : Step2/nbody.s

# target to generate assembly for a file
Step2/nbody.cu.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody2.dir/build.make CMakeFiles/nbody2.dir/Step2/nbody.cu.s
.PHONY : Step2/nbody.cu.s

Step3/main.o: Step3/main.cu.o
.PHONY : Step3/main.o

# target to build an object file
Step3/main.cu.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody3.dir/build.make CMakeFiles/nbody3.dir/Step3/main.cu.o
.PHONY : Step3/main.cu.o

Step3/main.i: Step3/main.cu.i
.PHONY : Step3/main.i

# target to preprocess a source file
Step3/main.cu.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody3.dir/build.make CMakeFiles/nbody3.dir/Step3/main.cu.i
.PHONY : Step3/main.cu.i

Step3/main.s: Step3/main.cu.s
.PHONY : Step3/main.s

# target to generate assembly for a file
Step3/main.cu.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody3.dir/build.make CMakeFiles/nbody3.dir/Step3/main.cu.s
.PHONY : Step3/main.cu.s

Step3/nbody.o: Step3/nbody.cu.o
.PHONY : Step3/nbody.o

# target to build an object file
Step3/nbody.cu.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody3.dir/build.make CMakeFiles/nbody3.dir/Step3/nbody.cu.o
.PHONY : Step3/nbody.cu.o

Step3/nbody.i: Step3/nbody.cu.i
.PHONY : Step3/nbody.i

# target to preprocess a source file
Step3/nbody.cu.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody3.dir/build.make CMakeFiles/nbody3.dir/Step3/nbody.cu.i
.PHONY : Step3/nbody.cu.i

Step3/nbody.s: Step3/nbody.cu.s
.PHONY : Step3/nbody.s

# target to generate assembly for a file
Step3/nbody.cu.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody3.dir/build.make CMakeFiles/nbody3.dir/Step3/nbody.cu.s
.PHONY : Step3/nbody.cu.s

Step4/main.o: Step4/main.cu.o
.PHONY : Step4/main.o

# target to build an object file
Step4/main.cu.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody4.dir/build.make CMakeFiles/nbody4.dir/Step4/main.cu.o
.PHONY : Step4/main.cu.o

Step4/main.i: Step4/main.cu.i
.PHONY : Step4/main.i

# target to preprocess a source file
Step4/main.cu.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody4.dir/build.make CMakeFiles/nbody4.dir/Step4/main.cu.i
.PHONY : Step4/main.cu.i

Step4/main.s: Step4/main.cu.s
.PHONY : Step4/main.s

# target to generate assembly for a file
Step4/main.cu.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody4.dir/build.make CMakeFiles/nbody4.dir/Step4/main.cu.s
.PHONY : Step4/main.cu.s

Step4/nbody.o: Step4/nbody.cu.o
.PHONY : Step4/nbody.o

# target to build an object file
Step4/nbody.cu.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody4.dir/build.make CMakeFiles/nbody4.dir/Step4/nbody.cu.o
.PHONY : Step4/nbody.cu.o

Step4/nbody.i: Step4/nbody.cu.i
.PHONY : Step4/nbody.i

# target to preprocess a source file
Step4/nbody.cu.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody4.dir/build.make CMakeFiles/nbody4.dir/Step4/nbody.cu.i
.PHONY : Step4/nbody.cu.i

Step4/nbody.s: Step4/nbody.cu.s
.PHONY : Step4/nbody.s

# target to generate assembly for a file
Step4/nbody.cu.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/nbody4.dir/build.make CMakeFiles/nbody4.dir/Step4/nbody.cu.s
.PHONY : Step4/nbody.cu.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... gen"
	@echo "... h5Helper"
	@echo "... nbody0"
	@echo "... nbody1"
	@echo "... nbody2"
	@echo "... nbody3"
	@echo "... nbody4"
	@echo "... nbodyCpu"
	@echo "... Commons/gen.o"
	@echo "... Commons/gen.i"
	@echo "... Commons/gen.s"
	@echo "... Commons/h5Helper.o"
	@echo "... Commons/h5Helper.i"
	@echo "... Commons/h5Helper.s"
	@echo "... Cpu/main.o"
	@echo "... Cpu/main.i"
	@echo "... Cpu/main.s"
	@echo "... Cpu/nbody.o"
	@echo "... Cpu/nbody.i"
	@echo "... Cpu/nbody.s"
	@echo "... Step0/main.o"
	@echo "... Step0/main.i"
	@echo "... Step0/main.s"
	@echo "... Step0/nbody.o"
	@echo "... Step0/nbody.i"
	@echo "... Step0/nbody.s"
	@echo "... Step1/main.o"
	@echo "... Step1/main.i"
	@echo "... Step1/main.s"
	@echo "... Step1/nbody.o"
	@echo "... Step1/nbody.i"
	@echo "... Step1/nbody.s"
	@echo "... Step2/main.o"
	@echo "... Step2/main.i"
	@echo "... Step2/main.s"
	@echo "... Step2/nbody.o"
	@echo "... Step2/nbody.i"
	@echo "... Step2/nbody.s"
	@echo "... Step3/main.o"
	@echo "... Step3/main.i"
	@echo "... Step3/main.s"
	@echo "... Step3/nbody.o"
	@echo "... Step3/nbody.i"
	@echo "... Step3/nbody.s"
	@echo "... Step4/main.o"
	@echo "... Step4/main.i"
	@echo "... Step4/main.s"
	@echo "... Step4/nbody.o"
	@echo "... Step4/nbody.i"
	@echo "... Step4/nbody.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

