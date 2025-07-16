default: all

# host compiler
HOST_COMPILER ?= g++
# -fsycl 
# -fsycl-targets=spir64_gen
# -fdebug-info-for-profiling -gline-tables-only # icpx # g++


NAMESPACE_DIR := ../src/cpp

CPP_DIR_SRC := ../src/cpp
CPP_DIR_OBJ := ../obj/cpp
CPP_DIR_DEP := ../obj/cpp

CPP_SRCS_LOCAL := $(wildcard $(CPP_DIR_SRC)/*.cpp)
CPP_OBJS_LOCAL := $(patsubst $(CPP_DIR_SRC)/%.cpp,$(CPP_DIR_OBJ)/%.o,$(CPP_SRCS_LOCAL))
CPP_DEPS_LOCAL := $(patsubst $(CPP_DIR_SRC)/%.cpp,$(CPP_DIR_DEP)/%.d,$(CPP_SRCS_LOCAL))


CPP_DIR_SRC_REMOTE := ../test_src_remote/cpp
CPP_DIR_OBJ_REMOTE := ../test_obj_remote/cpp
CPP_DIR_DEP_REMOTE := ../test_obj_remote/cpp

CPP_SRCS_REMOTE := $(wildcard $(CPP_DIR_SRC_REMOTE)/*.cpp)
CPP_OBJS_REMOTE := $(patsubst $(CPP_DIR_SRC_REMOTE)/%.cpp,$(CPP_DIR_OBJ_REMOTE)/%.o,$(CPP_SRCS_REMOTE))
CPP_DEPS_REMOTE := $(patsubst $(CPP_DIR_SRC_REMOTE)/%.cpp,$(CPP_DIR_DEP_REMOTE)/%.d,$(CPP_SRCS_REMOTE))


CPP_DIR_SRC_MAIN := ../src/main_cpp
CPP_DIR_OBJ_MAIN := ../obj/main_cpp
CPP_DIR_DEP_MAIN := ../obj/main_cpp

CPP_SRCS_LOCAL_MAIN := $(wildcard $(CPP_DIR_SRC_MAIN)/*.cpp)
CPP_OBJS_LOCAL_MAIN := $(patsubst $(CPP_DIR_SRC_MAIN)/%.cpp,$(CPP_DIR_OBJ_MAIN)/%.o,$(CPP_SRCS_LOCAL_MAIN))
CPP_DEPS_LOCAL_MAIN := $(patsubst $(CPP_DIR_SRC_MAIN)/%.cpp,$(CPP_DIR_DEP_MAIN)/%.d,$(CPP_SRCS_LOCAL_MAIN))


INC_SUBDIR = # linear_memory small_matrix helper storage 
INC_FLAGS  = $(addprefix -I$(HOME)/github/common_headers/, $(INC_SUBDIR))

INC_FLAGS_HOST = -I$(CPP_DIR_SRC) -I$(CPP_DIR_SRC_REMOTE) -I$(CPP_DIR_SRC_MAIN) $(INC_FLAGS) -I$(NAMESPACE_DIR)

CPPFLAGS_HOST = -O3 -Wall -fmessage-length=0 -std=c++20 -march=native # -fdebug-info-for-profiling -gline-tables-only # -Wconversion -Warith-conversion
# [2024/02/03]
# NOTE: Below is for intel.
# CPPFLAGS_HOST = -O3 -Wall -fmessage-length=0 -std=c++20 -march=native -no-ftz -fp-model precise # -fdebug-info-for-profiling -gline-tables-only # -Wconversion -Warith-conversion
CPPFLAGS_HOST += $(INC_FLAGS_HOST)


MACROFLAGS =


# DPPFLAGS: flags for dependencies
DPPFLAGS = -MMD -MP -MF "$(@:%.o=%.d)" -MT "$(@)"
# NOTE: For some reason, we need a space after "-MF" and "-MT" for it to work with 
# 		nvcc, while with gcc it is fine without space.

LINKER_FLAGS = # -lstdc++fs # -Xsycl-target-backend "-device pvc"

# Additional compiler flag
CPPFLAGS_HOST += # -g  # -g is needed for using gdb


CPP_EXE_DIR := .
CPP_EXE_LIST := $(patsubst $(CPP_DIR_SRC_MAIN)/%.cpp,$(CPP_EXE_DIR)/%.exe,$(CPP_SRCS_LOCAL_MAIN))


all : $(CPP_EXE_LIST)
# 	@echo 'executables: $(CPP_EXE_LIST)'
# 	@echo 'object files: $(CPP_OBJS_LOCAL_MAIN)'


CPP_EXE_NAME := $(patsubst $(CPP_DIR_SRC_MAIN)/%.cpp,$(CPP_EXE_DIR)/%,$(CPP_SRCS_LOCAL_MAIN))
$(CPP_EXE_NAME) : $(CPP_EXE_DIR)/% : $(CPP_EXE_DIR)/%.exe


# [2022/12/10]
# NOTE: When linking against (static) libraries, be sure to place the (user compiled) 
# 		object files first. Otherwise, strange linking error may occur.

# Linking
$(CPP_EXE_LIST) : $(CPP_EXE_DIR)/%.exe : $(CPP_DIR_OBJ_MAIN)/%.o $(CPP_OBJS_LOCAL) $(CPP_OBJS_REMOTE)
	@echo 'Invoking linker for target: $@'
	@echo 'main object to be linked: $<'
	@echo ' '
	$(HOST_COMPILER) -o $@ $^ $(LINKER_FLAGS)
	@echo 'Finished building target: $@'
	@echo ' '

# compiling .cpp
$(CPP_OBJS_LOCAL_MAIN) : $(CPP_DIR_OBJ_MAIN)/%.o: $(CPP_DIR_SRC_MAIN)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking Compiler'
	$(HOST_COMPILER) $(CPPFLAGS_HOST) $(DPPFLAGS) $(MACROFLAGS) -o $@ -c $<
	@echo 'Finished building: $<'
	@echo ' '

$(CPP_DIR_OBJ)/%.o: $(CPP_DIR_SRC)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking Compiler'
	${HOST_COMPILER} $(CPPFLAGS_HOST) $(DPPFLAGS) $(MACROFLAGS) -o "$@" -c "$<"
	@echo 'Finished building: $<'
	@echo ' '

$(CPP_DIR_OBJ_REMOTE)/%.o: $(CPP_DIR_SRC_REMOTE)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking Compiler'
	${HOST_COMPILER} $(CPPFLAGS_HOST) $(DPPFLAGS) $(MACROFLAGS) -o "$@" -c "$<"
	@echo 'Finished building: $<'
	@echo ' '


-include $(CPP_OBJS_LOCAL:.o=.d) $(CPP_OBJS_LOCAL_MAIN:.o=.d) $(CPP_OBJS_REMOTE:.o=.d)

# Other Targets
RM := rm -rf

clm : clean_main
# [2023/01/31] NOTE: clm is shorthand for clean_main

clean_main::
	-$(RM) $(CPP_DEPS_LOCAL_MAIN) $(CPP_OBJS_LOCAL_MAIN) $(CPP_EXE_LIST)
	-@echo ' '

clean:
	-$(RM) $(CPP_DEPS_LOCAL) $(CPP_OBJS_LOCAL) \
		   $(CPP_DEPS_REMOTE) $(CPP_OBJS_REMOTE) \
		   $(CPP_DEPS_LOCAL_MAIN) $(CPP_OBJS_LOCAL_MAIN) \
		   $(CPP_EXE_LIST)
	-@echo ' '

.PHONY: all clean dependents

