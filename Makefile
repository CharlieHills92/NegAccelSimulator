# Simple Makefile for ibsimu simulation project
# Uses main.cpp directly with NEW_IMPLEMENTATION

# Compiler and tools
CXX = g++
CXX_STD = c++11
RM = rm -f
MKDIR_P = mkdir -p

# Directories
BUILD_DIR = build
IBSIMU_ROOT = $(CURDIR)

# Use environment PKG_CONFIG_PATH if set, otherwise use default
PKG_CONFIG_PATH ?= $(IBSIMU_ROOT)/lib/pkgconfig
export PKG_CONFIG_PATH

# Compiler flags
CXXFLAGS = -std=$(CXX_STD) -Wall -Wextra -g -O2 -MMD -MP -DNEW_IMPLEMENTATION

# Always add ibsimu include path directly
CXXFLAGS += -I$(IBSIMU_ROOT)/include/ibsimu-1.0.6dev

# Library flags (try pkg-config first, fallback to direct paths)
LIBS := $(shell pkg-config --libs ibsimu-1.0.6dev 2>/dev/null || echo "-L$(IBSIMU_ROOT)/lib -libsimu-1.0.6dev -lcairo -lglib-2.0 -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas -llapack -lgsl -lgslcblas -lm")

# Include flags from pkg-config if available
CXXFLAGS += $(shell pkg-config --cflags ibsimu-1.0.6dev 2>/dev/null || echo "")

# Add cairo, GTK3, and glib include flags (needed by ibsimu headers)
CXXFLAGS += $(shell pkg-config --cflags cairo 2>/dev/null || echo "")
CXXFLAGS += $(shell pkg-config --cflags gtk+-3.0 2>/dev/null || echo "")
CXXFLAGS += $(shell pkg-config --cflags glib-2.0 2>/dev/null || echo "")

# Linker flags
LDFLAGS = -Wall -g

# Source files for main simulation (using refactored managers)
SOURCES = main.cpp ManageSimulation_New.cpp \
          SimulationParameters.cpp FileManager.cpp GeometryManager.cpp \
          FieldManager.cpp ParticleManager.cpp DiagnosticsManager.cpp \
          ScanManager.cpp \
          TransverseData.cpp my_diagnostics.cpp globals.cpp funct.cpp \
          cross_sections.cpp THCallback.cpp StrippingUtils.cpp

# Object files  
OBJECTS = $(SOURCES:%.cpp=$(BUILD_DIR)/%.o)

# Dependency files
DEPS = $(OBJECTS:.o=.d)

# Default target
all: runtest_new_v2

# Build main executable
runtest_new_v2: $(OBJECTS) | $(BUILD_DIR)
	@echo "Linking runtest_new_v2..."
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS) $(LIBS)
	@echo "Successfully built runtest_new_v2"

# Object file compilation
$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create build directory
$(BUILD_DIR):
	$(MKDIR_P) $(BUILD_DIR)

# Include dependency files
-include $(DEPS)

# Clean targets
clean:
	$(RM) $(BUILD_DIR)/*.o $(BUILD_DIR)/*.d

distclean: clean
	$(RM) runtest_new_v2
	$(RM) -r $(BUILD_DIR)

# Test targets
test: runtest_new_v2
	@echo "Running basic test..."
	./runtest_new_v2

# Test for grid power analysis
TEST_OBJECTS = $(filter-out $(BUILD_DIR)/main.o, $(OBJECTS))
test_grid_power: $(TEST_OBJECTS) test_grid_power.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

# Test for MTF grid power analysis
test_mtf_grid_power: $(TEST_OBJECTS) test_mtf_grid_power.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

# Test for loading and tracing particles with secondary generation
test_load_and_trace: $(TEST_OBJECTS) test_load_and_trace.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

.PHONY: all clean distclean test test_grid_power test_mtf_grid_power test_load_and_trace

# Help target
help:
	@echo "Available targets:"
	@echo "  all         - Build runtest_new_v2 executable (default)"
	@echo "  runtest_new_v2 - Build the simulation executable using main.cpp"
	@echo "  clean       - Remove object and dependency files"
	@echo "  distclean   - Remove all generated files"
	@echo "  test        - Build and run basic test"
	@echo "  test_grid_power - Build and run grid power analysis test"
	@echo "  test_load_and_trace - Build and run field loading and particle tracing test"
	@echo "  help        - Show this help message"
