CXX			= gcc
BIN_NAME	= mmult
CXXFLAGS	= -std=c99 -Wall -Wextra -fopenmp -lpapi -L/share/apps/papi/5.4.1/lib -I/share/apps/papi/5.4.1/include
####################################################################################
SRC_DIR		= src
BIN_DIR		= bin
BUILD_DIR	= build
SRC			= $(wildcard $(SRC_DIR)/*.c)
OBJ			= $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRC))
DEPS		= $(patsubst $(BUILD_DIR)/%.o,$(BUILD_DIR)/%.d,$(OBJ))
vpath %.c $(SRC_DIR) 
####################################################################################

.DEFAULT_GOAL = all

$(BUILD_DIR)/%.d: %.c
	$(CXX) -M $(CXXFLAGS) $< -o $@

$(BUILD_DIR)/%.o: %.c
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(BIN_DIR)/$(BIN_NAME): $(DEPS) $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)

all: $(BIN_DIR)/$(BIN_NAME)

clean:
	rm -f $(BUILD_DIR)/* $(BIN_DIR)/*
	clear
