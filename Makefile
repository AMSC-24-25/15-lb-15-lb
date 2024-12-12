
INC_DIR = include
SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = build

CC = g++
CFLAGS = -Wall -I$(INC_DIR) -fopenmp

TARGET = $(BIN_DIR)/lbm

SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

all: $(TARGET)
	./$(TARGET)

$(TARGET): $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $^

$(OBJ_FILES): $(SRC_FILES)
	$(CC) $(CFLAGS) -c $< -o $@


