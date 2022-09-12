# Basics
SHELL = /bin/sh
REPORT = report/main.tex

# C build setup
CC = gcc
SRCDIR = src
BUILD_DIR = build
TEST_DIR = test
TEST_TARGET = bin/test.elf
TARGET = bin/mq.elf

# Files
SRCEXT := c
SRCS :=	$(shell find $(SRCDIR)/c -type f -name *.$(SRCEXT))
OBJ := $(patsubst $(SRCDIR)/c/%,$(BUILD_DIR)/%,$(SRCS:.$(SRCEXT)=.o))
TEST_SRCS := $(shell find $(TEST_DIR) -type f -name *.$(SRCEXT)) $(filter-out src/c/mq.c,$(SRCS))
TEST_OBJ := $(patsubst $(TEST_DIR)/%,$(BUILD_DIR)/%,$(TEST_SRCS:.$(SRCEXT)=.o))

# GCC flags
SAN :=\
			 -fsanitize=address \
			 -fsanitize=leak \
			 -fsanitize=undefined
DEBUG := $(SAN) -Wall -Wextra -O0
OPT := -O3 -Os
CFLAGS := -g 
LIB := 
INC := -Iinc

VPATH = src/c:test

# Targets
all: $(TARGET)

$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) -o $@ -c $(CFLAGS) $<

$(TARGET): $(OBJ)
	@echo "Not written yet"

tests: CFLAGS += $(DEBUG)
tests: $(TEST_TARGET) 

$(TEST_TARGET): $(TEST_OBJ)
	@echo "Not written yet"

pdf:
	latexmk -pdf -silent -cd $(REPORT)

pdfclean:
	latexmk -C -cd $(REPORT)

clean:
	$(RM) -rf $(BUILD_DIR) $(TARGET) $(TEST_TARGET) bin

$(BUILD_DIR):
	mkdir -p $@

.PHONY: clean all pdf tests pdfclean
