# Basics
SHELL = /bin/sh
REPORT = report/main.tex
NSIZE = 32
MSIZE = 32

# C build setup
CC = gcc
SRCDIR = src
BUILD_DIR = build
BIN_DIR = bin
TEST_DIR = test
TEST_TARGET = bin/test
TARGET = bin/mq.so

# Files
SRCEXT := c
SRCS :=	$(shell find $(SRCDIR)/c -type f -name *.$(SRCEXT))
OBJ := $(patsubst $(SRCDIR)/c/%,$(BUILD_DIR)/%,$(SRCS:.$(SRCEXT)=.o))
TEST_SRCS := $(shell find $(TEST_DIR) -type f -name *.$(SRCEXT))
TEST_OBJ := $(patsubst $(TEST_DIR)/%,$(BUILD_DIR)/%,$(TEST_SRCS:.$(SRCEXT)=.o)) $(OBJ)

# GCC flags
SAN :=\
			 -fsanitize=address \
			 -fsanitize=leak \
			 -fsanitize=undefined
DEBUG := $(SAN) -Wall -Wextra -O0 -g -D_DEBUG
OPT := -O3
LDFLAGS := -lm
INC := -Iinc
CFLAGS := $(INC) -DP$(MSIZE) -DV$(NSIZE) $(LIB)

VPATH = src/c:test

all: $(TARGET)

$(BUILD_DIR)/%.o: %.c binom.h sage | $(BUILD_DIR) $(BIN_DIR)
	$(CC) -o $@ $(CFLAGS) -c $<

$(TARGET): CFLAGS += $(OPT) -fPIC
$(TARGET): $(OBJ)
	$(CC) -shared -o $@ $^

tests: CFLAGS += $(DEBUG)
tests: LDFLAGS += -fsanitize=address -fsanitize=undefined -fsanitize=leak
tests: $(TEST_TARGET)

$(TEST_TARGET): $(TEST_OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

sage:
	make -C src/sage 

binom.h:
	./binom.py $(NSIZE) $(MSIZE) > inc/binom.h

pdf:
	latexmk -pdf -silent -cd $(REPORT)

pdfclean:
	latexmk -C -cd $(REPORT)
	rm *.bbl *.xml

clean:
	$(RM) -rf $(BUILD_DIR) $(BIN_DIR)
	make -C src/sage/ clean

$(BIN_DIR):
	mkdir -p $@

$(BUILD_DIR):
	mkdir -p $@

.PHONY: clean all pdf tests pdfclean
