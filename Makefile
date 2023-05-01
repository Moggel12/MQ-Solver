# Basics
SHELL = /bin/sh
REPORT = report/main.tex
REGSIZE = 32


# C build setup
CC = gcc
SRCDIR = src
BUILD_DIR = build
BIN_DIR = bin
TEST_DIR = test
STD_DIR = $(BUILD_DIR)/standard
VEC_DIR = $(BUILD_DIR)/vectorzied
TEST_TARGET = bin/test
TARGET = bin/mq.so
VEC_TARGET = bin/mq_vec.so

# Files
SRCEXT := c
ifeq ($(REGSIZE), 128)
WORKDIR := $(SRCDIR)/c/vectorized
SRCS := $(shell find $(SRCDIR)/c/vectorized -type f -name '*.$(SRCEXT)')
INTSIZE := 16
CFLAGS += -mavx -DINT$(INTSIZE)
$(shell ./binom.py 64 64 > inc/binom.h)
else ifeq ($(REGSIZE), 256)
WORKDIR := $(SRCDIR)/c/vectorized
SRCS := $(shell find $(SRCDIR)/c/vectorized -type f -name '*.$(SRCEXT)')
INTSIZE := 16
CFLAGS += -mavx2 -DINT$(INTSIZE)
$(shell ./binom.py 64 64 > inc/binom.h)
else 
WORKDIR := $(SRCDIR)/c/standard
SRCS :=	$(shell find $(SRCDIR)/c/standard -type f -name '*.$(SRCEXT)')
INTSIZE := $(REGSIZE)
$(shell ./binom.py $(INTSIZE) $(INTSIZE) > inc/binom.h)
endif

SRCS += $(patsubst %,$(SRCDIR)/c/%,$(shell ls -pA src/c | grep -v /))

TEST_SRCS := $(shell find $(TEST_DIR) -type f -name *.$(SRCEXT))
OBJ := $(patsubst $(WORKDIR)/%,$(BUILD_DIR)/%,$(filter $(WORKDIR)/%, $(SRCS:.$(SRCEXT)=.o))) $(patsubst $(SRCDIR)/c/%,$(BUILD_DIR)/%, $(filter-out $(WORKDIR)/%, $(SRCS:.$(SRCEXT)=.o)))
TEST_OBJ := $(patsubst $(TEST_DIR)/%,$(BUILD_DIR)/%,$(TEST_SRCS:.$(SRCEXT)=.o)) $(OBJ)

$(shell echo "$(REGSIZE) $(INTSIZE)" > src/sage/.compile_config)

# GCC flags
SAN :=\
			 -fsanitize=address \
			 -fsanitize=leak \
			 -fsanitize=undefined
DEBUG := $(SAN) -Wall -Wextra -O0 -g -D_DEBUG
OPT := -O3
LDFLAGS := -lm
INC := -Iinc
CFLAGS += $(INC) -DREG$(REGSIZE) $(LIB)

VPATH = src/c:$(WORKDIR):test

$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR) $(BIN_DIR)
	$(CC) -o $@ $(CFLAGS) -c $<

$(TARGET): CFLAGS += $(OPT) -fPIC
$(TARGET): $(OBJ) | sage $(STD_DIR)
	$(CC) -shared -o $@ $^

tests: CFLAGS += $(DEBUG)
tests: LDFLAGS += -fsanitize=address -fsanitize=undefined -fsanitize=leak
tests: $(TEST_TARGET) | sage

$(TEST_TARGET): $(TEST_OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

sage:
	make -C src/sage 

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

$(STD_DIR):
	mkdir -p $(BUILD_DIR)/$(STD_DIR)

$(VEC_DIR):
	mkdir -p $(BUILD_DIR)/$(VEC_DIR)

.PHONY: clean standard vector pdf tests pdfclean sage
