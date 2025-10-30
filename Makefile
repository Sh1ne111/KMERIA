BIN_DIR := ./bin
SRC_DIR := src
OBJ_DIR := obj
INCLUDE_DIR := ./include
LIB_DIR := ./lib

#compiler
CC  := gcc
CXX := g++
CFLAGS := -Wall -O3 -I $(INCLUDE_DIR)
CXXFLAGS := $(CFLAGS) -std=c++14 -fPIE
LDFLAGS := -L $(LIB_DIR) -pie
LIBS := -lz -lhts -lpthread -lbz2
PROG := kmeria

ifneq ($(asan),)
    CFLAGS += -fsanitize=address
    CXXFLAGS += -fsanitize=address
    LDFLAGS += -fsanitize=address
endif

C_SOURCES := $(wildcard $(SRC_DIR)/*.c)
CPP_SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
SOURCES := $(C_SOURCES) $(CPP_SOURCES)

C_OBJS := $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(C_SOURCES))
CPP_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(CPP_SOURCES))
OBJS := $(C_OBJS) $(CPP_OBJS)

$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))

.PHONY: all clean install package

all: $(PROG)

# rule for linking
$(PROG): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -DUSE_BGZF -c $< -o $@

$(OBJ_DIR)/main.o: $(SRC_DIR)/main.c
$(OBJ_DIR)/kcount.o: $(SRC_DIR)/kcount.c $(SRC_DIR)/kstring.h $(SRC_DIR)/kthread.h $(SRC_DIR)/ketopt.h $(SRC_DIR)/kseq.h $(SRC_DIR)/kvec.h $(SRC_DIR)/ksort.h $(SRC_DIR)/util.h
$(OBJ_DIR)/kdump.o: $(SRC_DIR)/kdump.c $(SRC_DIR)/kdump.h $(SRC_DIR)/kstring.h
$(OBJ_DIR)/kmatrix.o: $(SRC_DIR)/kmatrix.cpp $(SRC_DIR)/kmatrix.h
$(OBJ_DIR)/kfilter.o: $(SRC_DIR)/kfilter.cpp $(SRC_DIR)/kfilter.h
$(OBJ_DIR)/kmtob.o: $(SRC_DIR)/kmtob.cpp $(SRC_DIR)/kmtob.h
$(OBJ_DIR)/kbtog.o: $(SRC_DIR)/kbtog.cpp $(SRC_DIR)/kbtog.h
$(OBJ_DIR)/ksketch.o: $(SRC_DIR)/ksketch.c
$(OBJ_DIR)/kassoc.o: $(SRC_DIR)/kassoc.cpp $(SRC_DIR)/kassoc.h
$(OBJ_DIR)/fkr.o: $(SRC_DIR)/fkr.c
$(OBJ_DIR)/fkrtgs.o: $(SRC_DIR)/fkrtgs.c
$(OBJ_DIR)/kmc_file.o: $(SRC_DIR)/kmc_file.cpp $(SRC_DIR)/kmc_file.h $(SRC_DIR)/kmer_defs.h
$(OBJ_DIR)/mmer.o: $(SRC_DIR)/mmer.cpp $(SRC_DIR)/mmer.h $(SRC_DIR)/kmer_defs.h
$(OBJ_DIR)/kmer_api.o: $(SRC_DIR)/kmer_api.cpp $(SRC_DIR)/kmer_api.h $(SRC_DIR)/kmer_defs.h
$(OBJ_DIR)/kstring.o: $(SRC_DIR)/kstring.c $(SRC_DIR)/kstring.h
$(OBJ_DIR)/kthread.o: $(SRC_DIR)/kthread.c $(SRC_DIR)/kthread.h
$(OBJ_DIR)/kbam.o: $(SRC_DIR)/kbam.c
$(OBJ_DIR)/kaddp.o: $(SRC_DIR)/kaddp.c
$(OBJ_DIR)/functions.o: $(SRC_DIR)/functions.c

DEPFILES := $(OBJS:.o=.d)
-include $(DEPFILES)

$(OBJ_DIR)/%.d: $(SRC_DIR)/%.c
	@$(CC) $(CFLAGS) -MM -MT $(@:.d=.o) -MF $@ $<

$(OBJ_DIR)/%.d: $(SRC_DIR)/%.cpp
	@$(CXX) $(CXXFLAGS) -MM -MT $(@:.d=.o) -MF $@ $<

install: $(PROG)
	cp $(PROG) $(BIN_DIR)

clean:
	rm -rf $(PROG) $(OBJ_DIR)/*.o $(OBJ_DIR)/*.d $(SRC_DIR)/*.dSYM KMERIA.tar.gz

distclean: clean
	rm -rf $(OBJ_DIR)

package: clean
	tar -zcf KMERIA.tar.gz Makefile bin/ src/ include/ lib/ scripts/ external_tools/ kmeria_env.yml
