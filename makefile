EXEDIR := run
OBJDIR := bin
SRCDIR := src
INCDIR := inc
MAKEDIR := bin
LIBFILE := $(OBJDIR)/libStatObj.a

CXX := $(shell root-config --cxx)
EXTRA_WARNINGS := -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k -Winit-self -Winvalid-pch -Wlong-long -Wmissing-format-attribute -Wmissing-include-dirs -Wmissing-noreturn -Wpacked -Wpointer-arith -Wredundant-decls -Wstack-protector -Wswitch-default -Wswitch-enum -Wundef -Wunused -Wvariadic-macros -Wwrite-strings -Wctor-dtor-privacy -Wnon-virtual-dtor -Wsign-promo -Wsign-compare #-Wunsafe-loop-optimizations -Wfloat-equal -Wsign-conversion -Wunreachable-code
CXXFLAGS := -isystem $(shell root-config --incdir) -Wall -Wextra -pedantic -Werror -Wshadow -Woverloaded-virtual -Wold-style-cast $(EXTRA_WARNINGS) $(shell root-config --cflags) -O2 -I $(INCDIR)
LD := $(shell root-config --ld)
LDFLAGS := $(shell root-config --ldflags) -lGenVector
LDLIBS := $(shell root-config --libs) -lMinuit -lRooStats -lTreePlayer 

EXECUTABLES := $(addprefix $(EXEDIR)/, $(addsuffix .exe, $(notdir $(basename $(wildcard $(SRCDIR)/*.cxx))))) 
OBJECTS := $(addprefix $(OBJDIR)/, $(addsuffix .o, $(notdir $(basename $(wildcard $(SRCDIR)/*.cpp)))))

FIND_DEPS = $(CXX) $(CXXFLAGS) -MM -MG -MF $@ $<
EXPAND_DEPS = perl -pi -e 's|$*.o|$(OBJDIR)/$*.o $(MAKEDIR)/$*.d|g' $@
GET_DEPS = $(FIND_DEPS) && $(EXPAND_DEPS)
COMPILE = $(CXX) $(CXXFLAGS) -o $@ -c $<
LINK = $(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

vpath %.cpp $(SRCDIR)
vpath %.cxx $(SRCDIR)
vpath %.hpp $(INCDIR)
vpath %.o $(OBJDIR)
vpath %.exe $(EXEDIR)
vpath %.d $(MAKEDIR)

all: $(EXECUTABLES)

-include $(addsuffix .d,$(addprefix $(MAKEDIR)/,$(notdir $(basename $(wildcard $(SRCDIR)/*.cpp)))))
-include $(addsuffix .d,$(addprefix $(MAKEDIR)/,$(notdir $(basename $(wildcard $(SRCDIR)/*.cxx)))))


$(LIBFILE): $(OBJECTS)

$(MAKEDIR)/%.d: $(SRCDIR)/%.cpp
	$(GET_DEPS)

$(MAKEDIR)/%.d: $(SRCDIR)/%.cxx
	$(GET_DEPS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMPILE)

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	$(COMPILE)

$(OBJDIR)/%.a:
	ar rcsv $@ $^

$(EXEDIR)/generate_tree_classes.exe: $(OBJDIR)/generate_tree_classes.o
	$(LINK)

$(EXEDIR)/%.exe: $(OBJDIR)/%.o $(LIBFILE)
	$(LINK)

# Auto-generated code
.SECONDARY: dummy_nano_tree.all dummy_corrections_tree.all dummy_pico_tree.all dummy_baby_tree.all
.PRECIOUS: generate_tree_classes.o 

$(SRCDIR)/nano_tree.cpp $(INCDIR)/nano_tree.hpp: dummy_nano_tree.all
dummy_nano_tree.all: $(EXEDIR)/generate_tree_classes.exe 
	./$< 

$(SRCDIR)/corrections_tree.cpp $(INCDIR)/corrections_tree.hpp: dummy_corrections_tree.all
dummy_corrections_tree.all: $(EXEDIR)/generate_tree_classes.exe 
	./$< 

$(SRCDIR)/pico_tree.cpp $(INCDIR)/pico_tree.hpp: dummy_pico_tree.all
dummy_pico_tree.all: $(EXEDIR)/generate_tree_classes.exe 
	./$< 

$(SRCDIR)/baby_tree.cpp $(INCDIR)/baby_tree.hpp: dummy_baby_tree.all
dummy_baby_tree.all: $(EXEDIR)/generate_tree_classes.exe 
	./$< 

.DELETE_ON_ERROR:
