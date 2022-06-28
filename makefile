EXEDIR := run
OBJDIR := bin
ZOBJDIR := zlib-1.2.12
SRCDIR := src
INCDIR := inc
ZLIBDIR := include
MAKEDIR := bin
LIBFILE := $(OBJDIR)/libStatObj.a

CXX := $(shell root-config --cxx)
EXTRA_WARNINGS := -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k -Winit-self -Winvalid-pch -Wlong-long -Wmissing-format-attribute -Wmissing-include-dirs -Wmissing-noreturn -Wpacked -Wpointer-arith -Wredundant-decls -Wstack-protector -Wundef -Wvariadic-macros -Wwrite-strings -Wctor-dtor-privacy -Wnon-virtual-dtor -Wsign-promo -Wsign-compare -Wno-unused-variable # -Wunsafe-loop-optimizations -Wfloat-equal -Wsign-conversion -Wunreachable-code -Wswitch-enum -Wunused
CXXFLAGS := -isystem $(shell root-config --incdir) -Wall -Wextra -pedantic -Werror -Wold-style-cast $(EXTRA_WARNINGS) -pthread -std=c++17 -m64 -I/cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw-patch/CMSSW_10_2_11_patch1/external/slc6_amd64_gcc700/bin/../../../../../../../slc6_amd64_gcc700/lcg/root/6.12.07-gnimlf5/include -O2 -I $(INCDIR) -I $(ZLIBDIR) # -Wshadow -Woverloaded-virtual
LD := $(shell root-config --ld)
LDFLAGS := $(shell root-config --ldflags) -lGenVector
LDLIBS := -L/homes/psiddire/Central/nano2pico/lib -lz $(shell root-config --libs) -lMinuit -lRooStats -lTreePlayer
# root-config --cflags

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
vpath %.h $(ZLIBDIR)
vpath %.o $(OBJDIR)
vpath %.o $(ZOBJDIR)
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
