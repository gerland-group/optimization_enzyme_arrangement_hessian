EXECUTABLES:=optimise 
SYS_LIBS=
DEBUG:=0 
OPT:=2 
MAKEFLAGS=

##**Nothing below here should need to be changed**
CC=g++
CPPFLAGS=-ggdb$(DEBUG) -O$(OPT) -c -Wall -DHAVE_INLINE -W -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -fomit-frame-pointer -flto -fopenmp
LDFLAGS=-ggdb$(DEBUG) -flto -fopenmp

DEPDIR:=.deps
MAINS:=$(addsuffix .cpp,$(EXECUTABLES))
SOURCES:=$(filter-out $(MAINS),$(wildcard *.cpp))
OBJECTS:=$(SOURCES:.cpp=.o)
DEPS:=$(addprefix ${DEPDIR}/,${SOURCES:.cpp=.Po} ${MAINS:.cpp=.Po})
.SUFFIXES: .cpp .o .Po .h

all: $(EXECUTABLES) $(DEPS)

$(EXECUTABLES): %: %.o $(OBJECTS)
	${CC} ${LDFLAGS} $(OBJECTS) $@.o ${SYS_LIBS} -o $@

clean:
	rm -f *.o $(EXECUTABLES)

$(DEPDIR)/%.Po: %.cpp
	[ -d $(DEPDIR) ] || mkdir $(DEPDIR);\
		$(CC) $(CPPFLAGS) -M -MF $@ $<

include $(DEPS)

