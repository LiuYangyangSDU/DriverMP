CXX=g++
CXXFLAGS=-Ofast -fopenmp

OBJDIR=.obj
SRCDIR=src
INCDIR=include

SRCS=$(wildcard $(SRCDIR)/*.cpp)

HEADERS=$(wildcard $(INCDIR)/*.h)

OBJS=$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

TARGET=DriverMP

all: $(TARGET)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c -o $@ $<

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -rf $(OBJDIR) $(TARGET)