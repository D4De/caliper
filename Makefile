CXX=g++
CFLAGS=-c -O3
LDFLAGS=-lm

CALIPER_SRC=src/caliper.cpp src/utils.cpp
THERMAL_SRC=src/thermal.cpp src/utils.cpp
CALIPER_OBJ=$(CALIPER_SRC:.cpp=.o)
THERMAL_OBJ=$(THERMAL_SRC:.cpp=.o)
CALIPER_EX=caliper
THERMAL_EX=dummythermal

all: $(CALIPER_EX) $(THERMAL_EX)

$(CALIPER_EX): $(CALIPER_OBJ)  
	$(CXX) $(CALIPER_OBJ) $(LDFLAGS) -o $@

$(THERMAL_EX): $(THERMAL_OBJ)  
	$(CXX) $(THERMAL_OBJ) $(LDFLAGS) -o $@

.cpp.o:
	$(CXX) $(CFLAGS) $< -o $@
	
clean:
	rm -rf src/*.o $(CALIPER_EX) $(THERMAL_EX) 
