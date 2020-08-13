CXXFLAGS := $(CXXFLAGS) -std=c++11 -pedantic -Wall -Wextra
COMMON_DEFINES := -D POPS_TEST
INCLUDES := -Iinclude/

all: test_date test_raster test_simulation test_treatments test_spread_rate test_statistics test_scheduling test_deterministic test_model test_quarantine test_random

test_date: tests/test_date.cpp include/pops/date.hpp
	g++ $(CXXFLAGS) $(INCLUDES) $(COMMON_DEFINES) tests/test_date.cpp -o test_date

test_raster: tests/test_raster.cpp include/pops/*.hpp
	g++ $(CXXFLAGS) $(INCLUDES) $(COMMON_DEFINES) tests/test_raster.cpp -o test_raster

test_simulation: tests/test_simulation.cpp include/pops/*.hpp
	g++ $(CXXFLAGS) $(INCLUDES) $(COMMON_DEFINES) tests/test_simulation.cpp -o test_simulation

test_treatments: tests/test_treatments.cpp include/pops/*.hpp
	g++ $(CXXFLAGS) $(INCLUDES) $(COMMON_DEFINES) tests/test_treatments.cpp -o test_treatments

test_spread_rate: tests/test_spread_rate.cpp include/pops/*.hpp
	g++ $(CXXFLAGS) $(INCLUDES) $(COMMON_DEFINES) tests/test_spread_rate.cpp -o test_spread_rate

test_statistics: tests/test_statistics.cpp include/pops/*.hpp
	g++ $(CXXFLAGS) $(INCLUDES) $(COMMON_DEFINES) tests/test_statistics.cpp -o test_statistics

test_scheduling: tests/test_scheduling.cpp include/pops/*.hpp
	g++ $(CXXFLAGS) $(INCLUDES) $(COMMON_DEFINES) tests/test_scheduling.cpp -o test_scheduling

test_deterministic: tests/test_deterministic.cpp include/pops/*.hpp
	g++ $(CXXFLAGS) $(INCLUDES) $(COMMON_DEFINES) tests/test_deterministic.cpp -o test_deterministic

test_model: tests/test_model.cpp include/pops/*.hpp
	g++ $(CXXFLAGS) $(INCLUDES) $(COMMON_DEFINES) tests/test_model.cpp -o test_model

test_quarantine: tests/test_quarantine.cpp include/pops/*.hpp
	g++ $(CXXFLAGS) $(INCLUDES) $(COMMON_DEFINES) tests/test_quarantine.cpp -o test_quarantine

test_random: tests/test_random.cpp include/pops/*.hpp
	g++ $(CXXFLAGS) $(INCLUDES) $(COMMON_DEFINES) tests/test_random.cpp -o test_random

test:
	./test_date
	./test_raster
	./test_simulation
	./test_treatments
	./test_spread_rate
	./test_statistics
	./test_scheduling
	./test_deterministic
	./test_model
	./test_quarantine
	./test_random

doc:
	doxygen

clean:
	rm -f test_date
	rm -f test_raster
	rm -f test_simulation
	rm -f test_treatments
	rm -f test_spread_rate
	rm -f test_statistics
	rm -f test_scheduling
	rm -f test_deterministic
	rm -f test_model
	rm -f test_quarantine
	rm -f test_random
