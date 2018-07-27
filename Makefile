all: test_simulation test_date

test_date: test_date.cpp date.hpp
	g++ -std=c++11 -pedantic -Wall -Wextra -D POPSS_TEST test_date.cpp -o test_date

test_simulation: test_simulation.cpp *.hpp
	g++ -std=c++11 -pedantic -Wall -Wextra -D POPSS_TEST test_simulation.cpp -o test_simulation

clean:
	rm -f test_date
	rm -f test_simulation
