all: generateStats.exe

generateStats.exe: AMELLIE_utils.o generateStats.o
	g++ -g AMELLIE_utils.o generateStats.o -o generateStats.exe

generateStats.o: generateStats.cpp AMELLIE_utils.hpp
	g++ -g -W -Wall -std=c++11 -c generateStats.cpp
	
AMELLIE_utils.o: AMELLIE_utils.cpp AMELLIE_utils.hpp
	g++ -g -W -Wall -std=c++11 -c AMELLIE_utils.cpp

clean :
	rm *.o generateStats.exe
