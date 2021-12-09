all: generateStats.exe

generateStats.exe: generateStats.cpp AMELLIE_utils.hpp AMELLIE_utils.cpp
	g++ -g -W -Wall -std=c++1y generateStats.cpp AMELLIE_utils.hpp  AMELLIE_utils.cpp -o generateStats.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -L${RATROOT}/lib -lRATEvent_Linux
	

clean :
	rm *.o generateStats.exe
