all: generateStats.exe getRegions.exe getTrackingInfo.exe

generateStats.exe: generateStats.cpp AMELLIE_utils.hpp AMELLIE_utils.cpp
	g++ -g -W -Wall -std=c++1y generateStats.cpp AMELLIE_utils.hpp  AMELLIE_utils.cpp -o generateStats.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
	
getRegions.exe: getRegions.cpp AMELLIE_utils.hpp AMELLIE_utils.cpp
	g++ -g -W -Wall -std=c++1y getRegions.cpp AMELLIE_utils.hpp  AMELLIE_utils.cpp -o getRegions.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux

getTrackingInfo.exe: getTrackingInfo.cpp
	g++ -g -std=c++1y getTrackingInfo.cpp -o getTrackingInfo.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
	
clean :
	rm *.exe
