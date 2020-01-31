CC = g++ -std=gnu++11

CFLAGS =  -g -pedantic -Wall -fcheck=all -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused -fopenmp -c -O3 $(DEBUG)
#CFLAGS =  -g -Wall -Wshadow -fopenmp -c -O3 $(DEBUG)
#-Werror
LFLAGS =  -g -Wall -Wshadow -fopenmp $(DEBUG)
MFLAGS =  -g -Wall -O3 $(DEBUG) 

grain3d : grain3d.o cnode.o cbody.o functions.o globals.o topchanges.o     
	$(CC) $(LFLAGS) grain3d.o cbody.o cnode.o functions.o globals.o topchanges.o  -o grain3d

grain3d.o : grain3d.cpp 
	$(CC) $(CFLAGS) grain3d.cpp

cnode.o : cnode.cpp
	$(CC) $(CFLAGS) cnode.cpp

cbody.o : cbody.cpp
	$(CC) $(CFLAGS) cbody.cpp

functions.o : functions.cpp
	$(CC) $(CFLAGS) functions.cpp

globals.o : globals.cpp
	$(CC) $(CFLAGS) globals.cpp

topchanges.o : topchanges.cpp
	$(CC) $(CFLAGS) topchanges.cpp

