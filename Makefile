fitz: main.o fitzgibbon.o
	g++ main.o fitzgibbon.o -o fitz

main.o: main.cpp
	g++ -c main.cpp

fitzgibbon.o: fitzgibbon.cpp fitzgibbon.h
	g++ -c fitzgibbon.cpp 

clean: 
	rm *.o 

