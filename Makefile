build:
	g++ -O2 -c -std=c++17 *.cpp */*.cpp -Wno-narrowing
	g++ -o ml *.o -lpng
	./ml


clean: 
	rm -f *.o
