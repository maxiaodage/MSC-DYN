
output: main.o
		g++ main.o -o output -fopenmp

main.o:main.C
		g++ -c main.C  -fopenmp -I ./
clean:
		rm *.o output
