EXEC = genePhoto

all : $(EXEC)

genePhoto.o : genePhoto.cpp
	g++ -c genePhoto.cpp

$(EXEC) : genePhoto.o
	g++ genePhoto.o -o $@

clean :
	rm -f *.o

mrproper: clean
	rm $(EXEC)



