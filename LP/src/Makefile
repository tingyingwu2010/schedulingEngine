all : plPhoto

plPhoto.o : plPhoto.cpp
	g++ -c plPhoto.cpp

plPhoto : plPhoto.o
	g++ plPhoto.o -o plPhoto

clean :
	rm -f *.o plPhoto
