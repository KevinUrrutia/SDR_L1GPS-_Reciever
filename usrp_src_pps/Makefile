CC=g++
CFLAGS = -lboost_program_options -lboost_iostreams -pthread -luhd -lstdc++fs -o logger
SOURCES=main.cpp
EXECUTABLE=IQ_Recorder

$(EXECUTABLE): $(SOURCES)
	$(CC) $(SOURCES) $(CFLAGS)
#	$(CC) main.o -o logger $(CFLAGS)

.PHONY : clean
clean :
	rm logger
