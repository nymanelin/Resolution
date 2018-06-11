TARGET=resolution
INCS=-I$(shell root-config --incdir)
LIBS=$(shell root-config --libs)
all: $(TARGET)

.PHONY: clean
clean:
	rm -f $(TARGET)

$(TARGET):
	g++ -O3 $(INCS) $(LIBS) -std=c++11 -o $(TARGET) $(TARGET).cpp

