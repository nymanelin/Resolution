TARGET=resolution
INCS=-I$(shell root-config --incdir)
LIBS=$(shell root-config --libs)
all: $(TARGET)

.PHONY: clean
clean:
	rm -f $(TARGET)

$(TARGET):
	g++  $(TARGET).cpp -o $(TARGET) -O3 $(INCS) $(LIBS) -std=c++11 

