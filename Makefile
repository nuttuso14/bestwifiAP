CC = g++
CFLAGS = -g -c
STANDARD = -std=c++11


all: bestWifiAp 

bestWifiAp: bestWifiAp.o
	$(CC) -o $@ $? 

bestWifiAp.o:	bestWifiAp.cpp
	$(CC) $(STANDARD) $(CFLAGS) -o $@ bestWifiAp.cpp 

clean:
	rm *.o *.txt bestWifiAp 

