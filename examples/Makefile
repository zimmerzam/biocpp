CC = g++46
CFLAGS = -msse2 -O2 -Wall -std=c++0x -static

EXECUTABLE = container ids pdb fasta_align iterate dpss zimm_align surface_area contact_map gaussian_network_model

all: $(EXECUTABLE)

$(EXECUTABLE):
	$(CC) $(CFLAGS) $@.cpp -o $@

clean:
	rm -f $(EXECUTABLE) *~
