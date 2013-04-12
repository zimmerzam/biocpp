CC = g++
CFLAGS = -msse2 -O2 -Wall -std=c++0x -static

EXECUTABLE = container ids pdb fasta_align iterate dpss zimm_align surface_area, fuzzy_contact_map

all: $(EXECUTABLE)

$(EXECUTABLE):
	$(CC) $(CFLAGS) $@.cpp -o $@

clean:
	rm -f $(EXECUTABLE) *~
