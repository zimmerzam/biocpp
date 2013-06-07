CC = g++46
C_INC_PATH = -I /usr/include/
CFLAGS = -msse2 -O2 -Wall -std=c++0x -static

EXECUTABLE = container ids pdb fasta_align iterate dpss zimm_align surface_area contact_map gaussian_network_model anisotropic_network_model mutations alpha_variance

all: $(EXECUTABLE)

$(EXECUTABLE):
	$(CC) $(C_INC_PATH) $(CFLAGS) $@.cpp -o $@

clean:
	rm -f $(EXECUTABLE) *~
