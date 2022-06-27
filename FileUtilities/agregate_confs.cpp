#include <iostream>

int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cerr << "./agregate_confs [NUM_MOLECULES]" << std::endl;
		return 1;
	}
	char* endPtr;
	const uint64_t num_mols = strtoul(argv[1], &endPtr, 10);
	uint64_t agregated_confs = 0;
	for (uint64_t i = 0; i < num_mols; i++) {
		agregated_confs += 4;
		uint64_t* address = &agregated_confs;
		const unsigned long size_buffer = sizeof(agregated_confs);
		std::cout.write(reinterpret_cast<char*>(address), size_buffer);
	}
}
