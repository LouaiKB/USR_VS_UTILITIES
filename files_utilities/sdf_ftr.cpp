#include <iostream>

int main(int argc, char* argv[])
{
	size_t num_bytes = 0;
	for (std::string line; std::getline(std::cin, line);)
	{
		if (line == "$$$$")
		{
			num_bytes = std::cin.tellg();
			std::cout.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
		}
	}
}