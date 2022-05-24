#include <iostream>

using namespace std;

int main()
{
	size_t num_bytes;
	for (string line; getline(cin, line);)
	{
		num_bytes = cin.tellg();
		cout.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
	}
}