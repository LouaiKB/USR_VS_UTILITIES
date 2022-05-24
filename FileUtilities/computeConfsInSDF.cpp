#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	long counter = 0;
	for (string line; getline(cin, line);)
	{
		if (line == "$$$$")
		{
			counter++;
			cout << counter << endl;
		}
	}
	cout << "SDF contains " << counter << " conformer" << endl;
}