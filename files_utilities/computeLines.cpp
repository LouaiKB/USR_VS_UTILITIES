#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	long num_lines = 0;
	string line;
	for (line; getline(cin, line);)
	{
		num_lines++;
		cout << num_lines << " " << line <<  endl;
	}
	cout << "Number of generated compounds are: " << num_lines << endl;
}