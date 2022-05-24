#include <iostream>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <vector>

using namespace std;
using namespace boost::filesystem;

template <typename T>
inline vector<T> read(const path src)
{
	boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	cout << "Reading " << src << " of " << num_bytes << " bytes!" << endl;
	vector<T> buf;
	buf.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(buf.data()), num_bytes);
	return buf;
}

int main(int argc, char* argv[])
{
	uint64_t agregated_confs = 0;
	const size_t num_mols = 14646562; // for examples
	for (size_t i = 0; i < num_mols; i++)
	{
		agregated_confs += 4;
		cout.write(reinterpret_cast<char*>(&agregated_confs), sizeof(agregated_confs));
	}
	// to read the molecules
	const auto mconfs = read<size_t>("mconfs.u64");
	for (size_t i = 0; i < 100; i++)
	{
		cout << mconfs[i] << endl;
	}
}