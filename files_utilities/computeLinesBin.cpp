#include <iostream>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>

using namespace std;
using namespace boost::filesystem;
using boost_ofstream = boost::filesystem::ofstream;
using boost_ifstream = boost::filesystem::ifstream;

template<typename T>
static inline vector<T> read(const path src)
{
	boost_ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	cout << "Reading " << src << " of " << num_bytes << " bytes" << endl;
	vector<T> buf;
	buf.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(buf.data()), num_bytes);
	return buf;
}

template<typename T>
void print_arr(const T& arr)
{
	cout << '[';
	for (size_t i = 0; i < arr.size(); i++)
	{
		(i != arr.size() - 1) ? cout << arr[i] << ',' : cout << arr[i] << ']' << endl;
	}
}

int main(int argc, char* argv[])
{
	vector<array<float, 4>> biased_4prop;
	vector<float> biased_feats;
	vector<array<int16_t, 5>> biased_5prop;
	const path bin_path = argv[1];
	const path output4 = "/home/kassa-baghdouche/Documents/USR-VS/USR_FILES/final/ligand_final_4properties.f32";
	const path output5 = "/home/kassa-baghdouche/Documents/USR-VS/USR_FILES/final/ligand_CB_final_5properties.i16";
	const path output_usrcat = "/home/kassa-baghdouche/Documents/USR-VS/USR_FILES/final/ligand_CB_final_usrcat.f64";

	boost_ofstream ofs(output4, ios::binary | ios::app);
	boost_ofstream ofs5(output5, ios::binary | ios::app);
	boost_ofstream usrcat(output_usrcat, ios::binary | ios::app);

	if (bin_path.extension() == ".f32")
	{
		vector<array<float, 4>> prop = read<array<float, 4>>(bin_path);
		cout << "This file contains " << prop.size() << " lines" << endl;
		for (size_t i = 0; i < 507692; i++)
		{
			auto& buf = prop[i];
			const size_t num_bytes = sizeof(buf);
			ofs.write(reinterpret_cast<char*>(buf.data()), num_bytes);
		}
	}
	else if (bin_path.extension() == ".i16")
	{
		vector<array<int16_t, 5>> prop = read<array<int16_t, 5>>(bin_path);
		cout << "This file contains " << prop.size() << " lines" << endl;
		for (size_t i = 0; i < 507692; i++)
		{
			auto& buf = prop[i];
			const size_t num_bytes = sizeof(buf);
			ofs5.write(reinterpret_cast<char*>(buf.data()), num_bytes);
		}
	}
	else if (bin_path.extension() == ".f64")
	{
		vector<array<float, 60>> prop = read<array<float, 60>>(bin_path);
		cout << "This file contains " << prop.size() << " lines" << endl;
		for (size_t i = 0; i < 21222940; i++)
		{
			auto& buf = prop[i];
			const size_t num_bytes = sizeof(buf);
			usrcat.write(reinterpret_cast<char*>(buf.data()), num_bytes);
		}
	}
}