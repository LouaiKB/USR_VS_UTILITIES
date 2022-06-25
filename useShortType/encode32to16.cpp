#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#define USRCAT_FEATURES 60

namespace DEBUG {};

using namespace std;
using namespace DEBUG;
using namespace boost::filesystem;
using boost_ifstream = boost::filesystem::ifstream;
using boost_ofstream = boost::filesystem::ofstream;
using Encode_int16_t = array<int16_t, USRCAT_FEATURES>;
using Encode_float32_t = array<_Float32, USRCAT_FEATURES>;

/**
 * @brief this namespace is created for operator overloading to print in stdout an array format
 * for debugging purposes 
 * 
 */

namespace DEBUG {
template<typename T, size_t size>
ostream& operator<<(ostream& stream, const array<T, size>& arr) {
    stream << '[';
    for (unsigned i = 0; i < arr.size(); i++) (i == arr.size() - 1) ? stream << arr[i] : stream << arr[i] << " ,";
    stream << ']' << endl;
    return stream;
}
template<typename T, size_t size>
ostream& operator<<(ostream& stream, const vector<array<T, size>>& vect) {
    stream << '[';
    for (size_t i = 0; i < vect.size(); i++) {
        const array<T, size> arr = vect[i];
        stream << '[';
        for (unsigned j = 0; j < arr.size(); j++) (j == arr.size() - 1) ? stream << arr[j] : stream << arr[j] << " ,";
        (i == vect.size() - 1) ? stream << ']' : stream << ']' << ", ";
    }
    stream << ']' << endl;
    return stream;
}
}

template<typename T>
static inline vector<T> get_buffer(const path& src) {
	boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	cerr << "Reading a file of " << num_bytes << " bytes." << endl;
	ifs.seekg(0);
	vector<T> buffer;
	buffer.resize(num_bytes / sizeof(T));
	ifs.read(reinterpret_cast<char*>(buffer.data()), num_bytes);
	return buffer;
}

/**
 * @brief Read a USRCAT features encoded with the float, and write the features in short type (2bytes)
 * 
 * @tparam T 
 * @tparam OT 
 * @param src_file 
 * @param output_file 
 */
template<typename T, typename OT>
static void read_and_transform(const path& src_file, const path& output_file) {
    boost_ifstream ifs(src_file, ios::in | ios::binary | ios::ate); // ios::ate to seek to the end of file, for getting the number of bytes
    boost_ofstream ofs(output_file, ios::binary);
    const size_t num_bytes = ifs.tellg();
    cerr << "Reading file of " << num_bytes << " bytes." << endl;
    ifs.seekg(0, ios::beg);
    vector<T> buffer_array;
    OT buffer_output;
    unsigned int counter = 0;
    // while (!ifs.eof()) {
    cout << "Start reading" << endl;
    ifs.read(reinterpret_cast<char*>(buffer_array.data()), num_bytes);
    cout << "finish reading" << endl;
    // cerr << "Processing molecule N° " << counter++ << endl;
    for (const auto& buf : buffer_array) {
        for (int i = 0; i < USRCAT_FEATURES; i++) {
            buffer_output[i] = (int16_t)(buf[i] * 1000);
        }
        cout << buffer_output;
        ofs.write(reinterpret_cast<char*>(buffer_output.data()), sizeof(OT));
    }
    // cout << buffer_array;
    // cout << "--------------------------------\n";
    cout << buffer_output;
    // }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "./encode32to16 [INPUT] [OUTPUT]" << endl;
        return 1;
    }
    const path input32_t = argv[1];
    const path output16_t = argv[2];
    const auto features32_t = get_buffer<Encode_float32_t>(input32_t);
    boost_ofstream ofs(output16_t, ios::app | ios::binary);
    Encode_int16_t features16_t;
    for (const auto& buf : features32_t) {
        for (int i = 0; i < USRCAT_FEATURES; i++) {
            features16_t[i] = (int16_t)(buf[i] * 1000);
        }
        cout << features16_t;
        ofs.write(reinterpret_cast<char*>(features16_t.data()), sizeof(Encode_int16_t));
    }
    // read_and_transform<Encode_float32_t, Encode_int16_t>(input32_t, output16_t);
}


