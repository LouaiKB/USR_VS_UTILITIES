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
}


