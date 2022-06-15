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
    boost_ifstream ifs(src_file, ios::in |ios::binary | ios::ate); // ios::ate to seek to the end of file, for getting the number of bytes
    boost_ofstream ofs(output_file, ios::app | ios::binary);
    const size_t num_bytes = ifs.tellg();
    cerr << "Reading file of " << num_bytes << " bytes." << endl;
    ifs.seekg(0);
    T buffer_array;
    OT buffer_output;
    unsigned int counter = 0;
    while (!ifs.eof()) {
        ifs.read(reinterpret_cast<char*>(buffer_array.data()), sizeof(T));
        cerr << "Processing molecule N° " << counter++ << endl;
        for (int i = 0; i < USRCAT_FEATURES; i++) {
            buffer_output[i] = (int16_t)(buffer_array[i] * 1000);
        }
        cout << buffer_output;
        ofs.write(reinterpret_cast<char*>(buffer_output.data()), sizeof(OT));
    }
}

int main(int argc, char* argv[]) {
    const path input32_t = argv[1];
    const path output16_t = argv[2];
    read_and_transform<Encode_float32_t, Encode_int16_t>(input32_t, output16_t);
}


