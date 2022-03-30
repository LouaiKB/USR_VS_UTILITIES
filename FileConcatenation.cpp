/**
 * @file FileConcatenation.cpp
 * @author Louai KASSA BAGHDOUCHE
 * @brief Concatenation of large f64 binary  files
 * @version 0.1
 * @date 2022-03-30
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <iostream>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <chrono>

using namespace std;
using namespace std::chrono;
using namespace boost::filesystem;
using namespace boost::posix_time;
using boost_ofstream = boost::filesystem::ofstream;
using boost_ifstream = boost::filesystem::ifstream;

auto local_time()
{
    return to_simple_string(microsec_clock::local_time()) + " ";
}

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cerr << "Usage: ./file_concatenation [first_file] [second_file]";
        return 1;
    }

    cout << local_time() << "Process started!" << endl;
    const auto started = system_clock::now();

    const auto first_file = argv[1];
    const auto second_file = argv[2];

    boost_ofstream of_original(first_file, ios::binary | ios::app);
    boost_ifstream if_file(second_file, ios::binary);

    of_original.seekp(0, ios::end);
    of_original << if_file.rdbuf();

    cout << local_time() << "Process finished!" << endl;
    const auto completed = system_clock::now();
    const auto runtime = (completed - started).count() * 1e-9;

    cout << "Concatenation of the the two files completed in " << runtime << " seconds" << endl;
}