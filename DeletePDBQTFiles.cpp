#include <iostream>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>

using namespace std;
using namespace boost::filesystem;

int main(int argc, char* argv[])
{
    const path pdbqt_folder = argv[1];
    vector<path> pdbqt_files(2818);
    const string pdbqt_extension = ".pdbqt";

    for (const auto& entry : recursive_directory_iterator(pdbqt_folder))
    {
        if (entry.path().extension() == pdbqt_extension)
        {
            pdbqt_files.push_back(entry.path());
        }
    }
    for (auto& p : pdbqt_files)
    {
        cout << "removing " << p << endl;
        remove(p);
    }    
}