/**
 * @file GenerateSmiFiles.cpp
 * @author Louai KASSA BAGHDOUCHE
 * @brief Generate SMI files from a folder of PDBQT files
 * @version 0.1
 * @date 2022-03-30
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>

using namespace std;
using namespace boost::filesystem;
using boost_ofstream = boost::filesystem::ofstream;
using boost_ifstream = boost::filesystem::ifstream;

// function to remove whitespaces
inline auto stripWhiteSpaces(string& str)
{
    str.erase(remove(str.begin(), str.end(), ' '), str.end());
}


int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cerr << "Usage: ./id_smiles [PDBQT FOLDER] [OUTPUT SMI FILE]" << endl;
        return 0;
    }

    // read CLI arguments
    const path pdbqt_folder = argv[1];
    const path smi_file = argv[2];

    // ofstream instances
    boost_ofstream smilesfile(smi_file, ios::app);

    // constants
    const string pdbqt_extension = ".pdbqt";
    string line, compound, smiles;
    int pos;
    size_t counter = 0;

    for (const auto& entry : recursive_directory_iterator(pdbqt_folder))
    {
        if (entry.path().extension() == pdbqt_extension)
        {
            cout << "Porcessing molecule Number° " << counter << endl;
            boost_ifstream ifs(entry.path());
            while (getline(ifs, line))
            {
                if (line.find("Compound:") != string::npos)
                {
                    pos = line.find(':');
                    compound = line.substr(pos + 1);
                    stripWhiteSpaces(compound);
                }

                if (line.find("SMILES:") != string::npos)
                {
                    pos = line.find(':');
                    smiles = line.substr(pos + 1);
                    stripWhiteSpaces(smiles);
                    smilesfile << compound << '\t' << smiles << '\n';
                    break;
                }
            }
            counter++;
        }
    }
    cout << "Process finished, generating id and smiles for " << counter << " pdbqt files." << endl;
}