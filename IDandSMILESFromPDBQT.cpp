/**
 * @file IDandSMILESFromPDBQT.cpp
 * @author Louai KASSA BAGHDOUCHE
 * @brief Fetch IDs and SMILES from a folder of PDBQT files
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

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        cerr << "Usage: ./id_smiles [PDBQT FOLDER] [REAL ID FILE] [SMILES FILE]" << endl;
        return 1;
    }
    // read CLI arguments
    const path pdbqt_folder = argv[1];
    const path realid_file = argv[2];
    const path smiles_file = argv[3];

    // ofstream instances
    boost_ofstream realidsfile(realid_file, ios::app);
    boost_ofstream smilesfile(smiles_file, ios::app);

    // lambda function to remove whitespaces
    auto stripWhiteSpaces = [](string& str) {
        str.erase(remove(str.begin(), str.end(), ' '), str.end());
    };

    // constants
    const string pdbqt_extension = ".pdbqt";

    // Variables
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
                    realidsfile << compound << '\n';
                }

                if (line.find("SMILES:") != string::npos)
                {
                    pos = line.find(':');
                    smiles = line.substr(pos + 1);
                    stripWhiteSpaces(smiles);
                    smilesfile << smiles << '\n';
                    break;
                }
            }
            counter++;
        }
    }
    cout << "Process finished, generating id and smiles for " << counter << " pdbqt files." << endl;
}
