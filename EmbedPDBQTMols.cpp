/**
 * @file EmbedPDBQTMols.cpp
 * @author Louai KASSA BAGHDOUCHE
 * @brief A C++ code to generate conformers from multiple pdbqt files, by fetching the SMILES and ID from PDBQT
 * @date 2022-03-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FragCatalog/FragFPGenerator.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <fstream>
#include <string>

using namespace std;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;
using namespace boost::filesystem;
using boost_ifstream = boost::filesystem::ifstream;
using boost_ofstream = boost::filesystem::ofstream;


int main(int argc, char* argv[])
{
  if (argc != 3)
  {
    cerr << "Usage: ./EmbedPDBQTMols [PDBQT FOLDER] [CONFORMERS SDF]" << endl;
    return 1;
  }

  // Obtain the files from the CLI argument
  const auto pdbqt_folder = argv[1];
  const auto conformers_file = argv[2];

  // lambda function to remove whitespaces
  auto stripWhiteSpaces = [](string& str) {
      str.erase(remove(str.begin(), str.end(), ' '), str.end());
  };

  // Inialize the output file of conformers
  boost_ofstream conf_file(conformers_file, ios::app);

  // Initalize constants
  const string pdbqt_extension = ".pdbqt";
  
  // Initialize variables
  string line, compound, smiles;
  int pos;
  size_t counter = 0;
  EmbedParameters params(srETKDGv3);
  params.randomSeed = 209;
  SDWriter writer(&conf_file);

  // Search for pdbqt files into the pdbqt folder
  for (const auto& entry : recursive_directory_iterator(pdbqt_folder))
  {
    if (entry.path().extension() == pdbqt_extension)
    {
      cout << "Processing molecule Number° " << counter << endl;
      cout << "Path of the molecule " << entry.path() << endl;
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
          /**
           * @brief Some compounds of the REAL library strangely contain a 'q' character in their 
           * SMILES, which causes a Segfault error raised from the RDKit API.
           * Here is a check of the goodness of the SMILES molecule, to avoid the Segfault error
           */
          if (smiles.find('q') != string::npos)
          {
            cerr << "Incorrect format of the SMILES" << endl;
            break;
          }
          else
          {
            try
            {
              const unique_ptr<ROMol> smi_ptr(SmilesToMol(smiles));
              const unique_ptr<ROMol> mol_ptr(addHs(*smi_ptr));
              auto& mol = *mol_ptr;
              mol.setProp("_Name", compound);
              const auto confIds = EmbedMultipleConfs(mol, 4, params);
              if (confIds.empty())
              {
                cerr << "Error, in parsing molecule. Conformers not generated!" << endl;
                break;
              }
              cout << confIds.size() << " Conformers of " << compound << '\t' << smiles << " are succefully generated!" << endl;
              // Writing conformers
              for (const auto confId : confIds)
              {
                writer.write(mol, confId);
              }
              counter++;
              break;
            }
            catch (const MolSanitizeException& e)
            {
              cerr << "Kekulization problem!" << endl;
              break;
            }
          }
        }
      }
    }
  }
  cout << "Process completed!" << endl;
}