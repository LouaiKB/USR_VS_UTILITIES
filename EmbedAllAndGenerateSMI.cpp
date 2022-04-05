/**
 * @file EmbedPDBQTMols.cpp
 * @author Louai KASSA BAGHDOUCHE
 * @brief A C++ code to generate conformers from multiple pdbqt files, by fetching the SMILES and ID from PDBQT and generate SMI files
 * @date 2022-03-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <future>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>

using namespace std;
using namespace std::chrono;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;
using namespace RDKit::Descriptors;
using namespace boost::filesystem;
using boost_ifstream = boost::filesystem::ifstream;
using boost_ofstream = boost::filesystem::ofstream;

// function to remove whitespaces
inline auto stripWhiteSpaces(string& str)
{
  str.erase(remove(str.begin(), str.end(), ' '), str.end());
}

int main(int argc, char* argv[])
{
  if (argc != 4)
  {
    cerr << "Usage: ./EmbedAllAndGenerateSMIL [PDBQT FOLDER] [CONFORMERS SDF] [OUTPUT SMI]" << endl;
    return 1;
  }

  // Obtain the files from the CLI argument
  const auto pdbqt_folder = argv[1];
  const auto conformers_file = argv[2];
  const auto smi_file = argv[3];

  // generate smiles.txt and realid.txt and properties
  int p = string(smi_file).find('.');
  path only_smiles_txt = string(smi_file).substr(0, p) + "_only_smiles.txt";
  path only_id_txt = string(smi_file).substr(0, p) + "_only_id.txt";
  p = string(conformers_file).find('.');
  path prop4_txt = string(conformers_file).substr(0, p) + "_4properties.txt";
  path prop5_txt = string(conformers_file).substr(0, p) + "_5properties.txt";

  // Inialize the output files
  boost_ofstream conf_file(conformers_file, ios::app);
  boost_ofstream smifile(smi_file, ios::app);
  boost_ofstream smilesfile(only_smiles_txt, ios::app);
  boost_ofstream id_file(only_id_txt, ios::app);
  boost_ofstream realprop4(prop4_txt, ios::app);
  boost_ofstream realprop5(prop5_txt, ios::app);

  // Initalize constants
  const string pdbqt_extension = ".pdbqt";
  
  // Initialize variables
  string line, compound, smiles, next_line;
  int pos;
  size_t counter = 0;
  EmbedParameters params(srETKDGv3);
  params.randomSeed = 209;
  params.numThreads = thread::hardware_concurrency() - 10;
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
          /**
           * @brief Some compounds of the REAL library strangely contain a 'q' and 'r' characters in their 
           * SMILES, and other compounds have a splitted SMILES in two lines. These problems cause a Segfault error raised from the RDKit API.
           * Here is a check of the goodness of the SMILES molecule, to avoid the Segfault error
           */
          pos = line.find(':');
          smiles = line.substr(pos + 1);
          stripWhiteSpaces(smiles);
          if (getline(ifs, next_line))
          {
            if (next_line.find("REMARK") != string::npos)
            {
              cout << "No problem with the SMILES!" << endl;
            }
            else
            {
              cerr << "The SMILES is splited, trying to fix the SMILES" << endl;
              smiles = smiles + next_line;
              cout << "The fixed SMILES: " << smiles << endl;
            }
          }
          if (smiles.find('q') != string::npos || smiles.find('r') != string::npos)
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
              // generate 4 chemical properties for the molecules
              realprop4 << calcExactMW(mol) << '\t'
                        << calcClogP(mol) << '\t'
                        << calcTPSA(mol) << '\t'
                        << calcLabuteASA(mol) << '\n';
              // generate 5 chemical properties for the molecules
              realprop5 << mol.getNumHeavyAtoms() << '\t'
                        << calcNumHBD(mol) << '\t'
                        << calcNumHBA(mol) << '\t'
                        << calcNumRotatableBonds(mol) << '\t'
                        << calcNumRings(mol) << '\n';
              // generate conformers
              const auto confIds = EmbedMultipleConfs(mol, 4, params);
              if (confIds.empty())
              {
                cerr << "Error, in parsing molecule. Conformers not generated!" << endl;
                break;
              }
              // Check if the molecule has 4 conformers
              if (confIds.size() == 4)
              {
                cout << confIds.size() << " Conformers of " << compound << '\t' << smiles << " are succefully generated!" << endl;
                id_file << compound << '\n';
                smilesfile << smiles << '\n';
                smifile << compound << '\t' << smiles << '\n';
                // Writing conformers in the output SDF
                for (const auto confId : confIds)
                {
                  writer.write(mol, confId);
                }
                counter++;
                break;
              }
              else
              {
                cout << "Molecule doesn't have 4 conformers!" << endl;
                break;
              }
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
  cout << "Process completed! for " << counter << " compounds" << endl;
}