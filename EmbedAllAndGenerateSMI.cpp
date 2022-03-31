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
#include <GraphMol/FragCatalog/FragFPGenerator.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>

using namespace std;
using namespace std::chrono;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;
using namespace boost::filesystem;
using boost_ifstream = boost::filesystem::ifstream;
using boost_ofstream = boost::filesystem::ofstream;

// function to remove whitespaces
inline auto stripWhiteSpaces(string& str)
{
  str.erase(remove(str.begin(), str.end(), ' '), str.end());
}

// Inline function to generate conformers
inline auto EmbedConformers(ROMol& mol, EmbedParameters& params)
{
  return EmbedMultipleConfs(mol, 4, params);
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

  // generate smiles.txt and realid.txt
  int p = string(smi_file).find('.');
  path only_smiles_txt = string(smi_file).substr(0, p) + "_only_smiles.txt";
  path only_id_txt = string(smi_file).substr(0, p) + "_only_id.txt";

  // Inialize the output file of conformers
  boost_ofstream conf_file(conformers_file, ios::app);
  boost_ofstream smifile(smi_file, ios::app);
  boost_ofstream smilesfile(only_smiles_txt, ios::app);
  boost_ofstream id_file(only_id_txt, ios::app);

  // Initalize constants
  const string pdbqt_extension = ".pdbqt";
  
  // Initialize variables
  string line, compound, smiles, next_line;
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
            id_file << compound << '\n';
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
              // to optimize performance, only conformers that are generated in less than 2500 ms are kept
              future<INT_VECT> future_conformers = async(launch::deferred, EmbedConformers, ref(mol), ref(params));
              auto status = future_conformers.wait_for(milliseconds(2500));
              if (status == future_status::timeout)
              {
                cerr << "Conformers take too long to be generated" << endl;
                break;
              }
              const auto confIds = future_conformers.get();
              if (confIds.empty())
              {
                cerr << "Error, in parsing molecule. Conformers not generated!" << endl;
                break;
              }
              cout << confIds.size() << " Conformers of " << compound << '\t' << smiles << " are succefully generated!" << endl;
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