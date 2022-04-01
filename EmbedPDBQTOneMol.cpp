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
#include <future>
#include <fstream>
#include <string>
#include <chrono>
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

class Thread_job
{
public:
  explicit Thread_job(ROMol, EmbedParameters&);
  ~Thread_job();
  inline void run(SDWriter&);
  inline INT_VECT embedMolecules();
  const inline ROMol& getMolecule() const;
  inline bool conformersGenerated();
  const inline int getNumConformers() const;

private:
  const unsigned int m_NumConfs = 4;
  ROMol m_Mol;
  INT_VECT m_ConfIds;
  const EmbedParameters m_Params;
  condition_variable _cv;
  mutex _mu;
};

Thread_job::Thread_job(ROMol mol, EmbedParameters& params) : m_Mol(mol), m_Params(params) {}
Thread_job::~Thread_job() {}

inline INT_VECT Thread_job::embedMolecules()
{
  return EmbedMultipleConfs(m_Mol, m_NumConfs, m_Params);
}

inline void Thread_job::run(SDWriter& sdf_writer)
{

  thread thread_job([this]() {
    m_ConfIds = embedMolecules();
    _cv.notify_one();
  });
  thread_job.detach();
  {
    unique_lock<mutex> l(_mu);
    if (_cv.wait_for(l, 1.5s) == cv_status::timeout)
      throw runtime_error("Timeout");
  }
  if (m_ConfIds.size() != 4)
  {
    m_ConfIds.clear();
    return;
  }
  // writing in the sdf file
  for (const auto confId : m_ConfIds)
  {
    sdf_writer.write(m_Mol, confId);
  }
}

const inline ROMol& Thread_job::getMolecule() const
{
  return m_Mol;
}

inline bool Thread_job::conformersGenerated()
{
  return m_ConfIds.empty();
}

const inline int Thread_job::getNumConformers() const
{
  return m_ConfIds.size();
}

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
  if (argc != 3)
  {
    cerr << "Usage: ./EmbedPDBQTMols [PDBQT FOLDER] [CONFORMERS SDF]" << endl;
    return 1;
  }

  // Obtain the files from the CLI argument
  const auto pdbqt_file = argv[1];
  const auto conformers_file = argv[2];

  // Inialize the output file of conformers
  boost_ofstream conf_file(conformers_file, ios::app);

  // Initalize constants
  const string pdbqt_extension = ".pdbqt";
  
  // Initialize variables
  string line, compound, smiles, next_line;
  int pos;
  EmbedParameters params(srETKDGv3);
  params.randomSeed = 209;
  SDWriter writer(&conf_file);
  INT_VECT confIds;
  future<INT_VECT> future_confs;

  boost_ifstream ifs(pdbqt_file);
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
       * @brief Some compounds of the REAL library strangely contain a 'q' and 'r' characters in their 
       * SMILES, and other compounds have a splitted SMILES in two lines. These problems cause a Segfault error raised from the RDKit API.
       * Here is a check of the goodness of the SMILES molecule, to avoid the Segfault error
       */
      if (getline(ifs, next_line))
      {
        if (next_line.find("REMARK") != string::npos)
        {
          cout << "No problem with the SMILES." << endl;
        }
        else
        {
          cerr << "The SMILES is splited, trying to fix the SMILES.." << endl;
          cout << next_line << endl;
          smiles = smiles + next_line;
          cout << "The SMILES IS : " << smiles << endl;
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
          {
            Thread_job thread(mol, params);
            try
            {
              thread.run(writer);
              if (thread.conformersGenerated())
              {
                cout << "Error conformers not generated" << endl;
                break;
              } 
              else
              {
                cout << thread.getNumConformers() << " Conformers of " << compound << '\t' << smiles << " are succefully generated!" << endl;
                break;
              }
            }
            catch(const std::runtime_error& e)
            {
              std::cerr << e.what() << '\n';
              break;
            }
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
  cout << "Process completed!" << endl;
}