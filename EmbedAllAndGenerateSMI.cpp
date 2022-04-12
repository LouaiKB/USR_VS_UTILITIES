/**
 * @file EmbedPDBQTMols.cpp
 * @author Louai KASSA BAGHDOUCHE
 * @brief A C++ code to generate conformers and their USRCAT Features from multiple pdbqt files, by fetching the SMILES and ID 
 * from PDBQT and generate SMI files
 * @date 2022-03-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <string>
#include <thread>
#include <array>
#include <vector>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>

using namespace std;
using namespace std::chrono;
using namespace RDGeom;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;
using namespace RDKit::Descriptors;
using namespace boost::filesystem;
using boost_ifstream = boost::filesystem::ifstream;
using boost_ofstream = boost::filesystem::ofstream;

static const size_t s_Num_references = 4;
static const size_t s_Num_subsets = 5;
static const array<string, s_Num_subsets> s_SubsetSMARTS
{{
  "[!#1]", // heavy
  "[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
  "[a]", // aromatic
  "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
  "[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
}};

static array<unique_ptr<ROMol>, s_Num_subsets> SubsetMols;
static array<vector<int>, s_Num_subsets> subsets;
static array<Point3D, s_Num_references> references;
static array<vector<float>, s_Num_references> dista;

// Compute the Euclidien distance
template<typename T>
static inline float dist2(const T& p0, const T& p1)
{
	const auto d0 = p0[0] - p1[0];
	const auto d1 = p0[1] - p1[1];
	const auto d2 = p0[2] - p1[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

static array<float, 60> usrcat_features(ROMol& mol, int index)
{
  array<float, 60> features;
  // Wrap SMARTS strings to ROMol objects.
  for (size_t k = 0; k < s_Num_subsets; ++k)
  {
    SubsetMols[k].reset(reinterpret_cast<ROMol*>(SmartsToMol(s_SubsetSMARTS[k])));
  }
  const auto num_points = mol.getNumHeavyAtoms();
  const auto& conformer = mol.getConformer(index);
  for (size_t k = 0; k < s_Num_subsets; ++k)
  {
    vector<vector<pair<int, int>>> matchVect;
    SubstructMatch(mol, *SubsetMols[k], matchVect);
    const auto num_matches = matchVect.size();
    auto& subset = subsets[k];
    subset.resize(num_matches);
    for (size_t j = 0; j < num_matches; ++j)
    {
      subset[j] = matchVect[j].front().second;
    }
  }
  const auto& subset0 = subsets.front();
  // assert(subset0.size() == num_points);

  for (auto& ref : references)
  {
    ref.x = ref.y = ref.z = 0;
  }
  auto& ctd = references[0];
  auto& cst = references[1];
  auto& fct = references[2];
  auto& ftf = references[3];
  for (const auto sub : subset0)
  {
    const auto& a = conformer.getAtomPos(sub);
    ctd += a;
  }
  ctd /= num_points;
  float cst_dist = numeric_limits<float>::max();
  float fct_dist = numeric_limits<float>::lowest();
  float ftf_dist = numeric_limits<float>::lowest();
  for (const auto sub : subset0)
  {
    const auto& a = conformer.getAtomPos(sub);
    const auto this_dist = dist2(a, ctd);
    if (this_dist < cst_dist)
    {
      cst = a;
      cst_dist = this_dist;
    }
    if (this_dist > fct_dist)
    {
      fct = a;
      fct_dist = this_dist;
    }
  }
  for (const auto sub : subset0)
  {
    const auto& a = conformer.getAtomPos(sub);
    const auto this_dist = dist2(a, fct);
    if (this_dist > ftf_dist)
    {
      ftf = a;
      ftf_dist = this_dist;
    }
  }
  // Precalculate the distances of heavy atoms to the reference points, given that subsets[1 to 4] are subsets of subsets[0].
  for (size_t ref = 0; ref < s_Num_references; ++ref)
  {
    const auto& reference = references[ref];
    auto& distp = dista[ref];
    distp.resize(num_points);
    for (size_t p = 0; p < num_points; ++p)
    {
      distp[subset0[p]] = sqrt(dist2(conformer.getAtomPos(subset0[p]), reference));
    }
  }
  // loop over pharmacophoric subsets and reference points.
  size_t qo = 0;
  for (const auto& subset : subsets)
  {
    const auto n = subset.size();
    for (size_t ref = 0; ref < s_Num_references; ++ref)
    {
      // Load distances from precalculated ones
      const auto& distp = dista[ref];
      vector<float> dists(n);
      for (size_t a = 0; a < n; ++a)
      {
        dists[a] = distp[subset[a]];
      }
      // Compute moments
      array<float, 3> m{};
      if (n > 2)
      {
        const auto v = 1.0 / n;
        for (size_t j = 0; j < n; ++j)
        {
          const auto d = dists[j];
          m[0] += d;
        }
        m[0] *= v;
        for (size_t j = 0; j < n; ++j)
        {
          const auto d = dists[j] - m[0];
          m[1] += d * d;
        }
        m[1] = sqrt(m[1] * v);
        for (size_t j = 0; j < n; ++j)
        {
          const auto d = dists[j] - m[0];
          m[2] += d * d * d;
        }
        m[2] = cbrt(m[2] * v);
      }
      else if (n == 2)
      {
        m[0] = 0.5 * (dists[0] + dists[1]);
        m[1] = 0.5 * fabs(dists[0] - dists[1]);
      }
      else if (n == 1)
      {
        m[0] = dists[0];
      }
      for (const auto e : m)
      {
        features[qo++] = e;
      }
    }
  }
  return features;
}

// function to remove whitespaces
static inline auto stripWhiteSpaces(string& str)
{
  str.erase(remove(str.begin(), str.end(), ' '), str.end());
}

template<typename T>
static inline void write_to_binary(T& buf, boost_ofstream& ofs)
{
  // pass the ofstream object by reference, to avoid declaring it each time
  const size_t num_bytes = sizeof(buf);
  ofs.write(reinterpret_cast<char*>(buf.data()), num_bytes);
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
  const path only_smiles_txt = string(smi_file).substr(0, p) + "_only_smiles.txt";
  const path only_id_txt = string(smi_file).substr(0, p) + "_only_id.txt";
  p = string(conformers_file).find('.');
  const path rfprop_file = string(conformers_file).substr(0, p) + "_4properties.f32";
  const path riprop_file = string(conformers_file).substr(0, p) + "_5properties.i16";
  const path usrcatf64_file = string(conformers_file).substr(0, p) + "_usrcat.f64";

  // Inialize the output files
  boost_ofstream conf_file(conformers_file, ios::app);
  boost_ofstream smifile(smi_file, ios::app);
  boost_ofstream smilesfile(only_smiles_txt, ios::app);
  boost_ofstream id_file(only_id_txt, ios::app);
  boost_ofstream rfprop(rfprop_file, ios::binary | ios::app);
  boost_ofstream riprop(riprop_file, ios::binary | ios::app);
  boost_ofstream usrcatf64(usrcatf64_file, ios::binary | ios::app);

  // Initalize constants
  const string pdbqt_extension = ".pdbqt";

  // Initialize variables
  string line, compound, smiles, next_line;
  int pos;
  vector<thread> thread_pool(3); // thread pool with 3 threads
  array<float, 60> features;
  array<float, 4> realfprop;
  array<int16_t, 5> realiprop;
  size_t processed_ligand = 0;
  size_t all_ligands = 0;
  EmbedParameters params(srETKDGv3);
  params.randomSeed = 209;
  params.numThreads = 8;
  params.useRandomCoords = true; // this parameter is used to make the process faster
  params.maxIterations = 3; // max iterations to 1
  SDWriter writer(&conf_file);

  // Search for pdbqt files into the pdbqt folder
  for (const auto& entry : recursive_directory_iterator(pdbqt_folder))
  {
    if (entry.path().extension() == pdbqt_extension)
    {
      std::cout << "Processing molecule Number° " << processed_ligand << '/' << all_ligands << endl;
      std::cout << "Path of the molecule " << entry.path() << endl;
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
              std::cout << "No problem with the SMILES!" << endl;
            }
            else
            {
              cerr << "The SMILES is splited, trying to fix the SMILES" << endl;
              smiles = smiles + next_line;
              std::cout << "The fixed SMILES: " << smiles << endl;
            }
          }
          if (smiles.find('q') != string::npos || smiles.find('r') != string::npos || smiles.find('s') != string::npos)
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
                // generate 4 chemical properties for the molecules
                realfprop[0] = calcExactMW(mol);
                realfprop[1] = calcClogP(mol);
                realfprop[2] = calcTPSA(mol);
                realfprop[3] = calcLabuteASA(mol);
                // add this task to the thread inside the threads vector
                thread_pool.emplace_back([&]() {
                  write_to_binary<array<float, 4>>(realfprop, rfprop);
                });
                //write_to_binary<array<float, 4>>(realfprop, rfprop_file);

                // generate 5 chemical properties for the molecules
                realiprop[0] = mol.getNumHeavyAtoms();
                realiprop[1] = calcNumHBD(mol);
                realiprop[2] = calcNumHBA(mol);
                realiprop[3] = calcNumRotatableBonds(mol);
                realiprop[4] = calcNumRings(mol);
                // add this task to the thread inside the threads vector
                thread_pool.emplace_back([&]() {
         	  write_to_binary<array<int16_t, 5>>(realiprop, riprop);
                });

                std::cout << confIds.size() << " Conformers of " << compound << '\t' << smiles << " are succefully generated!" << endl;
                id_file << compound << '\n';
                smilesfile << smiles << '\n';
                smifile << compound << '\t' << smiles << '\n';
                // Writing conformers in the output SDF
                for (const auto confId : confIds)
                {
                  writer.write(mol, confId);
                }

                // generate the USRCAT features for each conformer
                thread_pool.emplace_back([&]() {
                  for (int i = 0; i < 4; i++)
                  {
                    features = usrcat_features(mol, i);
                    write_to_binary<array<float, 60>>(features, usrcatf64);
                  }
                });

                /*for (int i = 0; i < 4; i++)
                {
                  features = usrcat_features(mol, i);
                  write_to_binary<array<float, 60>>(features, usrcatf64);
                }*/
                processed_ligand++;
                break;
              }
              else
              {
                std::cout << "Molecule doesn't have 4 conformers!" << endl;
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
      // join the threads
      for (thread& th : thread_pool)
      {
        th.join();
      }
      thread_pool.clear();
      all_ligands++;
    }
  }
  std::cout << "Process completed! for " << counter << " compounds" << endl;
}
