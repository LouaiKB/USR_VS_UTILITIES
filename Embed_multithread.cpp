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
#include <string>
#include <array>
#include <vector>
#include <future>
#include <mutex>
#include <chrono>
#include <thread>
#include <condition_variable>
#include <boost/asio/io_service.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/MolDescriptors.h>

using namespace std;
using namespace std::chrono;
using namespace RDKit;
using namespace RDGeom;
using namespace RDKit::Descriptors;
using namespace boost::asio;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;
using namespace boost::filesystem;
using boost_ifstream = boost::filesystem::ifstream;
using boost_ofstream = boost::filesystem::ofstream;

class io_service_pool : public io_service, public vector<future<void>>
{
public:
  //! Creates a number of threads to listen to the post event of an io service.
  explicit io_service_pool(const unsigned concurrency) : w(make_unique<work>(*this))
  {
    reserve(concurrency);
    for (unsigned i = 0; i < concurrency - 5; ++i)
    {
      emplace_back(async(launch::async, [&]()
      {
        run();
      }));
    }
  }

  //! Waits for all the posted work and created threads to complete, and propagates thrown exceptions if any.
  void wait()
  {
    w.reset();
    for (auto& f : *this)
    {
      f.get();
    }
  }
private:
	unique_ptr<work> w; //!< An io service work object, resetting which to nullptr signals the io service to stop receiving additional work.
};

//! Represents a thread safe counter.
template <typename T>
class safe_counter
{
public:
  //! Initializes the counter to 0 and its expected hit value to z.
  void init(const T z)
  {
    n = z;
    i = 0;
  }

  //! Increments the counter by 1 in a thread safe manner, and wakes up the calling thread waiting on the internal mutex.
  void increment()
  {
    lock_guard<mutex> guard(m);
    if (++i == n) cv.notify_one();
  }

  //! Waits until the counter reaches its expected hit value.
  void wait()
  {
    unique_lock<mutex> lock(m);
    if (i < n)
      cv.wait(lock);
  }

private:
        mutex m;
	condition_variable cv;
	T n; //!< Expected hit value.
	T i; //!< Counter value.
};


template<typename T>
float dist2(const T& p0, const T& p1)
{
  const auto d0 = p0[0] - p1[0];
  const auto d1 = p0[1] - p1[1];
  const auto d2 = p0[2] - p1[2];
  return d0 * d0 + d1 * d1 + d2 * d2;
}

array<float, 60> usrcat_generator(ROMol& mol, int index)
{
  const size_t num_references = 4;
  const size_t num_subsets = 5;
  const array<string, 5> SubsetSMARTS
  {{
    "[!#1]", // heavy
    "[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
    "[a]", // aromatic
    "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
    "[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
  }};
  array<unique_ptr<ROMol>, 5> SubsetMols;
  array<vector<int>, 5> subsets;
  array<Point3D, 4> references;
  array<vector<float>, 4> dista;
  array<float, 60> features;
  // Wrap SMARTS strings to ROMol objects.
  for (size_t k = 0; k < num_subsets; ++k)
  {
    SubsetMols[k].reset(reinterpret_cast<ROMol*>(SmartsToMol(SubsetSMARTS[k])));
  }
  const auto num_points = mol.getNumHeavyAtoms();
  const auto& conformer = mol.getConformer(index);
  for (size_t k = 0; k < num_subsets; ++k)
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
  for (size_t ref = 0; ref < num_references; ++ref)
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
    for (size_t ref = 0; ref < num_references; ++ref)
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

int main(int argc, char* argv[])
{
  if (argc != 4)
  {
    cerr << "Usage: ./EmbedAllAndGenerateSMIL [PDBQT FOLDER] [CONFORMERS SDF] [NUM_PDBQT]" << endl;
    return 1;
  }
  const auto started = system_clock::now();
  cout << "started!" << endl;

  // Obtain the files from the CLI argument
  const auto pdbqt_folder = argv[1];
  const string conformers_file = argv[2];
  const auto num_of_pdbqt = atoi(argv[3]);

  // create strings file from the input
  int p = conformers_file.find('.');
  string only_smiles = conformers_file.substr(0, p) + "_only_smiles.txt";
  string only_id = conformers_file.substr(0, p) + "_only_id.txt";
  string rfprop_file = conformers_file.substr(0, p) + "_4properties.f32";
  string riprop_file = conformers_file.substr(0, p) + "_5properties.i16";
  string usrcat_file = conformers_file.substr(0, p) + "_usrcat.f64";

  // Initalize constant
  const string pdbqt_extension = ".pdbqt";
  const size_t num_threads = thread::hardware_concurrency();
  const auto num_chunks = num_threads << 4;

  // Initialize an io service pool and create worker threads for later use
  io_service_pool io(num_threads);
  safe_counter<size_t> cnt;

  // Initialize array of properties
  array<float, 4> realfprop;
  array<int16_t, 5> realiprop;
  array<float, 60> usrcat_features;

  // Initialize vectors
  vector<path> pdbqts_vector;
  vector<future<void>> thread_pool;
  cout << "reservation .." << endl;
  pdbqts_vector.reserve(num_of_pdbqt);

  // Initializing of ofstreams
  boost_ofstream ofs(conformers_file, ios::app);
  boost_ofstream ofs_smi(only_smiles, ios::app);
  boost_ofstream ofs_id(only_id, ios::app);
  boost_ofstream ofs_rfprop(rfprop_file, ios::binary | ios::app);
  boost_ofstream ofs_riprop(riprop_file, ios::binary | ios::app);
  boost_ofstream ofs_usrcat(usrcat_file, ios::binary | ios::app);

  // Initialize SDwriter objects to write conformers to sdf
  SDWriter writer(&ofs);

  // Search for pdbqt files into the pdbqt folder
  cout << "Processing the vector ..." << endl;
  for (const auto& entry : recursive_directory_iterator(pdbqt_folder))
  {
    if (entry.path().extension() == pdbqt_extension)
    {
      pdbqts_vector.emplace_back(entry.path());
    }
  }

  mutex mu;
  size_t lig_count = 0;
  size_t num_files = pdbqts_vector.size();
  const size_t chunk_size = num_files / num_chunks;
  
  // Generate conformers and write properties in files in a multithreading way
  cnt.init(num_chunks);
  for (size_t l = 0; l < num_chunks; ++l)
  {
    io.post([&, l] ()
    {
      const auto chunk_beg = chunk_size * l;
      const auto chunk_end = chunk_size * l + chunk_size;
      string line, compound, smiles, next_line;
      int pos;
      EmbedParameters params(srETKDGv3);
      params.randomSeed = 209;
      params.numThreads = 8;
      params.useRandomCoords = true;
      params.maxIterations = 10;
      for (size_t counter = chunk_beg; counter < chunk_end; ++counter)
      {
        boost_ifstream pdbqt_ifs(pdbqts_vector[counter]);
        while (getline(pdbqt_ifs, line))
        {
          if (line.find("Compound:") != string::npos)
          {
            pos = line.find(':');
            compound = line.substr(pos + 2);
          }
          if (line.find("SMILES:") != string::npos)
          {
            pos = line.find(':');
            smiles = line.substr(pos + 2);
            if (getline(pdbqt_ifs, next_line))
            {
              if (next_line.find("REMARK") == string::npos)
              {
                // smiles is splitted, trying to fix it
                smiles = smiles + next_line;
              }
            }
            if (smiles.find('q') != string::npos || smiles.find('r') != string::npos || smiles.find('s') != string::npos)
            {
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
                  break;
                
                if (confIds.size() == 4)
                {
                  realfprop[0] = calcExactMW(mol);
                  realfprop[1] = calcClogP(mol);
                  realfprop[2] = calcTPSA(mol);
                  realfprop[3] = calcLabuteASA(mol);

                  realiprop[0] = mol.getNumHeavyAtoms();
                  realiprop[1] = calcNumHBD(mol);
                  realiprop[2] = calcNumHBA(mol);
                  realiprop[3] = calcNumRotatableBonds(mol);
                  realiprop[4] = calcNumRings(mol);

                  {
                    lock_guard<mutex> lock(mu);
                    std::cout << "Iteration N° " << counter << " And generated ligands: " << lig_count <<  endl;
                    std::cout << "Processing file: " << pdbqts_vector[counter] << endl;
                    cout << confIds.size() << " Conformers of " << compound << " : " << smiles << " are successfully generated!" << endl;
                    lig_count++;
                    for (const auto confId : confIds) {
                      writer.write(mol, confId);
                    }
                    const size_t num_bytes_realf = sizeof(realfprop);
                    const size_t num_bytes_reali = sizeof(realiprop);
                    ofs_rfprop.write(reinterpret_cast<char*>(realfprop.data()), num_bytes_realf);
                    ofs_riprop.write(reinterpret_cast<char*>(realiprop.data()), num_bytes_reali);
                    ofs_id << compound << '\n';
                    ofs_smi << smiles << '\n';
                    for (int i = 0; i < 4; i++)
                    {
                      usrcat_features = usrcat_generator(mol, i);
                      const size_t num_bytes = sizeof(usrcat_features);
                      ofs_usrcat.write(reinterpret_cast<char*>(usrcat_features.data()), num_bytes);
                    }
                  }

                }
              }
              catch (const AtomValenceException& e)
              {
                cerr << pdbqts_vector[counter] << " molecule raised an atom valence exception" << endl;
                break;
              }
              catch (const AtomKekulizeException& e)
              {
                cerr << pdbqts_vector[counter] << " molecule raised an atom kekulize exception" << endl;
                break;
              }
              catch (const AtomSanitizeException& e)
              {
                cerr << pdbqts_vector[counter] << " molecule raised an atom sanitize exception" << endl;
                break;
              }
            }
          }
        }
      }
      cnt.increment();
    });
  }
  cnt.wait();
  const auto finished = system_clock::now();
  const auto runtime = (finished - started).count() * 1e-9;

  cout << "process finished in " << setprecision(3) << runtime << " seconds" << endl; 

}
