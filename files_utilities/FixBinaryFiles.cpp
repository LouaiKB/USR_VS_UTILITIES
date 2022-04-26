#include <iostream>
#include <string>
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
using namespace RDGeom;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;
using namespace RDKit::Descriptors;
using namespace boost::filesystem;
using boost_ifstream = boost::filesystem::ifstream;
using boost_ofstream = boost::filesystem::ofstream;


template<typename T>
static bool EqualsArr(const T& arr1, T& arr2)
{
	const size_t size = arr1.size();
	for (int i = 0; i < size; ++i)
		if (arr1[i] != arr2[i])
			return false;
	return true;
}

template<typename T>
static inline void write_to_binary(T& buf, const path output_path)
{
  boost_ofstream ofs(output_path, ios::binary | ios::app);
  const size_t num_bytes = sizeof(buf);
  ofs.write(reinterpret_cast<char*>(buf.data()), num_bytes);
}

// function to remove whitespaces
inline auto stripWhiteSpaces(string& str)
{
  str.erase(remove(str.begin(), str.end(), ' '), str.end());
}


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

template <typename T>
inline vector<T> read(const path src)
{
	boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	vector<T> buf;
	buf.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(buf.data()), num_bytes);
	return buf;
}


int main(int argc, char* argv[])
{

  path ofs_usrcat = "/home/kassa-baghdouche/Documents/USR-VS/USR_FILES/final/usrcat.f64";
  path ofs_rfprop = "/home/kassa-baghdouche/Documents/USR-VS/USR_FILES/final/rfprop.f32";
  path ofs_riprop = "/home/kassa-baghdouche/Documents/USR-VS/USR_FILES/final/riprop.i16";

	string smiles, compound, line;
	int pos;
	array<float, 60> usrcat;
	array<float, 4> realfprop;
	array<int16_t, 5> realiprop;
	const auto pdbqt_file = argv[1];
	boost_ifstream ifs(pdbqt_file);
	EmbedParameters params(srETKDGv3);
	params.numThreads = 8;
	params.randomSeed = 209;
	params.useRandomCoords = true;
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
		
			try
			{
				const unique_ptr<ROMol> smi_ptr(SmilesToMol(smiles));
				const unique_ptr<ROMol> mol_ptr(addHs(*smi_ptr));
				auto& mol = *mol_ptr;
				cout << "embediiingé " << endl;
				const auto confIds = EmbedMultipleConfs(mol, 4, params);
				// generate 4 chemical properties for the molecules
				realfprop[0] = calcExactMW(mol);
				realfprop[1] = calcClogP(mol);
				realfprop[2] = calcTPSA(mol);
				realfprop[3] = calcLabuteASA(mol);

				// generate 5 chemical properties for the molecules
				realiprop[0] = mol.getNumHeavyAtoms();
				realiprop[1] = calcNumHBD(mol);
				realiprop[2] = calcNumHBA(mol);
				realiprop[3] = calcNumRotatableBonds(mol);
				realiprop[4] = calcNumRings(mol);
				usrcat = usrcat_features(mol, 3);
				break;
			}
			catch(const std::exception& e)
			{
				std::cerr << e.what() << '\n';
			}
		}
	}
	cout << "fetching the features" << endl;
	vector<array<float, 60>> usrcat_ft = read<array<float, 60>>("/home/kassa-baghdouche/Documents/USR-VS/USR_FILES/final/ligand_usrcat.f64");
	cout << "fetching the 5 properties" << endl;
	vector<array<int16_t, 5>> reali = read<array<int16_t, 5>>("/home/kassa-baghdouche/Documents/USR-VS/USR_FILES/final/ligand_5properties.i16");
	cout << "fetching the 4 properties" << endl;
	vector<array<float, 4>> realf = read<array<float, 4>>("/home/kassa-baghdouche/Documents/USR-VS/USR_FILES/final/ligand_4properties.f32");

	vector<array<float, 60>> new_usrcat;
	vector<array<int16_t, 5>> new_realiprop;
	vector<array<float, 4>> new_realfprop;

	for (const auto& arr : usrcat_ft)
	{
		if (EqualsArr(arr, usrcat))
		{
			cout << "features found!" << endl;
			new_usrcat.push_back(arr);
			break;
		}
		else
			new_usrcat.push_back(arr);
	}
	int coun = 0;
	for (const auto& a : reali)
	{
		if (EqualsArr(a, realiprop))
		{
			cout << "properties found" << endl;
      if (coun <= 1238)
        new_realiprop.push_back(a);
      else
        break;
			coun++;
		}
    else
      new_realiprop.push_back(a);
	}
	cout << "smae 5 prop features = " << coun << endl;
	coun = 0;
	for (const auto& a : realf)
	{
		if (EqualsArr(a, realfprop))
		{
			cout << "properties 4 found!" << endl;
      if (coun <= 9)
        new_realfprop.push_back(a);
      else
        break;
      
			coun++;
		}
    else
      new_realfprop.push_back(a);
	}
	cout << "same 4 prop = " << coun << endl;
  cout << "size usrcat = " << new_usrcat.size() << endl;
  cout << "size irpop = " << new_realiprop.size() << endl;
  cout << "size fprop = " << new_realfprop.size() << endl;

  for (auto& arr : new_usrcat)
  {
    write_to_binary<array<float, 60>>(arr, ofs_usrcat);
  }
  for (auto& arr : new_realiprop)
  {
    write_to_binary(arr, ofs_riprop);
  }
  for (auto& arr : new_realfprop)
  {
    write_to_binary(arr, ofs_rfprop);
  }
}