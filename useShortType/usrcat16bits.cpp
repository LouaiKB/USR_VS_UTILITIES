#include <iostream>
#include <array>
#include <vector>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/MolSupplier.h>

using namespace std;
using namespace RDKit;
using namespace RDGeom;
using namespace RDKit::Descriptors;
using namespace RDKit::MolOps;
using USRCAT_DESCRIPTORS = array<int16_t, 60>;

template<typename T>
static inline float dist2(const T& p0, const T& p1) {
  const auto d0 = p0[0] - p1[0];
  const auto d1 = p0[1] - p1[1];
  const auto d2 = p0[2] - p1[2];
  return d0 * d0 + d1 * d1 + d2 * d2;
}

static USRCAT_DESCRIPTORS usrcat_generator(ROMol& mol) {
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
  array<int16_t, 60> features;
  // Wrap SMARTS strings to ROMol objects.
  for (size_t k = 0; k < num_subsets; ++k) {
    SubsetMols[k].reset(reinterpret_cast<ROMol*>(SmartsToMol(SubsetSMARTS[k])));
  }
  const auto num_points = mol.getNumHeavyAtoms();
  const auto& conformer = mol.getConformer();
  for (size_t k = 0; k < num_subsets; ++k) {
    vector<vector<pair<int, int>>> matchVect;
    SubstructMatch(mol, *SubsetMols[k], matchVect);
    const auto num_matches = matchVect.size();
    auto& subset = subsets[k];
    subset.resize(num_matches);
    for (size_t j = 0; j < num_matches; ++j) {
      subset[j] = matchVect[j].front().second;
    }
  }
  const auto& subset0 = subsets.front();

  for (auto& ref : references) {
    ref.x = ref.y = ref.z = 0;
  }
  auto& ctd = references[0];
  auto& cst = references[1];
  auto& fct = references[2];
  auto& ftf = references[3];
  for (const auto sub : subset0) {
    const auto& a = conformer.getAtomPos(sub);
    ctd += a;
  }
  ctd /= num_points;
  float cst_dist = numeric_limits<float>::max();
  float fct_dist = numeric_limits<float>::lowest();
  float ftf_dist = numeric_limits<float>::lowest();
  for (const auto sub : subset0) {
    const auto& a = conformer.getAtomPos(sub);
    const auto this_dist = dist2(a, ctd);
    if (this_dist < cst_dist) {
      cst = a;
      cst_dist = this_dist;
    }
    if (this_dist > fct_dist) {
      fct = a;
      fct_dist = this_dist;
    }
  }
  for (const auto sub : subset0) {
    const auto& a = conformer.getAtomPos(sub);
    const auto this_dist = dist2(a, fct);
    if (this_dist > ftf_dist) {
      ftf = a;
      ftf_dist = this_dist;
    }
  }
  // Precalculate the distances of heavy atoms to the reference points, given that subsets[1 to 4] are subsets of subsets[0].
  for (size_t ref = 0; ref < num_references; ++ref) {
    const auto& reference = references[ref];
    auto& distp = dista[ref];
    distp.resize(num_points);
    for (size_t p = 0; p < num_points; ++p) {
      distp[subset0[p]] = sqrt(dist2(conformer.getAtomPos(subset0[p]), reference));
    }
  }
  // loop over pharmacophoric subsets and reference points.
  size_t qo = 0;
  for (const auto& subset : subsets) {
    const auto n = subset.size();
    for (size_t ref = 0; ref < num_references; ++ref) {
      // Load distances from precalculated ones
      const auto& distp = dista[ref];
      vector<float> dists(n);
      for (size_t a = 0; a < n; ++a) {
        dists[a] = distp[subset[a]];
      }
      // Compute moments
      array<float, 3> m{};
      if (n > 2) {
        const auto v = 1.0 / n;
        for (size_t j = 0; j < n; ++j) {
          const auto d = dists[j];
          m[0] += d;
        }
        m[0] *= v;
        for (size_t j = 0; j < n; ++j) {
          const auto d = dists[j] - m[0];
          m[1] += d * d;
        }
        m[1] = sqrt(m[1] * v);
        for (size_t j = 0; j < n; ++j) {
          const auto d = dists[j] - m[0];
          m[2] += d * d * d;
        }
        m[2] = cbrt(m[2] * v);
      }
      else if (n == 2) {
        m[0] = 0.5 * (dists[0] + dists[1]);
        m[1] = 0.5 * fabs(dists[0] - dists[1]);
      }
      else if (n == 1) {
        m[0] = dists[0];
      }
      for (const auto e : m) {
        int16_t val = e * 1000;
        features[qo++] = val;
      }
    }
  }
  return features;
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    cerr << "./usrcat16bits [INPUT]" << endl;
    return 1;
  }
  const auto conf = argv[1];
  SDMolSupplier sup(conf, true, false, true);
  while (!sup.atEnd()) {
    const unique_ptr<ROMol> mol_ptr(sup.next());
    auto& mol = *mol_ptr;
    cerr << mol.getProp<string>("_Name") << endl;
    auto usrcat = usrcat_generator(mol);
    cout.write(reinterpret_cast<char*>(usrcat.data()), sizeof(USRCAT_DESCRIPTORS));
  }
  return 0;
}