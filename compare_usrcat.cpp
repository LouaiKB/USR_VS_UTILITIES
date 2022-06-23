#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <thread>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolDrawing/DrawingToSVG.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <mongocxx/instance.hpp>
#include <mongocxx/pool.hpp>
#include <mongocxx/client.hpp>
#include <bsoncxx/json.hpp>
#include "safe_counter.hpp"
#include "io_service_pool.hpp"

using namespace std;
using namespace std::chrono;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDDepict;
using namespace RDKit::Drawing;
using namespace RDKit::MorganFingerprints;
using namespace RDGeom;
using namespace RDNumeric::Alignments;
using namespace MolTransforms;
using namespace boost::filesystem;
using namespace boost::gregorian;
using namespace boost::posix_time;
using namespace mongocxx;
using bsoncxx::builder::basic::kvp;


inline static auto local_time()
{
	return to_simple_string(microsec_clock::local_time()) + " ";
}

template <typename T>
inline vector<T> read(const path src)
{
	boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	std::cout << local_time() << "Reading " << src << " of " << num_bytes << " bytes" << endl;
	vector<T> buf;
	buf.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(buf.data()), num_bytes);
	return buf;
}

template <typename size_type>
class header_array
{
public:
	explicit header_array(path src)
	{
		src.replace_extension(".ftr");
		boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
		const size_t num_bytes = ifs.tellg();
		std::cout << local_time() << "Reading " << src << " of " << num_bytes << " bytes" << endl;
		hdr.resize(1 + num_bytes / sizeof(size_type));
		hdr.front() = 0;
		ifs.seekg(0);
		ifs.read(reinterpret_cast<char*>(hdr.data() + 1), num_bytes);
	}

	size_t size() const
	{
		return hdr.size() - 1;
	}

protected:
	vector<size_type> hdr;
};

template <typename size_type>
class string_array : public header_array<size_type>
{
public:
	explicit string_array(const path src) : header_array<size_type>(src)
	{
		boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
		const size_t num_bytes = ifs.tellg();
		std::cout << local_time() << "Reading " << src << " of " << num_bytes << " bytes" << endl;
		buf.resize(num_bytes);
		ifs.seekg(0);
		ifs.read(const_cast<char*>(buf.data()), num_bytes);
	}

	string operator[](const size_t index) const
	{
		const auto pos = this->hdr[index];
		const auto len = this->hdr[index + 1] - pos;
		return buf.substr(pos, len);
	}


protected:
	string buf;
};

template <typename size_type>
class stream_array : public header_array<size_type>
{
public:
	explicit stream_array(const path src) : header_array<size_type>(src), ifs(src, ios::binary)
	{
	}

	string operator[](const size_t index)
	{
		const auto pos = this->hdr[index];
		const auto len = this->hdr[index + 1] - pos;
		string buf;
		buf.resize(len);
		ifs.seekg(pos);
		ifs.read(const_cast<char*>(buf.data()), len);
		return buf;
	}

protected:
	boost::filesystem::ifstream ifs;
};

// To compute distance between two atoms
template<typename T>
double dist2(const T& p0, const T& p1)
{
	const auto d0 = p0[0] - p1[0];
	const auto d1 = p0[1] - p1[1];
	const auto d2 = p0[2] - p1[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

array<Point3D, 4> calcRefPoints(const ROMol& mol, const vector<int>& heavyAtoms)
{
	const auto num_points = heavyAtoms.size();
	assert(num_points == mol.getNumHeavyAtoms());
	const auto& conf = mol.getConformer();
	array<Point3D, 4> refPoints;
	for (auto& ref : refPoints)
	{
		assert(ref[0] == 0);
		assert(ref[1] == 0);
		assert(ref[2] == 0);
	}
	auto& ctd = refPoints[0];
	auto& cst = refPoints[1];
	auto& fct = refPoints[2];
	auto& ftf = refPoints[3];
	for (const auto i : heavyAtoms)
	{
		const auto& a = conf.getAtomPos(i);
		ctd += a;
	}
	ctd /= num_points;
	double cst_dist = numeric_limits<double>::max();
	double fct_dist = numeric_limits<double>::lowest();
	double ftf_dist = numeric_limits<double>::lowest();
	for (const auto i : heavyAtoms)
	{
		const auto& a = conf.getAtomPos(i);
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
	for (const auto i : heavyAtoms)
	{
		const auto& a = conf.getAtomPos(i);
		const auto this_dist = dist2(a, fct);
		if (this_dist > ftf_dist)
		{
			ftf = a;
			ftf_dist = this_dist;
		}
	}
	return refPoints;
}

int main(int argc, char* argv[])
{
	// Check the required number of command line arguments.
    if (argc != 2)
    {
	std::cout << "./compare [QUERY]" << endl;
	return 0;
    }

    const auto queries = argv[1];
    const size_t usr0 = 1;
    const size_t usr1 = 0;
    const size_t qnu0 = 60;
    const size_t qnu1 = 12;

    // Initialize constants.
    std::cout << local_time() << "Initializing" << endl;
    // const auto collection = "istar.usr2";
    const size_t num_usrs = 2;
    const array<string, 2> usr_names{{ "USR", "USRCAT" }};
    constexpr array<size_t, num_usrs> qn{{ 12, 60 }};
    constexpr array<double, num_usrs> qv{{ 1.0 / qn[0], 1.0 / qn[1] }};
    const size_t num_refPoints = 4;
    const size_t num_subsets = 5;
    const array<string, num_subsets> SubsetSMARTS
    {{
    	"[!#1]", // heavy
    	"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
    	"[a]", // aromatic
	"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
	"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
    }};
    const size_t num_hits = 100;

    // Wrap SMARTS strings to RWMol objects.
    array<unique_ptr<ROMol>, num_subsets> SubsetMols;
    for (size_t k = 0; k < num_subsets; ++k)
    {
    	SubsetMols[k].reset(reinterpret_cast<ROMol*>(SmartsToMol(SubsetSMARTS[k])));
    }

    // Read ZINC ID file.
    const string_array<size_t> zincids("~/Documents/USR-VS/USR_FILES/final/CHUNK4/ligand_BJ_id.txt");
    const auto num_ligands = zincids.size();
    std::cout << local_time() << "Found " << num_ligands << " database molecules" << endl;

    // Read SMILES file.
    const string_array<size_t> smileses("~/Documents/USR-VS/USR_FILES/final/CHUNK4/ligand_BJ_smiles.txt");
    assert(smileses.size() == num_ligands);
    const auto num_smiles = smileses.size();
    std::cout << local_time() << "Found " << num_smiles << " database smiles molecules" << endl;


    // Read property files of floating point types and integer types.
    const auto zfproperties = read<array<float, 4>>("~/Documents/USR-VS/USR_FILES/final/CHUNK4/ligand_BJ_4.f32");
    assert(zfproperties.size() == num_ligands);
    std::cout << "Number of f properties = " << zfproperties.size() << endl;
    const auto ziproperties = read<array<int16_t, 5>>("~/Documents/USR-VS/USR_FILES/final/CHUNK4/ligand_BJ_5.i16");
    assert(ziproperties.size() == num_ligands);
    std::cout << "Number of i properties = " << ziproperties.size() << endl;

    // Read cumulative number of conformers file.
    const auto mconfss = read<size_t>("~/Documents/USR-VS/USR_FILES/final/CHUNK4/mconfs.f64");
    const auto num_conformers = mconfss.back();
	assert(mconfss.size() == num_ligands);
    assert(num_conformers >= num_ligands);
    std::cout << local_time() << "Found " << num_conformers << " database conformers" << endl;

    // Read feature file.
    const auto features = read<array<int16_t, qn.back()>>("~/Desktop/usrcat.i16");
    cout << "SIZE OF FZATURES : " << features.size() << endl;
    assert(features.size() == num_ligands * 4);
    std::cout << local_time() << "Found " << features.size() << " usrcat features" << endl;
    std::cout << " features size = " << features.size() << " num ligands = " << num_ligands * 4 << endl;
    // Read ligand footer file and open ligand SDF file for seeking and reading.
    stream_array<size_t> ligands("/home/kassa-baghdouche/Documents/USR-VS/USR_FILES/final/CHUNK4/ligand_BJ.sdf");
    std::cout << local_time() << "Found " << ligands.size() << " conformers!" << endl;
    assert(ligands.size() == features.size());
    // Initialize variables.
    array<vector<int>, num_subsets> subsets;
    array<vector<double>, num_refPoints> dista;
    alignas(32) array<double, qn.back()> q;

    // Initialize vectors to store compounds' primary score and their corresponding conformer.
    vector<double> scores(num_ligands); // Primary score of molecules.
    vector<size_t> cnfids(num_ligands); // ID of conformer with the best primary score.
    const auto compare = [&](const size_t val0, const size_t val1) // Sort by the primary score.
    {
    	return scores[val0] < scores[val1];
    };

    // Initialize an io service pool and create worker threads for later use.
    const size_t num_threads = thread::hardware_concurrency();
    std::cout << local_time() << "Creating an io service pool of " << num_threads << " worker threads" << endl;
    io_service_pool io(num_threads);
    safe_counter<size_t> cnt;

    // Initialize the number of chunks and the number of molecules per chunk.
    const auto num_chunks = num_threads << 4;
    // const auto num_chunks = num_threads << 2;
    const auto chunk_size = 1 + (num_ligands - 1) / num_chunks;
    assert(chunk_size * num_chunks >= num_ligands);
    // assert(chunk_size >= num_hits);
    std::cout << local_time() << "Using " << num_chunks << " chunks and a chunk size of " << chunk_size << endl;
    vector<size_t> scase(num_ligands);
    vector<size_t> zcase(num_hits * (num_chunks - 1) + min(num_hits, num_ligands - chunk_size * (num_chunks - 1))); // The last chunk might have fewer than num_hits records.

    // Enter event loop.
    std::cout << local_time() << "Entering event loop" << endl;
    std::cout.setf(ios::fixed, ios::floatfield);
       
    SDMolSupplier sup(queries, true, false, true);
    unsigned int query_number = 1;

    boost::filesystem::ofstream hits_csv("hits_float_short.csv");
    hits_csv.setf(ios::fixed, ios::floatfield);
    hits_csv << "Query,Hit,USRCAT SCORE FLOAT rdkit" << '\n';
    string last_compound = "";
    while (!sup.atEnd())
    {
        const unique_ptr<ROMol> mol_ptr(sup.next());
        auto& qryMol = *mol_ptr;
        const auto num_atoms = qryMol.getNumAtoms();
        const auto num_heavy_atoms = qryMol.getNumHeavyAtoms();
        std::cout << local_time() << "Found " << num_atoms << " atoms and " << num_heavy_atoms << " heavy atoms" << endl;
        const string compound = qryMol.getProp<string>("_Name");
        if (last_compound == compound) continue;
        else 
        {
            last_compound = compound;
            cout << "Processing query N° " << query_number++ << endl;
            if (query_number == 10000) break;
        }
        
        // Classify atoms to pharmacophoric subsets.
        std::cout << local_time() << "Classifying atoms into subsets" << endl;
        for (size_t k = 0; k < num_subsets; ++k)
        {
            vector<vector<pair<int, int>>> matchVect;
            RDKit::SubstructMatch(qryMol, *SubsetMols[k], matchVect);
            const auto num_matches = matchVect.size();
            auto& subset = subsets[k];
            subset.resize(num_matches);
            for (size_t i = 0; i < num_matches; ++i)
            {
                subset[i] = matchVect[i].front().second;
            }
            std::cout << local_time() << "Found " << num_matches << " atoms for subset " << k << endl;
        }
        const auto& subset0 = subsets.front();
        assert(subset0.size() == num_heavy_atoms);
        // Calculate the four reference points.
        std::cout << local_time() << "Calculating " << num_refPoints << " reference points" << endl;
        const auto qryRefPoints = calcRefPoints(qryMol, subset0);
        const Point3DConstPtrVect qryRefPointv
        {{
            &qryRefPoints[0],
            &qryRefPoints[1],
            &qryRefPoints[2],
            &qryRefPoints[3],
        }};
        // Precalculate the distances of heavy atoms to the reference points, given that subsets[1 to 4] are subsets of subsets[0].
        std::cout << local_time() << "Calculating " << num_heavy_atoms * num_refPoints << " pairwise distances" << endl;
        const auto& qryCnf = qryMol.getConformer();
        for (size_t k = 0; k < num_refPoints; ++k)
        {
            const auto& refPoint = qryRefPoints[k];
            auto& distp = dista[k];
            distp.resize(num_atoms);
            for (size_t i = 0; i < num_heavy_atoms; ++i)
            {
                distp[subset0[i]] = sqrt(dist2(qryCnf.getAtomPos(subset0[i]), refPoint));
            }
        }
        // Loop over pharmacophoric subsets and reference points.
        // Generate USRCAT features
        std::cout << local_time() << "Calculating " << 3 * num_refPoints * num_subsets << " moments of USRCAT feature" << endl;
        size_t qo = 0;
        for (const auto& subset : subsets)
        {
            const auto n = subset.size();
            for (size_t k = 0; k < num_refPoints; ++k)
            {
                // Load distances from precalculated ones.
                const auto& distp = dista[k];
                vector<double> dists(n);
                for (size_t i = 0; i < n; ++i)
                {
                    dists[i] = distp[subset[i]];
                }

                // Compute moments.
                array<double, 3> m{};
                if (n > 2)
                {
                    const auto v = 1.0 / n;
                    for (size_t i = 0; i < n; ++i)
                    {
                        const auto d = dists[i];
                        m[0] += d;
                    }
                    m[0] *= v;
                    for (size_t i = 0; i < n; ++i)
                    {
                        const auto d = dists[i] - m[0];
                        m[1] += d * d;
                    }
                    m[1] = sqrt(m[1] * v);
                    for (size_t i = 0; i < n; ++i)
                    {
                        const auto d = dists[i] - m[0];
                        m[2] += d * d * d;
                    }
                    m[2] = cbrt(m[2] * v);
                }
                else if (n == 2)
                {
                    m[0] = 0.5 *     (dists[0] + dists[1]);
                    m[1] = 0.5 * fabs(dists[0] - dists[1]);
                }
                else if (n == 1)
                {
                    m[0] = dists[0];
                }
                for (const auto e : m)
                {
                    q[qo++] = e;
                }
            }
        }
        std::cout << local_time() << "Calculating " << num_ligands << " " << usr_names[usr0] << " scores" << endl;
        scores.assign(scores.size(), numeric_limits<double>::max());
        std::iota(scase.begin(), scase.end(), 0);
        cnt.init(num_chunks);
        for (size_t l = 0; l < num_chunks; ++l)
        {
            io.post([&,l]()
            {
                // Loop over molecules of the current chunk.
                const auto chunk_beg = chunk_size * l;
                const auto chunk_end = min(chunk_beg + chunk_size, num_ligands);
                for (size_t k = chunk_beg; k < chunk_end; ++k)
                {
                    // Loop over conformers of the current molecule and calculate their primary score.
                    auto& scorek = scores[k];
                    size_t j = k ? mconfss[k - 1] : 0;
                    for (const auto mconfs = mconfss[k]; j < mconfs; ++j)
                    {
                        const auto& d = features[j];
                        double s = 0;
                        for (size_t i = 0; i < qnu0; ++i)
                        {
                            float us = (float)d[i] / 1000;
                            float sc = *((float*)&us);
                            s += abs(q[i] - us);
                            if (s >= scorek) break;
                        }
                        if (s < scorek)
                        {
                            scorek = s;
                            cnfids[k] = j;
                        }
                    }
                }
                // Sort the scores of molecules of the current chunk.
                std::sort(scase.begin() + chunk_beg, scase.begin() + chunk_end, compare);

                // Copy the indexes of top hits of the current chunk to a global vector for final sorting.
                std::copy_n(scase.begin() + chunk_beg, min(num_hits, chunk_end - chunk_beg), zcase.begin() + num_hits * l);

                cnt.increment();
            });
        }
        cnt.wait();
        // Sort the top hits from chunks.
        std::cout << local_time() << "Sorting " << zcase.size() << " hits by " << usr_names[0] << " score" << endl;
        std::sort(zcase.begin(), zcase.end(), compare);
        // Obtain indexes to the first hit molecule and the hit conformer.
        const auto k = zcase[0];
        const auto j = cnfids[k];

        // Read SDF content of the hit conformer.
        const auto lig = ligands[j];

        // Construct a RDKit ROMol object.
        istringstream iss(lig);
        SDMolSupplier sup(&iss, false, true, false, true);
        assert(sup.length() == 1);
        assert(sup.atEnd());
        const unique_ptr<ROMol> hit_ptr(sup.next());
        auto& hitMol = *hit_ptr;

        // Find heavy atoms.
        vector<vector<pair<int, int>>> matchVect;
        RDKit::SubstructMatch(hitMol, *SubsetMols[0], matchVect);
        const auto num_matches = matchVect.size();
        assert(num_matches == hitMol.getNumHeavyAtoms());
        vector<int> hitHeavyAtoms(num_matches);
        for (size_t i = 0; i < num_matches; ++i)
        {
            hitHeavyAtoms[i] = matchVect[i].front().second;
            assert(hitHeavyAtoms[i] == i); // hitHeavyAtoms can be constructed using iota(hitHeavyAtoms.begin(), hitHeavyAtoms.end(), 0); because for RDKit-generated SDF molecules, heavy atom are always the first few atoms.
        }

        // Calculate the four reference points.
        const auto hitRefPoints = calcRefPoints(hitMol, hitHeavyAtoms);
        const Point3DConstPtrVect hitRefPointv
        {{
            &hitRefPoints[0],
            &hitRefPoints[1],
            &hitRefPoints[2],
            &hitRefPoints[3],
        }};

        // Calculate the secondary score of the saved conformer, which has the best primary score.
        const auto& d = features[j];
        double s = 0;
        for (size_t i = 0; i < qnu1; ++i)
        {
            s += abs(q[i] - d[i]);
        }

        const auto u0score = 1 / (1 + scores[k] * qv[usr0]); // Primary score of the current molecule.
        const auto u1score = 1 / (1 + s         * qv[usr1]); // Secondary score of the current molecule.
        size_t pos = zincids[k].find('\n');
        const auto zincid = zincids[k].substr(0, pos);
        hits_csv
            << compound
            << ',' << zincid
            << ',' << (usr1 ? u1score : u0score)
            << '\n'
        ;
    }
}

