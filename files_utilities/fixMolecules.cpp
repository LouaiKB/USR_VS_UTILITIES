#include <iostream>
#include <string>
#include <array>
#include <algorithm>
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
void print_arr(const T& arr)
{
	cout << '[';
	for (size_t i = 0; i < arr.size(); i++)
	{
		(i != arr.size() - 1) ? cout << arr[i] << ',' : cout << arr[i] << ']' << endl;
	}
}


int main(int argc, char* argv[])
{
	array<float, 4> realfprop;
	const path file = argv[1];
	const string ext = ".pdbqt";
	long counter = 0;
	string line;
	const array<string, 20> elements {{
		"PV-001925004350_3_T1.pdbqt",
		"PV-001924960689_3_T1.pdbqt",
		"PV-001924547224_2_T1.pdbqt",
		"PV-001925088837_2_T1.pdbqt",
		"PV-001924761490_3_T1.pdbqt",
		"PV-001924933847_12_T1.pdbqt",
		"PV-001925000277_1_T1.pdbqt",
		"PV-001924544716_1_T1.pdbqt",
		"PV-001924761104_4_T1.pdbqt",
		"PV-001924415292_6_T1.pdbqt",
		"PV-001924961746_1_T1.pdbqt",
		"PV-001925003222_2_T1.pdbqt",
		"PV-001924547164_2_T1.pdbqt",
		"PV-001924545172_4_T1.pdbqt",
		"PV-001924655398_2_T1.pdbqt",
		"PV-001924544931_1_T1.pdbqt",
		"PV-001925003198_2_T1.pdbqt",
		"PV-001924561317_1_T1.pdbqt",
		"PV-001924928144_3_T1.pdbqt",
		"PV-001924561275_1_T1.pdbqt"
	}};
	for (const auto& entry : recursive_directory_iterator(file))
	{
		if (entry.path().extension() == ext && find(begin(elements), end(elements), entry.path().filename()) != end(elements))
		{
			boost_ifstream ifs(entry.path());
			while (getline(ifs, line))
			{
				if (line.find("SMILES:") != string::npos)
				{
					int p = line.find(':');
					string smiles = line.substr(p + 2);
					const unique_ptr<ROMol> mol_ptr(SmilesToMol(smiles));
					const unique_ptr<ROMol> molhs(addHs(*mol_ptr));
					auto& mol = *molhs;
					realfprop[0] = calcExactMW(mol);
					realfprop[1] = calcClogP(mol);
					realfprop[2] = calcTPSA(mol);
					realfprop[3] = calcLabuteASA(mol);
					print_arr(realfprop);
					break;
				}
			}
		}
		counter++;
	}
}