#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FragCatalog/FragFPGenerator.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>

using namespace std;
using namespace RDKit;
using namespace RDGeom;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;
using namespace boost::filesystem;
using boost_ifstream = boost::filesystem::ifstream;
using boost_ofstream = boost::filesystem::ofstream;

int main(int argc, char* argv[])
{
	const auto confs_file = argv[1];
	boost_ofstream ofs(confs_file, ios::app);
	SDWriter writer(&ofs);
	string compound, smiles;
	for (string line; getline(cin, line);)
	{
		cout << line << endl;
		if (line.find("Compound:") != string::npos)
		{
			int pos = line.find(':');
			compound = line.substr(pos + 2);
			cout << compound << endl;
		}
		if (line.find("SMILES:") != string::npos)
		{
			int pos = line.find(':');
			smiles = line.substr(pos + 2);
			cout << smiles << endl;
			break;
		}
	}
	EmbedParameters params(srETKDGv3);
	params.useRandomCoords = true;
	params.randomSeed = 209;
	params.numThreads = 8;
	params.maxIterations = 5;
	const unique_ptr<ROMol> mol_ptr(SmilesToMol(smiles));
	const unique_ptr<ROMol> mol_hs(addHs(*mol_ptr));
	auto& mol = *mol_hs;
	mol.setProp("_Name", compound);
	const auto confIds = EmbedMultipleConfs(mol, 5, params);
	if (confIds.size() == 5)
	{
		cout << "5 conformers" << endl;
		for (int i = 4; i > 1; i--)
		{
			writer.write(mol, i);
		}
	}
	else
	{
		cout << "Molecules doesn't have 5 conformers" << endl;
	}

} 	