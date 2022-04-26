#include <iostream>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace std;
using namespace RDGeom;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;
using namespace RDKit::Descriptors;
using namespace boost::filesystem;
using boost_ofstream = boost::filesystem::ofstream;
using boost_ifstream = boost::filesystem::ifstream;


int main(int argc, char* argv[])
{
	const string sdf_path = argv[1];
	// const path output_confs = "/home/kassa-baghdouche/Documents/USR-VS/USR_FILES/final/ligands_fixed.sdf";
	// boost_ofstream ofs(output_confs, ios::app);
	// SDWriter writer(&cout);
	long counter = 0;
	int16_t c = 0;
	SDMolSupplier sup(sdf_path, true, false, true);
	while (!sup.atEnd())
	{
		// if (c != 4)
		// {
			const unique_ptr<ROMol> qry_ptr(sup.next());
			auto& mol = *qry_ptr;
			string compound = mol.getProp<string>("_Name");
			// if (compound == "PV-001924561275_1_T1")
			// {
				// writer.write(mol);
				// c++;
				counter++;
				// cout << counter << endl;
				// continue;
			// }
			// writer.write(mol);
			// counter++;
			cout << counter << " " << compound << endl;
	}
	// boost_ifstream ifs(sdf_path, ios::ate);
	// ifs.seekg(0);	
	// long counter = 0;
	// for (string line; getline(ifs, line);)
	// {
	// 	if (line == "$$$$")
	// 		counter++;
	// 	cout << counter << endl;
	// }
	cout << "Number of writtern compounds are " << counter << endl;
}