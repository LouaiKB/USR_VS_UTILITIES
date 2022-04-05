/**
 * @file encodeDescriptorsf32.cpp
 * @author Louai KASSA BAGHDOUCHE 
 * @brief generate chemical descriptors from smi files
 * @version 0.1
 * @date 2022-04-05
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include <string>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace std;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::Descriptors;
using namespace boost::filesystem;
using boost_ofstream = boost::filesystem::ofstream;

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		cerr << "usage: ./descriptors [INPUT SMI]" << endl;
		return 1;
	}
	const auto input_smi = argv[1];
	int pos = string(input_smi).find('.');
	path four_chem_desc = string(input_smi).substr(0, pos) + "_4properties.txt";
	path five_chem_desc = string(input_smi).substr(0, pos) + "_5properties.txt";
	boost_ofstream real4prop(four_chem_desc, ios::app);
	boost_ofstream real5prop(five_chem_desc, ios::app);
	SmilesMolSupplier smi(input_smi, "\t", 1, 0);
	size_t counter = 1;
	while (!smi.atEnd())
	{
		const unique_ptr<ROMol> mol_ptr(smi.next());
		const auto& mol = *mol_ptr;
		cout << "Processing molecule " << counter << endl;
		// generate 4 chemical properties for the molecules
		real4prop << calcExactMW(mol) << '\t'
							<< calcClogP(mol) << '\t'
							<< calcTPSA(mol) << '\t'
							<< calcLabuteASA(mol) << '\n';
		// generate 5 chemical properties for the molecules
		real5prop << mol.getNumHeavyAtoms() << '\t'
							<< calcNumHBD(mol) << '\t'
							<< calcNumHBA(mol) << '\t'
							<< calcNumRotatableBonds(mol) << '\t'
							<< calcNumRings(mol) << '\n';
		counter++;
	}
}