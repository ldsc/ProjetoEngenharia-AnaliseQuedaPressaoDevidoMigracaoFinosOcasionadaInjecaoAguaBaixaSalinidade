#include "CRocha.hpp"

CRocha::CRocha(double k, double fi, double _sigma_a0, double _sigma_am, double _N) : 
		permeabilidade{ k }, porosidade{ fi }, sigma_a0{ _sigma_a0 }, sigma_am{ _sigma_am }, N{ _N } {
	sigma_a0 = (1 - porosidade) * porcentagemCaulinita;
}

void CRocha::readFile(std::string path) {
	std::ifstream infile;
	infile.open(path);
	std::string temp;
	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// permeabilidade
	permeabilidade = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// porosidade
	porosidade = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// sigma_a0
	sigma_a0 = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// sigma_am
	sigma_am = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// N
	N = atof(temp.c_str());
}

void CRocha::printCRocha() {
	std::cout << "\n-----------------" << std::endl;
	std::cout << "Valores de CRocha" << std::endl;
	std::cout << "-----------------" << std::endl;
	std::cout << "permeabilidade: " << permeabilidade << std::endl;
	std::cout << "porosidade: " << porosidade << std::endl;
	std::cout << "sigma_a0: " << sigma_a0 << std::endl;
	std::cout << "sigma_am: " << sigma_am << std::endl;
	std::cout << "N: " << N << std::endl;
}