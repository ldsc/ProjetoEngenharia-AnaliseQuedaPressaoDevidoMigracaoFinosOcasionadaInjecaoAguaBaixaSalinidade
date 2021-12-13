#include "CParticulaFluido.hpp"

void CParticulaFluido::readFile(std::string path) {
	std::ifstream infile;
	infile.open(path);
	std::string temp;
	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// viscosidade
	viscosidade = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// lambdaPonte
	lambdaPonte = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// lambdaAdesao
	lambdaAdesao = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// parC
	parC = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// par_n
	par_n = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// velocidade
	velocidade = atof(temp.c_str());

}

void CParticulaFluido::printCParticulaFluido() {
	std::cout << "\n--------------------------" << std::endl;
	std::cout << "Valores de ParticulaFluido" << std::endl;
	std::cout << "--------------------------" << std::endl;
	std::cout << "viscosidade: " << viscosidade << std::endl;
	std::cout << "lambdaPonte: " << lambdaPonte << std::endl;
	std::cout << "lambdaAdesao: " << lambdaAdesao << std::endl;
	std::cout << "parC: " << parC << std::endl;
	std::cout << "par_n: " << par_n << std::endl;
	std::cout << "velocidade: " << velocidade << std::endl;
}