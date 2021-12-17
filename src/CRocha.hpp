#ifndef CROCHA_HPP
#define CROCHA_HPP

#include <fstream>
#include <iostream>
#include <string>

class CRocha{
protected:
	double permeabilidade;
	double porosidade;
	double sigma_am;
	double N;
	double beta = 0.5;
	double porcentagemCaulinita = 0.2;

	double sigma_a0;

	void readFile(std::string path);
public:
	//CRocha() {}
	// construtores por variáveis e por nome do arquivo txt
	CRocha(double k, double fi, double _sigma_a0, double _sigma_am, double _N);// : permeabilidade{ k }, porosidade{ fi }, sigma_a0{ _sigma_a0 }, sigma_am{ _sigma_am }, N{ _N } {}
	CRocha(std::string path) { 
		readFile(path); }

	double reducaoPermeabilidade(double sigma) { return 1 / (1 + beta * sigma); }

	// função para mostrar os parâmetros da classe
	void printCRocha();
};

#endif
