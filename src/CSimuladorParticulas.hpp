#ifndef CSIMULADORPARTICULAS_HPP
#define CSIMULADORPARTICULAS_HPP

#include "CGrid.hpp"
#include "CRocha.hpp"
#include "funcao.hpp"
#include "CGnuplot.hpp"
#include "metodosimpson.hpp"
#include "CParticulaFluido.hpp"

#include<vector>
#include<string>
#include<iostream>

class CSimuladorParticulas : public CParticulaFluido, CRocha {
private:
	size_t indiceTempoAtual = 0;
	size_t size_tempo;
	size_t size_malha;
	size_t numPontosIntegral=11;
	double start_x, end_x, start_t, end_t;

	std::vector<double> tempo;
	std::vector<double> malha;

	std::vector<CGrid*> resultados_ao_longo_do_tempo;

public:
	/// CONSTRUTORES
	CSimuladorParticulas(std::string pathParticulaFluido, std::string pathRocha);

	/// 'MAIN' metodo, ele que esta executando o objeto
	void run();

	/// CALCULOS
private:
	void readFile(std::string pathSimulador);
	double CalculoSigma_a(double x, double t);
	double CalculoDiffSigma_a(double x, double t);
	double CalculoLinhaZona(double x);
	double CalculoTb(double x);
	//double CalculoConcentracoes(double x, double t);
	//double CalculoConcentracoes_N_igual_1(double x, double t);
	//double CalculoConcentracoes_N_diferente_1(double x, double t);
	double CalculoDeltaPressao(double t);

	double CalculoSigma(double x, double t);
	double CalculoSigma_N_igual_1(double x, double t);
	double CalculoSigma_N_diferente_1(double x, double t);
	double funcao_sigma_n1(double x, double t);
	double funcao_sigma_n_dif_1(double x, double t);

public:
	/// metodos para SALVAR e APRESENTAR os resultados
	void printCSimuladorParticulas();
	void print_vector(std::vector<double> vetor);

	void saveInFile(std::vector<double> vector1, std::string name_vector1);
	void saveInFile(std::vector<double> vector1, std::vector<double> vector2, std::string name_vector1, std::string name_vector2);

	void plot(std::vector<double> vector1, std::vector<double> vector2, std::string name_vector1, std::string name_vector2);
};
#endif
