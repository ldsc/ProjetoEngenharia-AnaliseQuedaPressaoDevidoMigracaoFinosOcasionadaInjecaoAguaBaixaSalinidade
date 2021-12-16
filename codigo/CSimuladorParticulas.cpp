#include "CSimuladorParticulas.hpp"

CSimuladorParticulas::CSimuladorParticulas(std::string pathParticulaFluido, std::string pathRocha):CParticulaFluido(pathParticulaFluido), CRocha(pathRocha) {
	std::cout << "--------------------------------" << std::endl;
	std::cout << "----SIMULADOR DE PARTICULAS ----" << std::endl;
	std::cout << "--------------------------------" << std::endl;
	std::cout << "Deseja inserir os valores da simulacao manualmente? (s/n)" << std::endl;
	std::cout << "-->";
	std::string opt;
	std::cin >> opt;
	if (opt == "s") {
		std::cout << "Numero de tempos: ";
		std::cin >> size_tempo;
		std::cout << "Numero de espaco: ";
		std::cin >> size_malha;
		std::cout << "Posicao inicial testemunho: ";
		std::cin >> start_x;
		std::cout << "Posicao final testemunho: ";
		std::cin >> end_x;
		std::cout << "Tempo inicial: ";
		std::cin >> start_t;
		std::cout << "Tempo final: ";
		std::cin >> end_t;
	}
	else {
		std::cout <<"Infome o nome do arquivo: exemplo(Simulação.txt) ";
		std::string pathSimulador;
		std::cin >> pathSimulador;
		readFile(pathSimulador);
	}

	tempo = CGrid::linspace(start_t, end_t, size_tempo); /// nao eh objeto CGrid, eh vetor - metodo static
	malha = CGrid::linspace(start_x, end_x, size_malha);
	/// inicio as malhas
	resultados_ao_longo_do_tempo.resize(size_tempo);
	for (unsigned int i = 0; i < size_tempo; i++)
		resultados_ao_longo_do_tempo[i] = new CGrid(size_malha);
}

void CSimuladorParticulas::run() {
	int nx = resultados_ao_longo_do_tempo[0]->get_size();
	std::vector<double>linhazona(nx);
	std::vector<double>tb(nx);
	std::vector<double>delta_p(size_tempo);

	for (unsigned int j = 0; j < size_tempo; j++) { /// loop dos tempos
		for (unsigned int i = 0; i < nx; i++) { /// loop ao longo do testemunho
			resultados_ao_longo_do_tempo[j]->sigma_a[i] = CalculoSigma_a(malha[i], tempo[j]);
			resultados_ao_longo_do_tempo[j]->diff_sigma_a[i] = CalculoDiffSigma_a(malha[i], tempo[j]);
			//resultados_ao_longo_do_tempo[j]->concentracao[i] = CalculoConcentracoes(malha[i], tempo[j]);
			resultados_ao_longo_do_tempo[j]->sigma_s[i] = CalculoSigma(malha[i], tempo[j]);
			linhazona[i] = CalculoLinhaZona(malha[i]);
			tb[i] = CalculoTb(malha[i]);
		}
		delta_p[j] = CalculoDeltaPressao(tempo[j]);
		resultados_ao_longo_do_tempo[j]->saveGrid(malha, tempo[j]); /// os resultados são salvos em arquivo .txt
	}
	std::cout << "Vetor de tb: " << std::endl;
	print_vector(tb);
	saveInFile(tempo, "tempo");
	saveInFile(malha, "grid");
	plot(malha, resultados_ao_longo_do_tempo[size_tempo - 1]->get_sigma_s(), "grid", "sigma s");
	plot(tempo, delta_p, "tempo", "dp");
}

/// abaixo estao os calculos do simulador

double CSimuladorParticulas::CalculoSigma_a(double x, double t){
	if (N == 1) {
		if (t < (x * porosidade / velocidade))
			return sigma_a0;
		else {
			return sigma_am + (sigma_a0 - sigma_am) * exp(-parC * (t - x * porosidade / velocidade));
		}
	}	
	else {
		if (t < (x * porosidade / velocidade))
			return sigma_a0;
		else {
			return pow(sigma_am+pow(sigma_a0-sigma_am, 1- N) - parC *(1- N)*(t-x* porosidade / velocidade), N / (1- N));
		}
	}
}

double CSimuladorParticulas::CalculoDiffSigma_a(double x, double t) {
	if (N == 1) {
		if (t < (x * porosidade / velocidade))
			return 0.0;
		else {
			return -parC *(sigma_a0 - sigma_am) * exp(-parC * (t - x * porosidade / velocidade));
		}
	}
	else {
		if (t < (x * porosidade / velocidade))
			return 0.0;
		else {
			return -parC *pow( pow(sigma_a0 - sigma_am, 1 - N) - parC * (1 - N) * (t - x * porosidade / velocidade), N / (1 - N));
		}
	}
}

double CSimuladorParticulas::CalculoLinhaZona(double x) {
	return porosidade * x / velocidade;
}

double CSimuladorParticulas::CalculoTb(double x) {
	double cb = 0.0004; /// rocha
	double linhaZona = CalculoLinhaZona(x);

	double tb;
	if (N == 1)
		tb = linhaZona + (1/ parC) * log(parC *(sigma_a0 - sigma_am)*(1-exp(-lambdaPonte *x))/(cb* velocidade * lambdaPonte));
	else
		tb = 2 * linhaZona + pow(sigma_a0 - sigma_am, 1 - N) / (parC * (1 - N)) - (1 / (parC * (1 - N))) * pow(lambdaPonte * velocidade * cb / (parC * (1 - exp(-x * lambdaPonte))), (1 - N) / N);

	return tb;
}

//double CSimuladorParticulas::CalculoConcentracoes(double x, double t) {
//	if (N == 1)
//		return CalculoConcentracoes_N_igual_1(x, t);
//
//	else
//		return CalculoConcentracoes_N_diferente_1(x, t);
//}

//double CSimuladorParticulas::CalculoConcentracoes_N_igual_1(double x, double t) {
//	double concentracao;
//	if (t < CalculoTb(x))
//		concentracao = 0.0;
//	else {
//		if (C < cb)
//			concentracao = C * (sigma_a0 - sigma_am) * exp(-C * (t - CalculoTb(x))) * (1 - exp(-lambdaAdesao * x)) / (velocidade * lambdaAdesao);
//		else
//			concentracao = C * (sigma_a0 - sigma_am) * exp(-C * t - CalculoTb(x)) / (velocidade * (lambdaAdesao + lambdaPonte)) * (1 - (1 - velocidade * (lambdaAdesao + lambdaPonte) * cb / (C * (sigma_a0 - sigma_am) * exp(-C * (t - CalculoTb(x))))) / pow(1 - velocidade * lambdaAdesao * cb / (C * (sigma_a0 - sigma_am) * exp(-C * (t - CalculoTb(x)))), 1 + lambdaPonte / lambdaAdesao) * exp((lambdaAdesao + lambdaPonte) * x));
//	}
//	return concentracao;
//}

//double CSimuladorParticulas::CalculoConcentracoes_N_diferente_1(double x, double t) {
//	double concentracao;
//	if (t < CalculoTb(x))
//		concentracao = 0.0;
//	else {
//		if (C < cb)
//			concentracao = C * pow(pow(sigma_a0 - sigma_am, 1 - N) - C * (1 - N) * (t - 2 * CalculoTb(x)), N / (N - 1)) * (1 - exp(-x * lambdaAdesao)) / (lambdaAdesao * velocidade);
//		else {
//			double c1 = C * pow(pow(sigma_a0 - sigma_am, 1 - N) - C * (1 - N) * (t - 2 * CalculoTb(x)), N / (N - 1));
//			double esquerda = c1 / (velocidade * (lambdaAdesao + lambdaPonte));
//			double direita = 1-(1 - (velocidade*(lambdaAdesao+lambdaPonte))/c1)*exp(-x*(lambdaAdesao+lambdaPonte))/ pow(1-lambdaAdesao*velocidade*cb/c1,1+lambdaPonte/lambdaAdesao);
//			concentracao = esquerda * direita;
//		}
//	}
//	return concentracao;
//}

double CSimuladorParticulas::CalculoDeltaPressao(double t) {
	double integral = 0.0;
	double const1 = viscosidade * velocidade * end_x / permeabilidade;
	for (int i = 0; i < size_malha; i++) {
		integral += CalculoSigma(malha[i], t)*(malha[1]-malha[0]);
	}
	return -1*(const1 + const1 * beta * integral);
}

double CSimuladorParticulas::CalculoSigma(double x, double t) {
	if (N == 1)
		return CalculoSigma_N_igual_1(x, t);
	else
		return CalculoSigma_N_diferente_1(x, t);
}

double CSimuladorParticulas::CalculoSigma_N_igual_1(double x, double t) {
	double sigma = 0.0;
	double tb = CalculoTb(x);
	double C0 = x * porosidade / velocidade;
	if (t < C0)
		sigma = 0.0;
	else if (t >= C0 && t < tb) {
		sigma = funcao_sigma_n1(x, t);
	}
	else {
		double maxT = tb > C0 ? tb : C0;
		sigma = funcao_sigma_n1(x, maxT) + (sigma_a0 - sigma_am) * (1 - exp(-lambdaAdesao * x)) * exp(C * porosidade * x / velocidade) * (exp(-C * maxT) - exp(-C * t));
	}
	return sigma;
}

double CSimuladorParticulas::CalculoSigma_N_diferente_1(double x, double t) {
	double sigma = 0.0;
	double tb = CalculoTb(x);
	double C0 = x * porosidade / velocidade;
	if (t < C0)
		sigma = 0.0;
	else if (t >= C0 && t < tb) {
		sigma = funcao_sigma_n1(x, t);
	}
	else {
		double maxT = tb > C0 ? tb : C0;
		sigma = funcao_sigma_n1(x, maxT) +(1-exp(-x*lambdaAdesao))*(N-1)/(1-N)*pow(pow(sigma_a0-sigma_am, 1-N) - C*(1-N)*(t-2*porosidade*x/velocidade), 1/(1-N));
	}
	return sigma;
}

double CSimuladorParticulas::funcao_sigma_n1(double x, double t) {
	Funcao_Sigma_n_1 funcao(velocidade, lambdaAdesao, lambdaPonte, cb, C, sigma_a0, sigma_am, porosidade, x);
	MetodoIntegracaoNumerica1D* metodo = new MetodoSimpson(funcao);
	double integral = metodo->Integrar(porosidade * x / velocidade, t, numPontosIntegral);
	return C * (sigma_a0 - sigma_am) * ((1 - exp(-C * (t - porosidade * x / velocidade)) / C) + exp(-(lambdaAdesao + lambdaPonte - C * (x * porosidade / velocidade) * x)) * integral);
}

double CSimuladorParticulas::funcao_sigma_n_dif_1(double x, double t) {
	/// na linha abaixo, é criado a funcao relacionado a funcao sigma
	Funcao_Sigma_n_diferente_1 funcao(velocidade, lambdaAdesao, lambdaPonte, cb, C, sigma_a0, sigma_am, porosidade, x, N);

	/// na linha abaixo, é criado o método de simpson, e é enviada a função criada acima
	MetodoIntegracaoNumerica1D* metodo = new MetodoSimpson(funcao);
	// na linha abaixo, é executado o método para integrar
	double integral = metodo->Integrar(porosidade * x / velocidade, t, numPontosIntegral);
	return C * integral;
}

void CSimuladorParticulas::printCSimuladorParticulas() {
	printCParticulaFluido();
	printCRocha();
	std::cout << "\n--------------------" << std::endl;
	std::cout << "Classe da simulacao: " << std::endl;
	std::cout << "--------------------" << std::endl;

	std::cout << "Grid dos tempos:" << std::endl;
	print_vector(tempo);
}

void CSimuladorParticulas::saveInFile(std::vector<double> vector1, std::string name_vector1) {
	std::ofstream outdata; //save data
	outdata.open((name_vector1 + ".dat").c_str());
	outdata << "# "<< name_vector1 << std::endl;
	for (unsigned int i = 0; i < vector1.size(); i++)
		outdata << vector1[i] << std::endl;
	outdata.close();
}

void CSimuladorParticulas::saveInFile(std::vector<double> vector1, std::vector<double> vector2, std::string name_vector1, std::string name_vector2){
	if (vector1.size() != vector2.size()) {
		std::cout << "Nao foi possivel salvar os vetores, por terem tamanhos distintos!"<<std::endl;
		return;
	}

	std::ofstream outdata; //save data
	outdata.open((name_vector1+"_"+ name_vector2 + ".dat").c_str());
	outdata << "# " << name_vector1 << " " << name_vector2 << std::endl;
	for (unsigned int i = 0; i < vector1.size(); i++)
		outdata << vector1[i] << " " << vector2[i] << std::endl;
	outdata.close();
}

void CSimuladorParticulas::print_vector(std::vector<double> vetor) {
	std::cout << vetor[0]; /// este primeiro nao fica dentro do loop por causa do ' - '
	for (unsigned int i = 1; i < vetor.size(); i++)
		std::cout << " - " << vetor[i];
	std::cout << std::endl;
}

void CSimuladorParticulas::plot(std::vector<double> vector1, std::vector<double> vector2, std::string name_vector1, std::string name_vector2) {
	saveInFile(vector1, vector2, name_vector1, name_vector2);
	CGnuplot::plot((name_vector1 + "_" + name_vector2 + ".dat").c_str(), name_vector1, name_vector2, (name_vector1 + "_" + name_vector2 + ".png").c_str());
}

void CSimuladorParticulas::readFile(std::string pathSimulador) {

	std::ifstream infile;
	infile.open(pathSimulador);
	std::string temp;
	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// size tempo
	size_tempo = atoi(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// size malha
	size_malha = atoi(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// num pontos integral
	numPontosIntegral = atoi(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// start_x
	start_x = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// end_x
	end_x = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// start_t
	start_t = atof(temp.c_str());

	std::getline(infile, temp); // pular linha com texto
	std::getline(infile, temp);		// end_t
	end_t = atof(temp.c_str());
}
