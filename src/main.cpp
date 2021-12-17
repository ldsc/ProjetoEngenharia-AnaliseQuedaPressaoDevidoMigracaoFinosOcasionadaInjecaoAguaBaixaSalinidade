#include<iostream>
#include<string>

#include "CSimuladorParticulas.hpp"

using namespace std;

int main() {

	/////////////////////////////////////////////////////////
	CSimuladorParticulas simulacao("particulaFluido.txt", "rocha.txt");
	simulacao.printCSimuladorParticulas();
	simulacao.run();
	/////////////////////////////////////////////////////////
	/*delete x_par_ptr;
	delete y_par_ptr;
	delete objSeno2Dptr;
	delete metodo;*/
}
