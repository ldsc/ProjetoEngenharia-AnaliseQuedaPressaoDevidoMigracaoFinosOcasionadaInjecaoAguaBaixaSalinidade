#ifndef CGRID_HPP
#define CGRID_HPP

#include <string>
#include <vector>
#include <fstream> /// escrecrever em disco
#include <iomanip>     /// setw

class CGrid {
private:
	int _size;
public:
	std::vector<double> c_susp, sigma_a, sigma_s, diff_sigma_a, concentracao;

public:
	CGrid(int size_x) {
		_size = size_x;
		c_susp.resize(_size, 0.0);
		sigma_a.resize(_size, 0.0);
		sigma_s.resize(_size,0.0);
		diff_sigma_a.resize(_size,0.0);
		concentracao.resize(_size, 0.0);
	}

	int get_size() { return _size;}

	void saveGrid(std::vector<double> malha, double time);
	std::vector<double> get_sigma_s() { return sigma_s; }
	/// metodo para criar malhas, eh static para criar as malhas do espaco e do tempo
	std::vector<double> static linspace(double start, double end, int size);
};
#endif