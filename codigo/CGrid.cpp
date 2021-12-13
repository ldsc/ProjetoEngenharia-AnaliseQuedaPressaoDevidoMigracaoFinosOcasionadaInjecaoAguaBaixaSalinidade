#include "CGrid.hpp"

void CGrid::saveGrid(std::vector<double> malha, double time) {
	std::string path = "resultados_malha//malha_" + std::to_string(time) + ".txt";

	std::ofstream file(path);
	file << "     malha  -  sigma_a  -  sigma_s - diff_sigma_a\n";
	for (unsigned int i = 0; i < _size; i++) {
		file << std::setw(9) << malha[i]
			//<< std::setw(12) << c_susp[i]
			<< std::setw(12) << sigma_a[i]
			<< std::setw(12) << sigma_s[i]
			<< std::setw(15) << diff_sigma_a[i]
			//<< std::setw(15) << concentracao[i] << "\n";
	}
	file.close();
}

std::vector<double> CGrid::linspace(double start, double end, int size) {
	double delta = (end - start) / (1.0 * size - 1);
	std::vector<double> grid;

	for (int i = 0; i < size; i++)
		grid.push_back(start + delta * i);

	return grid;
}
