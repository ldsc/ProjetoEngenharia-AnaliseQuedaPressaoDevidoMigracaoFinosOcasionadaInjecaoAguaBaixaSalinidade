#include "metodosimpson.hpp"

double MetodoSimpson::Integrar(double a, double b, size_t numPontos) {
    // se numPontos for par, somo 1
    if (numPontos % 2 == 0) {
        std::cout << "\nWarning: numPontos passado eh par, somando 1 para integrar pelo metodo de Simpson!"<<std::endl;
        numPontos += 1;
    }

    double h = (b-a)/(numPontos-1);

    // método de Simpson
    double dx = (b - a) / (numPontos - 1);

    double resultado = (funcao(a) + funcao(b)) * h / 3;

    for (size_t i = 1; i < numPontos - 1; i++) {
        if (i % 2 == 0) // se i for par, multiplico por 2
            resultado += funcao(a + i * dx) * 2 * h / 3;
        else // se if for impar, multiplico por 4
            resultado += funcao(a + i * dx) * 4 * h / 3;
    }

    return resultado;
}
