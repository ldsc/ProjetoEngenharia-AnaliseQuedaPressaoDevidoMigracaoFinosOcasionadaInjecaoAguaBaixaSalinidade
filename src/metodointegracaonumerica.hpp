#ifndef METODOINTEGRACAONUMERICA_HPP_
#define METODOINTEGRACAONUMERICA_HPP_

#include <cstddef>

#include "funcao.hpp"

class MetodoIntegracaoNumerica1D {
public:
    MetodoIntegracaoNumerica1D(Funcao1x1& fx_in) : funcao(fx_in) {}
    virtual double Integrar(double a, double b, size_t numPontos) = 0; // Virtual pura

    Funcao1x1& funcao; // Funcao 1D a ser integrada
};
#endif