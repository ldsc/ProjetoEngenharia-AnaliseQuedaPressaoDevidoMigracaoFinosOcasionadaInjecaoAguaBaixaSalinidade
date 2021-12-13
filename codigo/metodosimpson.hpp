#ifndef METODOTRAPEZIO_HPP_
#define METODOTRAPEZIO_HPP_

#include "metodointegracaonumerica.hpp"
#include<iostream>

class MetodoSimpson : public MetodoIntegracaoNumerica1D {
public:
    // construtor
    MetodoSimpson(Funcao1x1& f_x_in) : MetodoIntegracaoNumerica1D(f_x_in) {}

    // a função abaixo resolve a integral pelo método de trapezio... a função é usada aqui dentro
    virtual double Integrar(double a, double b, size_t numPontos);
};
#endif