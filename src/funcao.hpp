#ifndef FUNCAO_HPP_
#define FUNCAO_HPP_

#include <math.h>
class Funcao1x1 {
public:
    Funcao1x1() {}
    Funcao1x1(double, double, double, double, double, double, double, double, double) {}
    Funcao1x1(double, double, double, double, double, double, double, double, double, double) {}
    virtual double operator()(double) = 0; // Funcao virtual pura 
};

class Funcao_Sigma_n_1 : public Funcao1x1 {
public:
    Funcao_Sigma_n_1(double _U, double _lambda_s, double _lambda_b, double _cb, double _C, double _sigma_a0, double _sigma_am, double _phi, double _x):
        U{ _U }, lambda_s{ _lambda_s }, lambda_b{ _lambda_b }, cb{ _cb }, C{ _C }, sigma_a0{ _sigma_a0 }, sigma_am{ _sigma_am }, phi{ _phi }, x{ _x }{}
   
    virtual double operator()(double t) {
        return (U * (lambda_s + lambda_b) * cb / (C * (sigma_a0 - sigma_am) * exp(C * phi * x / U)) - exp(-C * t)) / (1 - U * lambda_s * cb / (C * (sigma_a0 * sigma_am) * exp(C * phi * x / U)));
    }
private:
    double U, lambda_s, lambda_b, cb, C, sigma_a0, sigma_am, phi, x;
};

class Funcao_Sigma_n_diferente_1 : public Funcao1x1 {
public:
    Funcao_Sigma_n_diferente_1(double _U, double _lambda_s, double _lambda_b, double _cb, double _C, double _sigma_a0, double _sigma_am, double _phi, double _x, double _n) :
        U{ _U }, lambda_s{ _lambda_s }, lambda_b{ _lambda_b }, cb{ _cb }, C{ _C }, sigma_a0{ _sigma_a0 }, sigma_am{ _sigma_am }, phi{ _phi }, x{ _x }, n{ _n }{}

    virtual double operator()(double t) {
        double esquerda = pow(pow(sigma_a0 - sigma_am, 1 - n) - C * (1 - n) * (t - 2 * phi * x / U), n / (1 - n));
        double dir_sup = exp(x*(lambda_b+lambda_s)) * (1 - U*(lambda_s+lambda_b)*cb/ esquerda);
        double dir_inf = pow(1-lambda_s*U*cb/ esquerda, 1+lambda_b/lambda_s);
        return esquerda * (1 - dir_sup / dir_inf);
    }
private:
    double U, lambda_s, lambda_b, cb, C, sigma_a0, sigma_am, phi, x, n;
};
#endif