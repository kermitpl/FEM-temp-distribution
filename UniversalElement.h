//
// Created by Adrian on 15.11.2018.
//

#ifndef MES_MATRIXH_H
#define MES_MATRIXH_H

#include "BasicStructures.h"

class UniversalElement {
private:
    double matrixH[4][4];
    point integralPoints[4];
    double ksi[4], eta[4];
    double N[4][4];
    double dN_dKsi[4][4], dN_dEta[4][4];
    double jacobian[2][2][4], jacobianDetJ[2][2][4];
    double detJ[4];
    double k;
    double c;
    double ro;
    double alfa;
    double ambientT;
    double dN_dx[4][4], dN_dy[4][4];
    double dN_dx_dN_dx[4][4][4], dN_dy_dN_dy[4][4][4];
    double k_dN[4][4][4];
    double matrixC[4][4];
    double matrixHBCPartial[4][3][4][4];
    double vectorPPartial[4][3][4];
    double ksiBC[4][2], etaBC[4][2];
    double NBC[4][2][4];
    double sideLength[4];
    double matrixHBC[4][4];
    double vectorP[4];
public:
    void calculate(node * nodes[4], double **globalH, double **globalC, double *globalP);
    double ** getMatrixH();
    UniversalElement();
    UniversalElement(double alfa, double cw, double k, double ro, double ambientT);
};


#endif //MES_MATRIXH_H
