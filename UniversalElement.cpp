//
// Created by Adrian on 15.11.2018.
//

#include "UniversalElement.h"
#include<iostream>
#include<math.h>
#include"BasicStructures.h"
using namespace std;

void UniversalElement::calculate(node * nodes[4], double **globalH, double **globalC, double *globalP) {
    //cout<<"Beginning of universal element calculation"<<endl;

    //Setting integral points
    integralPoints[0].x=nodes[0]->x;
    integralPoints[0].y=nodes[0]->y;
    integralPoints[1].x=nodes[1]->x;
    integralPoints[1].y=nodes[1]->y;
    integralPoints[2].x=nodes[2]->x;
    integralPoints[2].y=nodes[2]->y;
    integralPoints[3].x=nodes[3]->x;
    integralPoints[3].y=nodes[3]->y;


    //Setting ksi and eta
    double almostSqrt3 = 1/sqrt(3);
    ksi[0]=(-1)*almostSqrt3;
    ksi[1]=almostSqrt3;
    ksi[2]=almostSqrt3;
    ksi[3]=(-1)*almostSqrt3;
    eta[0]=(-1)*almostSqrt3;
    eta[1]=(-1)*almostSqrt3;
    eta[2]=almostSqrt3;
    eta[3]=almostSqrt3;

    //Setting N's - functions of shape
    for (int i=0; i<4; i++)
    {
        N[0][i]=0.25*(1-ksi[i])*(1-eta[i]);
        N[1][i]=0.25*(1+ksi[i])*(1-eta[i]);
        N[2][i]=0.25*(1+ksi[i])*(1+eta[i]);
        N[3][i]=0.25*(1-ksi[i])*(1+eta[i]);
    }


    //Setting dN/dKsi and dN/dEta
    // i - number of shape function
    // j - number of integral point

    for (int j=0; j<4; j++)
    {
        dN_dKsi[0][j]=-0.25*(1-eta[j]);
        dN_dKsi[1][j]=0.25*(1-eta[j]);
        dN_dKsi[2][j]=0.25*(1+eta[j]);
        dN_dKsi[3][j]=-0.25*(1+eta[j]);
        dN_dEta[0][j]=-0.25*(1-ksi[j]);
        dN_dEta[1][j]=-0.25*(1+ksi[j]);
        dN_dEta[2][j]=0.25*(1+ksi[j]);
        dN_dEta[3][j]=0.25*(1-ksi[j]);
    }

    //Setting Jacobian, Jacobian divided by its determinant
    // i - number of integral point
    for (int i=0; i<4; i++)
    {
        //  dx/dKsi
        jacobian[0][0][i]=0;
        //  dy/dKsi
        jacobian[0][1][i]=0;
        //  dx/dEta
        jacobian[1][0][i]=0;
        //  dy/dEta
        jacobian[1][1][i]=0;
        for (int j=0; j<4; j++)
        {
            jacobian[0][0][i] += dN_dKsi[j][i]*integralPoints[j].x;
            jacobian[0][1][i] += dN_dKsi[j][i]*integralPoints[j].y;
            jacobian[1][0][i] += dN_dEta[j][i]*integralPoints[j].x;
            jacobian[1][1][i] += dN_dEta[j][i]*integralPoints[j].y;
        }
    }


    /*
    for (int i=0; i<4; i++) cout<<jacobian[0][0][i]<<"\t";
    cout<<endl;
    for (int i=0; i<4; i++) cout<<jacobian[0][1][i]<<"\t";
    cout<<endl;
    for (int i=0; i<4; i++) cout<<jacobian[1][0][i]<<"\t";
    cout<<endl;
    for (int i=0; i<4; i++) cout<<jacobian[1][1][i]<<"\t";
    cout<<endl<<endl;
     */

    // i - number of integral point
    for (int i=0; i<4; i++)
    {
        detJ[i]=jacobian[0][0][i]*jacobian[1][1][i]-jacobian[0][1][i]*jacobian[1][0][i];
        //cout<<"det["<<i<<"] = "<<detJ[i]<<endl;
        jacobianDetJ[0][0][i]=jacobian[1][1][i]/detJ[i];
        jacobianDetJ[0][1][i]=-jacobian[0][1][i]/detJ[i];
        jacobianDetJ[1][0][i]=-jacobian[1][0][i]/detJ[i];
        jacobianDetJ[1][1][i]=jacobian[0][0][i]/detJ[i];
    }

    //Setting dN/dx
    // i - number of integral point
    // j - number of shape function
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            dN_dx[i][j] = jacobianDetJ[0][0][i] * dN_dKsi[j][i] + jacobianDetJ[0][1][i] * dN_dEta[j][i];
            dN_dy[i][j] = jacobianDetJ[0][1][i] * dN_dKsi[j][i] + jacobianDetJ[1][1][i] * dN_dEta[j][i];
        }
    }

    //Multiplying for example {dN/dx}{dN/dx}T
    // i - number of integral point
    for (int i=0; i<4; i++)
    {
        for (int a=0; a<4; a++)
        {
            for (int b=0; b<4; b++)
            {
                dN_dx_dN_dx[a][b][i]=dN_dx[i][b]*dN_dx[i][a]*detJ[i];
                dN_dy_dN_dy[a][b][i]=dN_dy[i][b]*dN_dy[i][a]*detJ[i];
            }
        }
    }

    //Calculating "K*(     {dN/dx}{dN/dx}T  +  {dN/dy}{dN/dy}T)*DetJ
    // i - number of integral point
    for (int i=0; i<4; i++)
    {
        for (int a=0; a<4; a++)
        {
            for (int b=0; b<4; b++)
            {
                k_dN[a][b][i]=k*(dN_dx_dN_dx[a][b][i]+dN_dy_dN_dy[a][b][i]);
            }
        }
    }

    //Calculating final Matrix H
    for (int a=0; a<4; a++)
    {
        for (int b=0; b<4; b++)
        {
            matrixH[a][b]=0;
            for (int i=0; i<4; i++) matrixH[a][b] += k_dN[a][b][i];
            //cout<<matrixH[a][b]<<"\t";
        }
        //cout<<endl;
    }

    //Calculating Matrix C
    // i - number of integral point
    for (int a=0; a<4; a++)
    {
        for (int b=0; b<4; b++) matrixC[a][b]=0;
    }
    for (int i=0; i<4; i++)
    {
        for (int a=0; a<4; a++)
        {
            for (int b=0; b<4; b++)
            {
                //matrixC[a][b]+=detJ[i]*c*ro*N[i][a]*N[i][b];
                matrixC[a][b]+=detJ[i]*c*ro*N[a][i]*N[b][i];
            }
        }
    }

    //Setting ksi and eta tables for Matrix H Boundary Conditions
    // i - number of surface       j - number of integral point
    ksiBC[0][0]=-almostSqrt3; ksiBC[0][1]=almostSqrt3; etaBC[0][0]=-1; etaBC[0][1]=-1;
    ksiBC[1][0]=1; ksiBC[1][1]=1; etaBC[1][0]=-almostSqrt3; etaBC[1][1]=almostSqrt3;
    ksiBC[2][0]=almostSqrt3; ksiBC[2][1]=-almostSqrt3; etaBC[2][0]=1; etaBC[2][1]=1;
    ksiBC[3][0]=-1; ksiBC[3][1]=-1; etaBC[3][0]=almostSqrt3; etaBC[3][1]=-almostSqrt3;

    //Calculating N's for BC
    // i - number of surface    a - number of integral point    b - number of shape function
    for (int i=0; i<4; i++)
    {
        for (int a=0; a<2; a++)
        {
            NBC[i][a][0]=0.25*(1-ksiBC[i][a])*(1-etaBC[i][a]);
            NBC[i][a][1]=0.25*(1+ksiBC[i][a])*(1-etaBC[i][a]);
            NBC[i][a][2]=0.25*(1+ksiBC[i][a])*(1+etaBC[i][a]);
            NBC[i][a][3]=0.25*(1-ksiBC[i][a])*(1+etaBC[i][a]);
        }
    }

    //Setting side lengths
    sideLength[0]=integralPoints[1].x-integralPoints[0].x;
    sideLength[1]=integralPoints[2].y-integralPoints[1].y;
    sideLength[2]=integralPoints[2].x-integralPoints[3].x;
    sideLength[3]=integralPoints[3].y-integralPoints[0].y;

    //Calculating Matrix H BC Partial
    // i - number of surface    j - number of integral point
    for (int i=0; i<4; i++)
    {
        for (int a=0; a<4; a++)
        {
            for (int b=0; b<4; b++)
            {
                matrixHBCPartial[i][0][a][b]=alfa*NBC[i][0][a]*NBC[i][0][b];
                matrixHBCPartial[i][1][a][b]=alfa*NBC[i][1][a]*NBC[i][1][b];
                matrixHBCPartial[i][2][a][b]=(matrixHBCPartial[i][0][a][b]+matrixHBCPartial[i][1][a][b])*(sideLength[i]/2);
            }
            vectorPPartial[i][0][a]=alfa*ambientT*NBC[i][0][a];
            vectorPPartial[i][1][a]=alfa*ambientT*NBC[i][1][a];
            vectorPPartial[i][2][a]=(vectorPPartial[i][0][a]+vectorPPartial[i][1][a])*(sideLength[i]/2);
        }
    }
    /*
    for (int i=0; i<4; i++) {
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                cout << matrixHBCPartial[i][2][a][b] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    */

    //Calculating MatrixHBC
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            matrixHBC[i][j]=0;
        }
    }
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            if (nodes[0]->BC and nodes[1]->BC) matrixHBC[i][j]+=matrixHBCPartial[0][2][i][j];
            if (nodes[1]->BC and nodes[2]->BC) matrixHBC[i][j]+=matrixHBCPartial[1][2][i][j];
            if (nodes[2]->BC and nodes[3]->BC) matrixHBC[i][j]+=matrixHBCPartial[2][2][i][j];
            if (nodes[0]->BC and nodes[3]->BC) matrixHBC[i][j]+=matrixHBCPartial[3][2][i][j];
        }
    }
    for (int i=0; i<4; i++) vectorP[i]=0;

    for (int i=0; i<4; i++) {
        if (nodes[0]->BC and nodes[1]->BC) vectorP[0] += vectorPPartial[i][2][0];
        if (nodes[1]->BC and nodes[2]->BC) vectorP[1] += vectorPPartial[i][2][1];
        if (nodes[2]->BC and nodes[3]->BC) vectorP[2] += vectorPPartial[i][2][2];
        if (nodes[0]->BC and nodes[3]->BC) vectorP[3] += vectorPPartial[i][2][3];
    }
    for (int i=0; i<4; i++)
    {
        int ii=nodes[i]->ID;
        globalP[ii]+=vectorP[i];
    }

    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            matrixH[i][j]+=matrixHBC[i][j];
            int ii=nodes[i]->ID;
            int jj=nodes[j]->ID;
            globalH[ii][jj]+=matrixH[i][j];
            globalC[ii][jj]+=matrixC[i][j];

        }
    }

    //cout<<"End of universal element calculation. "<<endl;
}

UniversalElement::UniversalElement() {
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            matrixH[i][j]=0;
        }
    }
}

double **UniversalElement::getMatrixH() {
    double** a = new double*[4];
    for(int i = 0; i < 4; i++)
        a[i] = new double[4];
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
            a[i][j]=matrixH[i][j];
    }
    return a;
}

UniversalElement::UniversalElement(double alfa, double cw, double k, double ro, double ambientT) {
    this->k=k;
    this->alfa=alfa;
    this->c=cw;
    this->ro=ro;
    this->ambientT=ambientT;
}
