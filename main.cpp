#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>
#include "UniversalElement.h"
#include "BasicStructures.h"

using namespace std;

void showNodesGrid(node myNodes[], int nH, int nL);
void showNodesTemps(double *temps, int nH, int nL, int time);
void getDataFromFile(double &H, double &L, int &nH, int &nL, double &t0, int &time, int &step, double &ambientT, double &alfa, double &cw, double &k, double &ro);
void nodesCreation(node * myNodes, double H, double L, int nH, int nL, int nN, double t0);
void gridCreation(grid myGrid, node * myNodes, int nE, int nH, int nL);
double *gauss(double **A, int n);

int main() {
    //Settings of precision
    cout<<setprecision(3);
    cout<<fixed;

    double H, L, t0, alfa, cw, k, ro, ambientT;
    int nH, nL, time, step;

    //Data from file
    getDataFromFile(H, L, nH, nL, t0, time, step, ambientT, alfa, cw, k, ro);

    //Show data
    int nN = nH*nL, nE = (nH-1)*(nL-1);
    cout<<endl<<"----------------------------------------------------------------------------------------------"<<endl;
    cout << "H = "<<H<<"\t L = "<<L<<"\t nH = "<<nH<<"\t nL = "<<nL<<"\t t0 = "<<t0<<"\t nN = "<<nN<<"\t nE = "<<nE<<endl;
    cout << "Time = "<<time<<"\t step = "<<step<<"\t ambientT = "<<ambientT<<"\t alfa = "<<alfa<<"\t cw = "<<cw<<"\t k = "<<k<<"\t ro = "<<ro<<endl;
    cout<<"----------------------------------------------------------------------------------------------"<<endl<<endl;

    //Nodes creation
    node * myNodes = new node [nN];
    nodesCreation(myNodes, H, L, nH, nL, nN, t0);

    //Showing nodes grid
    showNodesGrid(myNodes, nH, nL);

    //Grid creation
    grid myGrid{};
    myGrid.elements = new element [nE];
    gridCreation(myGrid, myNodes, nE, nH, nL);

    //Showing grid
    //myGrid.show(nE);

    //Matrix H and C and Vector P - preparing for calculation
    double **globalH = new double * [nN];
    for (int i=0; i<nN; i++) globalH[i] = new double [nN];
    double **globalC = new double * [nN];
    for (int i=0; i<nN; i++) globalC[i] = new double [nN];
    for (int i=0; i<nN; i++)
    {
        for (int j=0; j<nN; j++) {
            globalH[i][j] = 0;
            globalC[i][j] = 0;
        }
    }
    double *globalP = new double [nN];
    for (int i=0; i<nN; i++) globalP[i]=0;
    double *globalPUnchanged = new double [nN];
    for (int i=0; i<nN; i++) globalPUnchanged[i]=0;

    //Calculation of global matrixes and vectors
    UniversalElement universalElement(alfa, cw, k, ro, ambientT);
    for (int i=0; i<nE; i++) {
        universalElement.calculate(myGrid.elements[i].eNodes, globalH, globalC, globalP);
    }

    //Preparing for loops with Gaussian Elimination
    double rest;
    double **gauss1 = new double * [nN];
    for (int i=0; i<nN; i++) gauss1[i] = new double [nN+1];

    //First loop
    for (int i=0; i<nN; i++)
    {
        rest=0;
        for (int j=0; j<nN; j++) {
            globalH[i][j]+=globalC[i][j]/step;
            rest +=(globalC[i][j]*t0)/step;
        }
        globalPUnchanged[i]=globalP[i];
        globalP[i]+=rest;
    }

    //Gauss
    for (int i=0; i<nN; i++){
        for (int j=0; j<nN; j++)
        {
            gauss1[i][j]=globalH[i][j];
        }
    }
    for (int i=0; i<nN; i++) gauss1[i][nN]=globalP[i];

    double *resultTemp=gauss(gauss1, nN);
    showNodesTemps(resultTemp, nH, nL, step);

    //Next steps
    for (int s=step; s<time; s+=step)
    {
        //Preparation
        for (int i=0; i<nN; i++)
        {
            globalP[i]=globalPUnchanged[i];
            rest=0;
            for (int j=0; j<nN; j++) {
                rest +=(globalC[i][j]*resultTemp[j])/step;
            }
            globalP[i]+=rest;
        }

        //Setting new matrix for gauss
        for (int i=0; i<nN; i++){
            for (int j=0; j<nN; j++)
            {
                gauss1[i][j]=globalH[i][j];
            }
        }
        for (int i=0; i<nN; i++) gauss1[i][nN]=globalP[i];

        resultTemp=gauss(gauss1, nN);
        showNodesTemps(resultTemp, nH, nL, s+step);
    }

    //Freeing memory
    delete [] gauss1;
    delete [] globalH;
    delete [] globalC;
    delete [] myNodes;
    delete [] myGrid.elements;

    cout<<endl<< "End of programme."<<endl;
    return 0;
}

void showNodesGrid(node myNodes[], int nH, int nL)
{
    for (int i=nH-1; i>=0; i--)
    {
        for (int j=0; j<=(nL-1); j++)
        {
            cout<<"{"<<i+j*nH<<"} ";
            myNodes[i+j*nH].show();
        }
        cout<<endl;
    }
}

void showNodesTemps(double *temps, int nH, int nL, int time)
{
    cout<<endl<<"Time: "<<time<<"  "<<"\t Max: "<<*max_element(temps, temps+nH*nL)<<"\t Min: "<<*min_element(temps, temps+nH*nL)<<endl;
    for (int i=nH-1; i>=0; i--)
    {
        for (int j=0; j<=(nL-1); j++)
        {
            cout<<temps[i+j*nH]<<"\t";
        }
        cout<<endl;
    }
    cout<<endl;
}


void getDataFromFile(double &H, double &L, int &nH, int &nL, double &t0, int &time, int &step, double &ambientT, double &alfa, double &cw, double &k, double &ro)
{
    fstream file;
    string path = "dane.txt";
    file.open(path, ios::in);
    if (file.good())
    {
        try
        {
            string temp;
            cout<<"Access to file granted."<<endl;
            file>>temp>>temp>>H;
            file>>temp>>temp>>L;
            file>>temp>>temp>>nH;
            file>>temp>>temp>>nL;
            file>>temp>>temp>>t0;
            file>>temp>>temp>>time;
            file>>temp>>temp>>step;
            file>>temp>>temp>>ambientT;
            file>>temp>>temp>>alfa;
            file>>temp>>temp>>cw;
            file>>temp>>temp>>k;
            file>>temp>>temp>>ro;
        }
        catch( exception & e)
        {
            cout<<"File corrupted. Using default data."<<endl;
        }
    }
    else cout<<"Access to file denied."<<endl;
}

void nodesCreation(node * myNodes, double H, double L, int nH, int nL, int nN, double t0)
{
    int currentRow=0;
    int currentColumn=0;
    for (int i=0; i<nN; i++)
    {
        if (currentColumn==0 or currentRow==0 or currentColumn==(nL-1) or currentRow==(nH-1)) myNodes[i].BC=true;
        myNodes[i].ID=i;
        myNodes[i].t0=t0;
        myNodes[i].x=0+currentColumn*(L/(nL-1));
        myNodes[i].y=0+currentRow*(H/(nH-1));
        currentRow++;
        if (currentRow%(nH)==0)
        {
            currentColumn++;
            currentRow=0;
        }
    }

}

void gridCreation(grid myGrid, node * myNodes, int nE, int nH, int nL)
{
    int j=0;
    int licznik=0;
    for (int i=0; i<nE; i++)
    {
        myGrid.elements[i].eNodes[0]=&myNodes[j];
        myGrid.elements[i].eNodes[1]=&myNodes[j+nH];
        myGrid.elements[i].eNodes[2]=&myNodes[j+nH+1];
        myGrid.elements[i].eNodes[3]=&myNodes[j+1];
        licznik++;
        j++;
        if (licznik==(nH-1)) {
            j++;
            licznik=0;
        }
    }
}

double *gauss(double **A, int n) {

    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    double *x = new double[n];
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}


