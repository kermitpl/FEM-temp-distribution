//
// Created by Adrian on 15.11.2018.
//

#ifndef MES_BASICSTRUCTURES_H
#define MES_BASICSTRUCTURES_H

#include <iostream>

using namespace std;

struct node{
    int ID;
    double x;
    double y;
    double t0;
    bool BC=false;
    void show()
    {
        //cout<<"ID="<<ID;
        //cout<<"["<<x<<" , "<<y<<"]"<<"("<<t0<<") \t";
        //cout<<t0<<"\t";
        if (BC) cout<<"1"<<"\t";
        else cout<<"0"<<"\t";
    }
};

struct element{
    node * eNodes [4];
};

struct grid{
    element * elements;
    void show(int nE)
    {
        for (int i=0; i<nE; i++)
        {
            cout<<"Element nr "<<i<<endl;
            for (int j=0; j<4; j++) {
                elements[i].eNodes[j]->show();
            }
            cout<<endl;
        }
    }
};

struct point{
    double x;
    double y;
};

#endif //MES_BASICSTRUCTURES_H
