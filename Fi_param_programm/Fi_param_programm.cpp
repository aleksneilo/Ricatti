//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
////      _____    _    _____     _____   _____       ////
////      \  __\  |_|   \  __\   |  __/   \  __\      ////
////       \ \    | |    \ \     | |__     \ \        ////
////      __\ \   | |   __\ \    | __/    __\ \       ////
////      \____\  |_|   \____\   |_|      \____\      ////
////                                                  ////
//////////////////////////////////////////////////////////
////////////////////////by R3ZZ///////////////////////////
//////////////////////////////////////////////////////////

#include <cstring>
#include "stdafx.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include "SFS.h"
#include <cmath>
#include <iomanip>
#include <thread>
#include <string>
#include <chrono>


#define pi 3.141592653589793
#define icom (complex<double>(0, 1.))

using namespace std;

//////// INITIALIZATION /////////

double Gb_I, Gb_FS, Gb_SF, R0A, R0AN, Ksi_F, Ksi_S, Ksi_N, Del0, Del0a, Del0b, Xi1, Xi2, T, w, w1, E1, E, ro_F, ro_N, ro_S, H, alphaT;//
double  L_SL, L_SR, L_S, L_N, L_F, L_S1, L_F1, L_S2, L_s;
int  w_obrez, iter, MODE;
int N_SL, N_SR, N_Mid, N_S, N_N, N_F, N_S1, N_F1, N_S2, N, N_CPR, NUM_Tech;
double epsG, epsDel, alpha;
double x;
double* Hi, * Ksii, * Roi, * Li, * Stype, * Tci, * Rbi;

int main()
{

    ///////// VARIABLES //////////////

        // Number of Layers

    NUM_Tech = 1;

    //   Number of points in grid for each layer
    N_Mid = 1000;// in self-sonst
    N_S = int(5 * N_Mid); N_F = int(0.5 * N_Mid); N_S1 = int(0.2 * N_Mid); N_F1 = int(0.6 * N_Mid); N_S2 = int(4.4 * N_Mid); N_N = int(5 * N_Mid);
    N = N_S + N_N;// +N_S1 + N_F1 + N_S2;// +N_N;//*NUM_Tech;
    L_S = N_S / (1. * N_Mid); L_F = N_F / (1. * N_Mid); L_S1 = N_S1 / (1. * N_Mid); L_F1 = N_F1 / (1. * N_Mid); L_S2 = N_S2 / (1. * N_Mid); L_N = N_N / (1. * N_Mid);
    //cout<<N_S<<" "<<N_F<<" "<<N_S1<<" "<<N_F1<<" "<<N_S2<<" "<<N<<endl;
    N_CPR = 1;

    Hi = new double[NUM_Tech];
    Ksii = new double[NUM_Tech]; // Normalized on S
    Roi = new double[NUM_Tech];  // Normalized on S
    Li = new double[NUM_Tech];
    Stype = new double[NUM_Tech];
    Tci = new double[NUM_Tech];
    Rbi = new double[NUM_Tech];
    double* DifDel, DifDelmax, DifDelmax_place, KIDP, KID;
    complex<double>* DelmaxP, * DelmaxAP, * DelPbuf, * DelAPbuf;
    DifDel = new double[126];
    int dd, findmax, stor, Ns = 1000000;
    // Initialization of arrays
    DelmaxP = new complex<double>[10 * N];
    DelmaxAP = new complex<double>[10 * N];
    DelPbuf = new complex<double>[10 * N];
    DelAPbuf = new complex<double>[10 * N];

    unsigned int n = std::thread::hardware_concurrency();
    std::cout << n << " concurrent threads are supported.\n";

    T = 0.5;
    //   Iteration constants11
    epsG = 1e-9;  // accuracy of normal Green function G iteration loop
    epsDel = 5e-6; // accuracy of pair potential \Delta iteration loop
    alpha = 0.59; // parameter in Delta loop. Just for increase of convergence
    iter = 0; // variable to count number of iterations in Delta loop
    MODE = 0;
    //   Material parameters
    Ksi_S = 1; // Coherence length of superconductor. Used for normalization of all other lengths and scales
    Ksi_F = 1;//3; // Coherence length of ferromagnet
    Ksi_N = 10;
    Xi1 = 0; // Phases of pair potential on the left (Xi1) and right (Xi2) electrodes
    Xi2 = -Xi1; // Josephson phase is their difference: Xi = Xi2 - Xi1
    w_obrez = int(30 / T); // Number of Matsubara frequencies used in calculation.
    H = 10; // Exchange field
    ro_S = 1;//520;  // Resistivity of superconductor
    ro_F = 1;//440;  // Resistivity of ferromagnet
    ro_N = 0.1;
    Gb_FS = 0.3; // boundary parameter of Ferromagnet/Superconductor interface. It depends from resistivity of the interface.
    R0A = Gb_FS * ro_F * Ksi_F;
    Gb_SF = Gb_FS * ro_F * Ksi_F / ro_S / Ksi_S; // it also deterimines Superconductor/Ferromagnet interface

    double DELS = 1, DELPmax, DELAPmax, dqmax, epsq, qbuf, roo, ksii, dss, Hext, IIP, II, IIPmax, IImax, I;
    int kk, iter1, itermax, itermax1, iterq, culc, dN;
    epsq = 1e-5;

    //_____________  Creation of OUTPUT file_______________________//
    string name1("6stLk(I)H="), name2("6stLk(X,I)H="), name3("2000stDSmax(X,ksi=ro,ds)H="), str1, str2, str3, str4, str5, str6, str7, str8, str9;// cheate one big file
    stringstream s1, s2, s3, s4, s5, s6, s7, s8, s9;
    s1 << H; s2 << T; s3 << L_S1; s4 << L_F; s5 << L_F1; s9 << L_S2; s6 << Gb_FS; s7 << Ksi_F; s8 << ro_F;
    s1 >> str1; s2 >> str2; s3 >> str3; s4 >> str4; s5 >> str5; s6 >> str6; s7 >> str7; s8 >> str8; s9 >> str9;
    name1.append(str1); name1.append("T="); name1.append(str2); name1.append("LS1="); name1.append(str3); name1.append("LF="); name1.append(str4); name1.append("LF1="); name1.append(str5); name1.append("Ls="); name1.append(str9); name1.append("gb="); name1.append(str6); name1.append("KsiF="); name1.append(str7); name1.append("roF="); name1.append(str8);                                            name1.append(".txt");
    name2.append(str1); name2.append("T="); name2.append(str2); name2.append("LS1="); name2.append(str3); name2.append("LF="); name2.append(str4); name2.append("LF1="); name2.append(str5); name2.append("Ls="); name2.append(str9); name2.append("gb="); name2.append(str6); name2.append("KsiF="); name2.append(str7); name2.append("roF="); name2.append(str8);                                            name2.append(".txt");
    name3.append(str1); name3.append("T="); name3.append(str2); name3.append("LS1="); name3.append(str3); name3.append("LF="); name3.append(str4); name3.append("LF1="); name3.append(str5); name3.append("Ls="); name3.append(str9); name3.append("gb="); name3.append(str6); name3.append("KsiF="); name3.append(str7); name3.append("roF="); name3.append(str8);                                            name3.append(".txt");
    const char* file1 = name1.c_str();
    const char* file2 = name2.c_str();
    const char* file3 = name3.c_str();
    ofstream fout1(file1);
    ofstream fout2(file2);
    //ofstream fout3(file3);
    double x0, x1, x2; char chh;
    ifstream file("11e5DELth(X,dF)H=100T=0.5LS1=0LF=0LF1=0Ls=0gb=0.3KsiF=2.5roF=1.txt");
    if (file.is_open()) cout << "good\n";
    else cout << "beed\n\n" << endl;



    //for (double hh = 1.; hh < 1.01; hh += 2.)
    //{
    Del0 = SelfConsZero();

    //for (Ksi_N = 1; Ksi_N < 10.01; Ksi_N += 19.)
        //for (N_N = int(3 * N_Mid); N_N<int(10.01 * N_Mid); N_N += int(10.1 * N_Mid))
    for (N_F = int(0.5 * N_Mid); N_F<int(1. * N_Mid); N_F += int(10.05 * N_Mid))
        //if((abs(L_F-0.5)<1e-5))//||(abs(L_F-0.6)<1e-5)||(abs(L_F-0.7)<1e-5)||(abs(L_F-0.8)<1e-5)||(abs(L_F-0.9)<1e-5)||(abs(L_F-1.1)<1e-5))
    {
        for (int i = 0; i < int(N_Mid * L_S); i++) { DelPbuf[i] = Del0; DelAPbuf[i] = Del0; }
        DifDelmax_place = 4.4;
        double** StorageDEL = new double* [7];  //ds    DELP    DELAP   ||DELP|-|DELAP||
        for (int i = 0; i < 7; i++) StorageDEL[i] = new double[150];

        L_S1 = N_S1 / (1. * N_Mid);
        N_F1 = N_F + int(0.1 * N_Mid);
        L_F = N_F / (1. * N_Mid);
        L_N = N_N / (1. * N_Mid);
        L_F1 = N_F1 / (1. * N_Mid);//*/
        R0A = Gb_FS * ro_F * Ksi_F;
        R0AN = Gb_FS * ro_N * Ksi_N;
        dN = int(0.0 * N_Mid); N_S2 = int((DifDelmax_place + 0.0) * N_Mid); L_S2 = N_S2 / (1. * N_Mid);//+55*4;
        DifDelmax = 0; stor = 0; kk = 0;

        for (int i = 0; i < 7; i++) for (int j = 0; j < 150; j++) StorageDEL[i][j] = 0.;
        while ((kk < 1))
        {
            findmax = 0; DifDelmax = 0;
            while ((findmax == 0) && (L_S2 > 1.))//&&(abs(DELS)>0.02))//
            {

                N = N_S + N_F + N_S1 + N_F1 + N_S2 + N_N;
                L_S2 = N_S2 / (1. * N_Mid);
                cout << N_S << " " << N_F << " " << N_S1 << " " << N_F1 << " " << N_S2 << " " << N_N << " " << N << endl;
                for (int i = 0; i < NUM_Tech; i++)
                {
                    if (i == 0)
                    {
                        Hi[i] = 0; Ksii[i] = Ksi_S; Roi[i] = ro_S; Li[i] = L_S; Stype[i] = 1.; Tci[i] = 1; Rbi[i] = R0A;
                    }
                    if ((i == 1))//||(i == 5))
                    {
                        Hi[i] = H; Ksii[i] = Ksi_F; Roi[i] = ro_F; Li[i] = L_F; Stype[i] = 0.;  Tci[i] = 0.; Rbi[i] = R0A;
                    }
                    if ((i == 2))//||(i == 4)||(i == 6))
                    {
                        Hi[i] = 0; Ksii[i] = Ksi_S; Roi[i] = ro_S; Li[i] = L_S1; Stype[i] = 1.; Tci[i] = 1; Rbi[i] = R0A;
                    }//*/
                    if ((i == 3) || (i == 7))
                    {
                        Hi[i] = H; Ksii[i] = Ksi_F; Roi[i] = ro_F; Li[i] = L_F1; Stype[i] = 0.;  Tci[i] = 0.; Rbi[i] = R0A;
                    }
                    if (i == 4)
                    {
                        Hi[i] = 0; Ksii[i] = Ksi_S; Roi[i] = ro_S; Li[i] = L_S2; Stype[i] = 1.; Tci[i] = 1; Rbi[i] = R0A;
                    }
                    if (i == 5)
                    {
                        Hi[i] = 0; Ksii[i] = Ksi_N; Roi[i] = ro_N; Li[i] = L_N; Stype[i] = 0.; Tci[i] = 0; Rbi[i] = R0AN;
                    }
                }
                complex<double>* Del, * DelP, * Delbuf, * G, * Fi, * Fi_old, * Fi1, * Fi1_old;
                double* Is, * q, * qP, * qbuf;
                Del = new complex<double>[N];
                DelP = new complex<double>[N];
                Delbuf = new complex<double>[N];
                Is = new double[N];
                q = new double[N];
                qP = new double[N];
                qbuf = new double[N];
                G = new complex<double>[N];
                //G1 = new complex<double>[10*N];
                Fi = new complex<double>[N];
                Fi_old = new complex<double>[N];
                Fi1 = new complex<double>[N];
                Fi1_old = new complex<double>[N];


                Del0 = SelfConsZero();
                culc = 1;
                for (int i = 0; i < 150; i++)
                    if (abs(StorageDEL[0][i] - L_S2) < 1e-5)
                    {
                        DelP[N_S + N_F + N_S1 + N_F1 + N_S2 / 2] = StorageDEL[1][i]; Del[N - N_N - 1] = StorageDEL[2][i]; iter1 = int(StorageDEL[3][i]); iter = int(StorageDEL[4][i]); IIP = StorageDEL[5][i]; II = StorageDEL[6][i];
                        culc = 0; // we found the last calculated at the ds lenght
                    }
                if (culc == 1)
                {

                    for (int i = 0; i < N; i++)
                    {
                        DelP[i] = get_type(i) * Del0;
                        if (i > N_S - 1) DelP[i] = -0. * get_type(i) * Del0;
                        Del[i] = get_type(i) * Del0;
                        if (i > N_S - 1) Del[i] = -0. * get_type(i) * Del0;
                        qP[i] = 0.;
                        q[i] = 0.;
                    }//*/


                    for (I = 0.0; I < 0.201; I += 0.001)
                        //if ((abs(I - 0.) < 1e-5)|| (abs(I - 0.1) < 1e-5) || (abs(I - 0.3) < 1e-5) )
                    {
                        ////////////////PARALLEL ORIENTATION///////////

                        Hi[1] = H;
                        dqmax = 1.; iterq = 0;
                        while ((dqmax > epsq))//&&(iterq==0))
                        {
                            for (int i = 1; i < N - 1; i++) qbuf[i] = qP[i];
                            //for (int i = 0; i < N; i++) G[i] = 1. + 0.1 * icom;
                                //for (int i = 0; i < N; i++) { DelP[i] = get_type(i) * DelPbuf[i]; Del[i] = get_type(i) * DelAPbuf[i];}
                            MODE = 0;
                            SelfConsParal(G, DelP, 0, qP, I);     // Self-consistent  calculation of pair potential.
                            iter1 = iter;
                            KIDP = KinInd(G, DelP, qP, I);
                            //for (int i = 1; i < N - 1; i++) cout << i << "  "  << I << "  "  << qP[i] << endl;
                            //for (int i = 0; i < N; i++) qP[i] = 0.1*qP[i]+0.9* qbuf[i];
                            iterq++;
                            dqmax = 0.;
                            for (int i = 1; i < N - 1; i++)
                            {
                                if (dqmax < abs(qbuf[i] - qP[i]))
                                    dqmax = abs(qbuf[i] - qP[i]);
                            }
                            cout << fixed << iterq << "  " << "dqmax=" << dqmax << "  " << KIDP << "  " << DelP[N - N_N - 1] << endl;
                        }
                        fout1 << fixed << I << "   " << L_N << "  " << Ksi_N << "  " << ro_N << "   " << real(DelP[0]) << "   " << real(DelP[N_S - 1]) << "   " << KIDP << "   " << qP[2] << endl;// "   " << IIP << "   " << II << "  " << abs(abs(II) - abs(IIP)) << "  " << culc << endl;// "  " << KIDP << "  " << KID << "  " << iter1 << "  " << iter << "  " << culc << endl;

                        /*double dGmax = 1; IIP = 0;
                        for (int iw = 0; iw < 1; iw++)
                        {
                            dGmax = 1;
                            w = -pi * T * (2. * iw + 1.);

                            while ((dGmax > epsG))
                            {
                                //for (int i=0; i<N; i++) Fi_old[i]= Fi[i];
                                Prog(Fi, G, DelP, w, qP, I);
                                /*for (int i = 0; i<N; i++)
                                {	G[i]= -conj(G[i]);
                                    Fi[i]=0.5*Fi[i]+0.5*Fi_old[i];
                                }
                                //for (int i=0; i<N; i++)	Fi1_old[i]= Fi1[i];
                                //Prog( Fi1, G, DelP, -w);
                                for (int i = 0; i < N; i++)
                                {
                                    //G[i]= -conj(G[i]);
                                    Fi1[i] = conj(Fi[i]);
                                }
                                Gcalc(G, &dGmax, Fi, Fi1, DelP, w, qP);
                                cout << dGmax << endl;
                                //cout<<fixed<<dGmax<<"  "<<abs(Fi[0]/get_wm(0,w))<<"  "<<abs(G[0]*Fi[0]/get_wm(0,w))<<"  "<<abs(G[N/2]*Fi[N/2]/get_wm(N/2,w))<<"  "<<abs(G[N-1]*Fi[N-1]/get_wm(N-1,w))<<endl;
                            }
                            for (int n = 0; n < N; ++n)	fout2 << fixed << n << "  " << I << "  " << real(DelP[n]) << "  " << qP[n] << "  " <<abs(Fi[n]) << "  " << real(Fi[n]) << "  " << imag(Fi[n]) << endl;
                        }/*/

                    }

                    Hi[1] = -H;
                    dqmax = 10.; iterq = 0;
                    for (I = 0.0; I < 0.201; I += 0.001)
                        //if ((abs(I - 0.) < 1e-5)|| (abs(I - 0.1) < 1e-5) || (abs(I - 0.3) < 1e-5) )
                    {
                        Hi[1] = -H;
                        dqmax = 1.; iterq = 0;
                        while ((dqmax > epsq) && (iterq == 0))
                        {
                            for (int i = 1; i < N - 1; i++) qbuf[i] = q[i];
                            //for (int i = 0; i < N; i++) G[i] = 1. + 0.1 * icom;
                                //for (int i = 0; i < N; i++) { DelP[i] = get_type(i) * DelPbuf[i]; Del[i] = get_type(i) * DelAPbuf[i];}
                            MODE = 0;
                            SelfConsParal(G, Del, 0, q, I);     // Self-consistent  calculation of pair potential.
                            iter1 = iter;
                            KID = KinInd(G, Del, q, I);
                            iterq++;
                            dqmax = 0.;
                            for (int i = 1; i < N - 1; i++)
                            {
                                if (dqmax < abs(qbuf[i] - q[i]))
                                    dqmax = abs(qbuf[i] - q[i]);
                            }

                            cout << fixed << iterq << "  " << "dqmax=" << dqmax << "  " << KID << "  " << Del[N - N_N - 1] << endl;

                        }
                        fout1 << fixed << I << "   " << L_N << "  " << Ksi_N << "  " << ro_N << "   " << real(Del[0]) << "   " << real(Del[N_S - 1]) << "   " << KID << "   " << q[2] << endl;// "   " << IIP << "   " << II << "  " << abs(abs(II) - abs(IIP)) << "  " << culc << endl;// "  " << KIDP << "  " << KID << "  " << iter1 << "  " << iter << "  " << culc << endl;


                        //for (int n = 0; n < N; ++n)	fout2 << fixed << n << "  " << I << "  " << real(DelP[n]) << "  " << real(Del[n])<< "  " << qP[n] << "  " << q[n] << endl;
                        //fout1 << fixed << I << "   " << L_S2 << "   " << L_N << "  " << Ksi_N << "  " << L_F << "  " << L_S1 << "  " << ro_N << "   " << real(DelP[N_S-1]) << "   " << real(Del[N_S-1]) << "   " << real(DelP[N - N_N - 1]) << "   " << real(Del[N - N_N - 1]) << "   " << abs(abs(real(Del[N - N_N - 1])) - abs(real(DelP[N - N_N - 1]))) << "  " << KIDP << "  " << KID << "  " << (KIDP - KID) << endl;// "   " << IIP << "   " << II << "  " << abs(abs(II) - abs(IIP)) << "  " << culc << endl;// "  " << KIDP << "  " << KID << "  " << iter1 << "  " << iter << "  " << culc << endl;
                        //fout1 << fixed << I << "   " << L_N << "  " << Ksi_N << "  " << ro_N << "   " << real(DelP[0]) << "   " << real(DelP[N_S - 1]) << "   " << KIDP  << endl;// "   " << IIP << "   " << II << "  " << abs(abs(II) - abs(IIP)) << "  " << culc << endl;// "  " << KIDP << "  " << KID << "  " << iter1 << "  " << iter << "  " << culc << endl;

                    }
                    //for (int n = 0; n < N; ++n)	fout2 << fixed << n << "  " << I << "  " << real(DelP[n]) << "  " << q[n] << endl
                    /*double dGmax = 1; IIP = 0;
                     for (int iw = 0; iw < 1; iw++)
                     {
                         dGmax = 1; iterG = 0;
                         w = -pi * T * (2. * iw + 1.);

                         while ((dGmax > epsG))
                         {
                             //for (int i=0; i<N; i++) Fi_old[i]= Fi[i];
                             Prog(Fi, G, DelP, w, q, I);
                             /*for (int i = 0; i<N; i++)
                             {	G[i]= -conj(G[i]);
                                 Fi[i]=0.5*Fi[i]+0.5*Fi_old[i];
                             }//*
                             //for (int i=0; i<N; i++)	Fi1_old[i]= Fi1[i];
                             //Prog( Fi1, G, DelP, -w);
                             for (int i = 0; i < N; i++)
                             {
                                 //G[i]= -conj(G[i]);
                                 Fi1[i] = conj(Fi[i]);
                             }
                             Gcalc(G, &dGmax, Fi, Fi1, DelP, w); iterG++;
                             //cout<<fixed<<dGmax<<"  "<<abs(Fi[0]/get_wm(0,w))<<"  "<<abs(G[0]*Fi[0]/get_wm(0,w))<<"  "<<abs(G[N/2]*Fi[N/2]/get_wm(N/2,w))<<"  "<<abs(G[N-1]*Fi[N-1]/get_wm(N-1,w))<<endl;
                         }
                         //IIP+=real(Del0/sqrt(w*w+Del0*Del0)*G[N-1]*Fi[N-1]/get_wm(N-1,w));
                         for (int n = 0; n < N; ++n)	fout2 << fixed << n << "  " << real(DelP[n]) << "  " << abs(G[n] * Fi[n] / get_wm(n, w)) << "  " << real(G[n] * Fi[n] / get_wm(n, w)) << "  " << imag(G[n] * Fi[n] / get_wm(n, w)) << "  " << arg(G[n] * Fi[n] / get_wm(n, w)) << "  " << abs(G[n]) << endl;
                         //fout1<<fixed<<iw<<"  "<<arg(G[N-1]*Fi[N-1]/get_wm(N-1,w))<<endl;
                     }//*/
                     // Hi[1] = -H;
                      //SelfConsParal(G, Del, 0, q, I);     // Self-consistent  calculation of pair potential.
                      //KID =  KinInd(G, Del, q, I);
                      //for (int i = 0; i < N; i++) { DelPbuf[i] = DelP[i]; DelAPbuf[i] = Del[i]; }
                      //fout1 << fixed << iterq << "  " << "dqmax=" << dqmax << "  " << KIDP << "  " << DelP[N - 1] << endl;
                      //for (int n = 0; n < N; ++n)	cout << fixed << n << "  " << KIDP << "  " << real(DelP[n]) << "  " << q[n] << endl;
                      /*iterq++;
                      dqmax = 0.;
                      for (int i = 1; i < N-1; i++)
                      {
                          if (dqmax < abs(qbuf[i] - q[i]))
                              dqmax = abs(qbuf[i] - q[i]);
                      }

                      cout << fixed << iterq << "  " << "dqmax=" << dqmax << "  " << KIDP<< "  " << DelP[N - 1] << endl;/*/
                      //}
                      //for (int n = 0; n < N; ++n)	fout2 << fixed << n << "  " << I << "  " << real(DelP[n]) << "  " << q[n] << endl;
                      //MODE=0;
                      //E_IS_calc(Is,G,DelP);
                      //for (int i=0; i<N; i++) cout<<i<<"  "<<Del[i]<<"  "<<get_type(i)<<"  "<<get_ksi(i)<<"  "<<get_wm(i,1)<<"  "<<Roi[Layer(i)]<<"  "<<Layer(i)<<endl;
                      //int iterG;
                       /*dGmax = 1; IIP = 0;
                      for(int iw=0; iw<1; iw++)
                      {       dGmax=1; iterG=0;
                              w=-pi*T*(2.*iw+1.);

                              while ((dGmax > epsG))
                              {
                                          //for (int i=0; i<N; i++) Fi_old[i]= Fi[i];
                                          Prog( Fi, G, Del, w,q,I);
                                          /*for (int i = 0; i<N; i++)
                                          {	G[i]= -conj(G[i]);
                                              Fi[i]=0.5*Fi[i]+0.5*Fi_old[i];
                                          }//*
                                          //for (int i=0; i<N; i++)	Fi1_old[i]= Fi1[i];
                                          //Prog( Fi1, G, DelP, -w);
                                          for (int i=0; i<N; i++)
                                          {
                                              //G[i]= -conj(G[i]);
                                              Fi1[i]=conj(Fi[i]);
                                          }
                                          Gcalc( G, &dGmax, Fi, Fi1, Del, w); iterG++;
                                          //cout<<fixed<<dGmax<<"  "<<abs(Fi[0]/get_wm(0,w))<<"  "<<abs(G[0]*Fi[0]/get_wm(0,w))<<"  "<<abs(G[N/2]*Fi[N/2]/get_wm(N/2,w))<<"  "<<abs(G[N-1]*Fi[N-1]/get_wm(N-1,w))<<endl;
                              }
                              //IIP+=real(Del0/sqrt(w*w+Del0*Del0)*G[N-1]*Fi[N-1]/get_wm(N-1,w));
                              for(int n=0; n<N;++n)	fout2<<fixed<<n<<"  "<<real(Del[n])<<"  "<<abs(G[n]*Fi[n]/get_wm(n,w))<<"  "<<real(G[n]*Fi[n]/get_wm(n,w))<<"  "<<imag(G[n]*Fi[n]/get_wm(n,w))<<"  "<<arg(G[n]*Fi[n]/get_wm(n,w))<<"  "<<abs(G[n])<<endl;
                              //fout1<<fixed<<iw<<"  "<<arg(G[N-1]*Fi[N-1]/get_wm(N-1,w))<<endl;
                      }//for(int n=0; n<N;++n)	fout2<<fixed<<n<<"  "<<real(Del[n])<<endl;/*/
                      //fout2<<fixed<<Xi2<<"  "<<abs(II)<<"  "<<real(II)<<endl;
                      //}

                      //////////////ANTIPARALLEL ORIENTATION//////////
                      //Hi[1] = -H; //Hi[7]=-H;

                      //DELS = real(DelP[N_S + N_F + N_S1 + N_F1 + N_S2 / 2]);
                      //StorageDEL[0][stor] = L_S2; StorageDEL[1][stor] = real(DelP[N_S + N_F + N_S1 + N_F1 + N_S2 / 2]); StorageDEL[2][stor] = real(Del[N - N_N - 1]); StorageDEL[3][stor] = double(iter1); StorageDEL[4][stor] = double(iter);// StorageDEL[5][stor] = IIP; StorageDEL[6][stor] = II;
                      //stor++;
                      //cout << "DELS=" << DELS << endl;
                  //}//*/

                  //fout1<<fixed<<L_F1<<"   "<<Xi2/pi<<"   "<<real(Del[0])<<"   "<<imag(Del[0])<<"   "<<real(Del[N-1])<<"   "<<imag(Del[N-1])<<"   "<<Is[N/2]<<endl;
                  //fout1 << fixed << L_S2 << "   " << L_N << "  " << Ksi_N << "  " << L_F << "  " << L_S1 << "  " << ro_F << "   " << real(DelP[N - N_N - 1]) << "   " << real(Del[N - N_N - 1]) << "   " << abs(abs(real(Del[N - N_N - 1])) - abs(real(DelP[N - N_N - 1]))) << "  " << KIDP << "  " << KID << "  " << (KIDP - KID)<<  endl;// "   " << IIP << "   " << II << "  " << abs(abs(II) - abs(IIP)) << "  " << culc << endl;// "  " << KIDP << "  " << KID << "  " << iter1 << "  " << iter << "  " << culc << endl;

                  /*if (abs(real(DelP[N_S + N_F + N_S1 + N_F1 + N_S2 / 2])) < 0.1)//(abs(abs(real(DelP[N -N_N- 1])) - abs(real(Del[N -N_N- 1]))) > DifDelmax)
                  {
                      DifDelmax = abs(abs(real(DelP[N - N_N - 1])) - 0.02 - 0. * abs(real(Del[N - N_N - 1])));
                      DifDelmax_place = L_S2;
                      itermax1 = iter1; itermax = iter;
                      //IIPmax = IIP; IImax = II;
                      DELPmax = real(DelP[N - N_N - 1]); DELAPmax = real(Del[N - N_N - 1]);
                      //for (int i = 0; i < N; i++) DelmaxP[i] = DelP[i];
                      //for (int i = 0; i < N; i++) DelmaxAP[i] = Del[i];
                      //for(int i=0; i<N; i+=100) cout<<DelmaxAP[i]<<"  "Del[i]<<endl;
                      findmax = 1;
                  }
                  else findmax = 0;
                  findmax=1;  //stop calc, we find max in this k-iteration
                  if (kk == 0) N_S2 -= dN;
                  if (kk == 1) N_S2 -= dN / 2;
                  if (kk == 2) N_S2 -= dN / 4;
                  if (kk == 3) N_S2 -= dN / 8;
                  //if(kk==3) N_S2-=1*50;
                  //if((kk==2)&&(findmax==1)) { cout<<"k=========="<<kk;
                  //for(int i=0; i<N; i++) fout2<<fixed<<i<<"  "<<L_N<<"   "<<Ksi_N<<" "<<real(DelP[i])<<" "<<real(Del[i])<<endl;
                  //for(int i=0; i<N; i++) fout3<<fixed<<i<<"  "<<ro_F<<"   "<<Ksi_F<<" "<<real(DelmaxAP[i])<<endl;}//*/
                    findmax = 1;
                    delete[] G;
                    delete[] Del;
                    delete[] DelP;
                    delete[] Delbuf;
                    delete[] Fi;
                    delete[] Fi_old;
                    delete[] Fi1;
                    delete[] Fi1_old;
                    delete[] q;
                    delete[] qP;
                    delete[] qbuf;
                }
                if (kk == 0) N_S2 += (3 * dN);
                if (kk == 1) N_S2 += (3 * dN) / 2;
                if (kk == 2) N_S2 += (3 * dN) / 4;
                //if(kk==0) N_S2+=(5*dN)/2;if(kk==1) N_S2+=(5*dN)/8;if(kk==2) N_S2+=(5*dN)/16;
                kk++;//*/
            }
        }
        //                        fout2 << fixed << L_N << "  " << L_S2 << "  " << Ksi_N << " " << DELS << " " << DifDelmax << "   " << DifDelmax_place << "   " << DELPmax << "   " << DELAPmax << endl;// "  " << abs(abs(IImax) - abs(IIPmax)) << "   " << IIPmax << "   " << IImax << "   " << itermax1 << "   " << itermax << endl;
                                //fout2<<fixed<<ro_F<<"   "<<Ksi_F<<" "<<DifDelmax<<"   "<<DifDelmax_place<<"   "<<DELPmax<<"   "<<DELAPmax<<"   "<<itermax1<<"   "<<itermax<<endl;
        for (int i = 0; i < 7; i++)
            delete[] StorageDEL[i];
        delete[] StorageDEL;

    }

    //}//}//}//}
    delete[] DelmaxP;
    delete[] DelmaxAP;
    delete[] DelPbuf;
    delete[] DelAPbuf;


    fout1 << "\n" << "H=" << H << ", T=" << T << ", L_S=" << L_S << ", L_F=" << L_F << ", L_S1=" << L_S1 << ", L_F1=" << L_F1 << ", L_S2=" << L_S2 << ", Gb_FS=" << Gb_FS << ", ksi_F=" << Ksi_F << ", ro_F=" << ro_F << ", R0A=" << R0A << ", Nmid=" << N_Mid << "\n";
    fout2 << "\n" << "H=" << H << ", T=" << T << ", L_S=" << L_S << ", L_F=" << L_F << ", L_S1=" << L_S1 << ", L_F1=" << L_F1 << ", L_S2=" << L_S2 << ", Gb_FS=" << Gb_FS << ", ksi_F=" << Ksi_F << ", ro_F=" << ro_F << ", R0A=" << R0A << ", Nmid=" << N_Mid << "\n";



    system("PAUSE");
    return EXIT_SUCCESS;
}
