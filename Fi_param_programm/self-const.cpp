
#include <cstring>
#include "stdafx.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include "SFS.h"
#include <thread>
#include <string>
#include <chrono>
#include <vector>
 #define pi 3.141592653589793
 #define icom (complex<double>(0, 1.))

 using namespace std;


//////// SelfCons - function to calculate pair potential ////////////

//   INPUT Requires Global variables: Scales and Del0 for initial pair potential
//   OUTPUT *G - pointer to array with Normal Green Function G
//   OUTPUT *Del - pointer to array with Pair Potential Delta

void SelfCons(complex<double> *G, complex<double> *Del, int Initial, double* q, double I)
    {



	double w, dGmax, dDelmax;
	int iterG;//, ii=N_S+N_F+N_S1/2;
     complex<double> *S1, *S2, *Fi, *Fi1, Delbuf;
	 double *Tc;
     complex<double> *Fi_old, *Fi1_old;

     Fi=new complex<double>[N];  // Parametrized Anomalous Green Function Fi(w) for positive Matsubara
     Fi1=new complex<double>[N]; // Parametrized Anomalous Green Function Fi(-w)  for negative Matsubara

     Fi_old=new complex<double>[N];  // Parametrized Anomalous Green Function Fi(w) for positive Matsubara
     Fi1_old=new complex<double>[N]; // Parametrized Anomalous Green Function Fi(-w) for negative Matsubara

     S1=new complex<double>[N]; // S1, S2 - Inner SUMS in self-consistent equation
	 S2 = new complex<double>[N];
	//double S2 ;                // S2 - SUM multiplied by Delta and S1 is free SUM
     dDelmax=1;

	 Tc = new double[N];
	 char TXTNAME1[40],buf1[40];
//  	fout.open("file.txt");

//	 double Xi3= 2*pi*(jX-1)/3;

     // Set Initial Delta And G values
	/*if (Initial == 0)
    	for (int i=0; i<N; i++)
        {
    	    G[i]=1.;
			Del[i]= get_type(i)*Del0;
			if(i>N_S-1) Del[i]= -0.*get_type(i)*Del0;//*exp(icom*0.01*pi/2.);+N_F+N_S1+N_F1
			Tc[i] = get_tc(i);
        }//*/
	if (Initial == 2)
		for (int i = 0; i<N; i++)
		{
			G[i] = 1.;
			Del[i] = 1e-2;
			Tc[i] = get_tc(i);
		}
    //for (int i=0; i<N; i++) cout<<i<<"  "<<Del[i]<<"  "<<get_type(i)<<"  "<<get_ksi(i)<<"  "<<get_wm(i,1)<<"  "<<Roi[Layer(i)]<<"  "<<Layer(i)<<endl;
    // Start self-consistent loop
   	iter=0;
    		//cout<<"1	"<<endl;
while ((dDelmax > epsDel)&&(iter<1000))
{
	   for (int i = 0; i < N; i++)
	   {
		   S1[i] = 0;
		   S2[i] = log(T/Tc[i]) / pi / (T);
	   }
        // Summarize over all Matsubara frequencies
        for(int iw=0; iw<w_obrez; iw++)
        {
            w=pi*T*(2.*iw+1.);
            dGmax=1;
            iterG=0;
		     // Iterative loop over G. Solve Nonlinear Usadel equation here
            while ((dGmax > epsG)&&(iterG<1000))
            {
                // Calculation of Anomalous Green functions Fi(w), Fi(-w)
				for (int i=0; i<N; i++)
				{
                	Fi_old[i]= Fi[i];
   //             	cout<<"Fi_old["<<i<<"]="<<Fi_old[i]<<" "<<'\n';
     //           	cout<<"iw="<<iw<<" "<<"T="<<T<<'\n';
            //    	}
				}
				Prog( Fi, G, Del, w, q, I);
                for (int i=0; i<N; i++)
                {
					G[i]= -conj(G[i]);
					Fi[i]=0.5*Fi[i]+0.5*Fi_old[i];//????
				//		for(int j=0;j<N;j++)

	//						cout<<"Fi["<<i<<"]="<<Fi[i]<<" "<<"G["<<i<<"]="<<G[i]<<'\n';
	//						cout<<"iw="<<iw<<" "<<"T="<<T<<'\n';
				//		}
				}

				for (int i=0; i<N; i++)
                {
					Fi1_old[i]= Fi1[i];
               //   	for(int j=0;j<N;j++)
        //          	cout<<"Fi1_old["<<i<<"]="<<Fi1_old[i]<<" "<<'\n';
      //            	cout<<"iw="<<iw<<" "<<"T="<<T<<'\n';
            	}
				Prog( Fi1, G, Del, -w, q, I);
            	for (int i=0; i<N; i++)
				{
					G[i]= -conj(G[i]);
					Fi1[i]=0.3*Fi1[i]+0.7*Fi1_old[i];//????
		//			cout<<"Fi1["<<i<<"]="<<Fi1[i]<<" "<<"G["<<i<<"]="<<G[i]<<'\n';
		//			cout<<"iw="<<iw<<" "<<"T="<<T<<'\n';
				}
			//	cout<<"T="<<T<<'\n';
                // Recalculate G by resulting Fi(w) and Fi(-w)
                Gcalc( G, &dGmax, Fi, Fi1, Del, w,q);

			//of
			//cout<<"dGmax= "<< dGmax<<'\n';
            //iterG++;
			}
			//cout<<w<<"  "<<dGmax<<endl;

         //   }

                // search summ of self-cons eq. for each point of the grid

                for (int i=0; i<N; i++)
				{
                    S1[i] =  S1[i]+ get_tc(i)*get_type(i)*(G[i]*Fi[i]/w + conj(G[i])*Fi1[i]/w);
					S2[i] += 2. / w;
				}

		//	fout2<<"Fi_old["<<iw<<"]="<<Fi_old[iw]<<" "<<"Fi["<<iw<<"]="<<Fi[iw]<<" "<<"Fi1["<<iw<<"]="<<Fi1_old[iw]<<" "<<"T="<<T<<" "<<"iw="<<iw<<" "<<"iter="<<iter<<'\n';
		//	cout<<"Fi_old["<<iw<<"]="<<Fi_old[iw]<<" "<<"Fi["<<iw<<"]="<<Fi[iw]<<" "<<"Fi1["<<iw<<"]="<<Fi1_old[iw]<<" "<<"T="<<T<<" "<<"iw="<<iw<<" "<<"iter="<<iter<<'\n';'\n';
				//	for(int j=0;j<N;j++)
		//		{
		//			fout<<"Fi_old["<<j<<"]="<<Fi_old[j]<<" "<<"Fi["<<j<<"]="<<Fi[j]<<" "<<"Fi1["<<j<<"]="<<Fi1_old[j]<<" "<<"T="<<T<<'\n';

		//			cout2<< iw << " "<< Fi_old[iw]<<" " <<Fi[iw]<<" "<< Fi1[iw]<<" "<<T<<'\n';
			//		fout<<"Fi_old["<<j<<"]="<<Fi_old[j]<<" "<<"Fi["<<j<<"]="<<Fi[j]<<" "<<"Fi1["<<j<<"]="<<Fi1_old[j]<<" "<<"T="<<T<<'\n';
			//		cout<< iw << " "<< Fi_old[iw]<<" " <<Fi[iw]<<" "<< Fi1[iw]<<" "<<T<<'\n';
       	//		}
        }
        // Find delta from S1 and S2 and check mismatch with previous step
        dDelmax=0;
        for (int i=0; i<N; i++)
        {
    //    	cout<<"4"<<'\n';
            Delbuf= Del[i];
            //Del[i]= S1[i]/S2;
            Del[i]= (S1[i]/S2[i]+Delbuf*(alpha-1))/alpha;
            if (dDelmax< abs(Del[i]-Delbuf))
                dDelmax = abs(Del[i]-Delbuf);
        }

	cout<<iter<<" "<< dDelmax<<" "<<real(Del[0])<<" "<<imag(Del[0])<<" "<<real(Del[N-1])<<" "<<imag(Del[N-1])<<endl;
        iter++;   // try it until convergence
  //      cout<<"iter= "<<iter<<" dDelmax="<< dDelmax<<'\n';
    //
	//	}
      }
      //cout<<iter<<" "<< dDelmax<<" "<<real(Del[0])<<endl;
      delete [] S1;
	  delete [] S2;
      delete [] Fi;
      delete [] Fi1;
	  delete[] Fi_old;
	  delete[] Fi1_old;
	  delete[] Tc;
    }

///////////////////////////////////////parallel programming
void SelfConsParal(complex<double> *G, complex<double> *Del, int Initial, double* q, double I)
{
     complex<double> Delbuf,S2, *SS;
     int *wth, Ns=0,//N_S+N_F+N_S1+N_F1+N_S2-2,
     NS1=0;//N_S+N_F+N_S1/2;

     SS=new complex<double> [N]; 			//sum of S1 in each x
     int amount_of_threads=8;
     wth=new int [amount_of_threads];
        for (int i=0; i<amount_of_threads-1; i++)
            wth[i]=int((0.3*i/7. + 0.69*i/7.*i/7.)*w_obrez);
        wth[1]=1; wth[amount_of_threads-1]=w_obrez;
        for (int i=0; i<amount_of_threads; i++) cout<<wth[i]<<endl;

     complex<double>** S1 = new complex<double> *[w_obrez];
     for (int i = 0; i < w_obrez; i++)
	   S1[i] = new complex<double>[N];

     //int imax; 
     int w0 = 0, w1 = 1;//, w2=2, w3=3, w4=5, w5=7, w6=9, w7=12, w8=15, w9=18, w10=22, w11=26, w12=30, w13=35, w14=39, w15=44, w16=49, w17=w_obrez;//54,w18=w_obrez;//for T=0.5
	double dDelmax=1, w;
     // Set Initial Delta And G values
	//if (Initial == 0)
    	/*for (int i = 0; i<N; i++)
        {
			Del[i]= get_type(i)*Del0;
            if ((i > N_S)) Del[i] = -0.;// *Hi[1] / abs(Hi[1]) * get_type(i) * Del0;//*exp(icom*Xi2*pi/2.);//*Hi[1]/abs(Hi[1])*Del0;//*exp(icom*Xi2);;//*exp(icom*0.01*pi/2.);+N_F+N_S1+N_F1
            //if (i > (N - N_N - 1)) Del[i] = -get_type(i) * Del0;// *exp(icom * Xi2 * pi);
        }//*/
	
    //for (int i=0; i<N; i++) cout<<fixed<<i<<"  "<<Del[i]<<"  "<<get_type(i)<<"  "<<get_ksi(i)<<"  "<<get_wm(i,1)<<"  "<<Roi[Layer(i)]<<"  "<<Rbi[Layer(i)]<<"  "<<Layer(i)<<endl;
    // Start self-consistent loop
   	iter=0;

while ((dDelmax > epsDel)&& (iter < 2000))//((dDelmax<0.1)||((iter<70)))&&(iter<3000))
{
    // Find delta from S1 and S2 and check mismatch with previous step
        //epsG=pow(10.,-6-3*iter/50.); if(iter>50)
        //epsG=1e-11;
     	for (int iw=0; iw<w_obrez; iw++) for (int j=0; j<N; j++)   S1[iw][j]=0.+0.*icom;
		for(int i=0;i<N;i++) SS[i]=0.+0.*icom;
		S2=log(T)/pi/T;
        for(int iw=0; iw<w_obrez; iw++)
		{	w=pi*T*(2.*iw+1.); //cual of S2
			S2+=2./w;
	  	}

	  	std::vector<std::thread> threads(amount_of_threads);
        for (int i=0; i<amount_of_threads-1; i++)
        {   //w0=wth[i];w1=wth[i+1];
            threads[i] = std::thread(Scalc,std::ref(Del), std::ref(S1),std::ref(wth[i]),std::ref(wth[i+1]), std::ref(q), std::ref(I));
        }

        for (int i=0; i<amount_of_threads-1; i++)
          threads[i].join();

		/*thread th1(Scalc,std::ref(Del), std::ref(S1),std::ref(w0 ),std::ref(w1));
	    thread th2(Scalc,std::ref(Del), std::ref(S1),std::ref(w1 ),std::ref(w2));
	    thread th3(Scalc,std::ref(Del), std::ref(S1),std::ref(w2 ),std::ref(w3));
	    thread th4(Scalc,std::ref(Del), std::ref(S1),std::ref(w3 ),std::ref(w4));
	    thread th5(Scalc,std::ref(Del), std::ref(S1),std::ref(w4 ),std::ref(w5));
	    thread th6(Scalc,std::ref(Del), std::ref(S1),std::ref(w5 ),std::ref(w6));
	    thread th7(Scalc,std::ref(Del), std::ref(S1),std::ref(w6 ),std::ref(w7));
	    thread th8(Scalc,std::ref(Del), std::ref(S1),std::ref(w7 ),std::ref(w8));
	    thread th9(Scalc,std::ref(Del), std::ref(S1),std::ref(w8 ),std::ref(w9));
	    thread th10(Scalc,std::ref(Del),std::ref(S1),std::ref(w9 ),std::ref(w10));
	    thread th11(Scalc,std::ref(Del),std::ref(S1),std::ref(w10),std::ref(w11));
	    thread th12(Scalc,std::ref(Del),std::ref(S1),std::ref(w11),std::ref(w12));
	    thread th13(Scalc,std::ref(Del),std::ref(S1),std::ref(w12),std::ref(w13));
	    thread th14(Scalc,std::ref(Del),std::ref(S1),std::ref(w13),std::ref(w14));
	    thread th15(Scalc,std::ref(Del),std::ref(S1),std::ref(w14),std::ref(w15));
	    thread th16(Scalc,std::ref(Del),std::ref(S1),std::ref(w15),std::ref(w16));
	    thread th17(Scalc,std::ref(Del),std::ref(S1),std::ref(w16),std::ref(w17));
	    //thread th18(Scalc,std::ref(Del),std::ref(S1),std::ref(w17),std::ref(w18));//*/
	    //th1.join();th2.join();th3.join();th4.join();th5.join();th6.join();th7.join();th8.join();th9.join();th10.join();th11.join();th12.join();th13.join();th14.join();th15.join();th16.join();th17.join();//th18.join();

	    for(int i=0;i<N;i++) for(int iw=0; iw<w_obrez; iw++)   SS[i]+=S1[iw][i];//

        dDelmax = 0; int imax;
	        for (int i=0; i<N; i++)
	        {
	            Delbuf= Del[i];
	            Del[i]= (SS[i]/S2+Delbuf*(alpha-1.))/alpha;
                if (dDelmax< abs(Del[i]-Delbuf))
	            {   dDelmax = abs(Del[i]-Delbuf);
	                imax=i;
	            }
	        }

            cout << iter << " " << T << "  " << dDelmax << "  " << real(Del[0]) << "  " << real(Del[N -N_N -1]) << endl;// +N_F + N_S1 + N_F1 + N_S2 / 2]) << "  " << real(Del[N - N_N - 1]) << endl;
	        iter++;   // try it until it converges
	        //for(int i=0;i<N;i++)    cout<<i<<"  "<<Del[i]<<endl;
	        //auto end = std::chrono::steady_clock::now();
            //auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
            //cout <<"t("<<iter<<")=" << elapsed_ms.count() << endl;
    }

delete [] SS;
for (int i = 0; i < w_obrez; i++)
    delete[] S1[i];
delete[] S1;

}


void Scalc(complex<double> *Del, complex<double> **S1, int &w_in, int &w_fin, double* q, double I)  //working in each thread
{
	int iterG=0, w_count=w_fin-w_in;
    double ww, dGmax;
    complex<double> *G,*G1, *Fi, *Fi_old, *Fi1, *Fi1_old, Gbuf,wm;
    G=new complex<double>[N];
    G1 = new complex<double>[N];
    Fi=new complex<double>[N];
    Fi_old=new complex<double>[N];
    Fi1=new complex<double>[N];
    Fi1_old=new complex<double>[N];

	for(int iw=w_in; iw<w_fin; iw++)
    {
        for (int i=0; i<N; i++)	    G[i]=1.+0.1*icom;
		//cout<<"iw="<<iw<<"  ID="<<this_thread::get_id()<<endl;
        ww=pi*T*(2.*iw+1.);
        dGmax=1.;
        iterG = 0; ;// int kk = 0; double dDmax_last = 0;
			 // Iterative loop over GG. Solve Nonlinear Usadel equation here
            while ((dGmax > epsG)&&(iterG<1000))//&&(kk==0))
            {   // Calculation of Anomalous Green functions Fi(w), Fi(-w)
                //for (int i = 0; i < N; i++) G1[i] = -conj(G[i]);
                //for (int i = 0; i < N; i++) Fi_old[i] = Fi[i];
                Prog( Fi, G, Del, ww, q, I);
                //for (int i = 0; i < N; i++) Fi[i] = 0.1*Fi[i] + 0.9*Fi_old[i];
                for (int i = 0; i < N; i++) Fi1[i] = conj(Fi[i]);
				//Prog( Fi1, G1, Del, -ww, q, I);
                Gcalc(G, &dGmax, Fi, Fi1, Del, ww,q);
                iterG++;
                //if(iw==1) cout<<Fi1[N/2]<<endl;
                //if(iw==0) cout<<"w="<<((ww/pi/T-1)/2.)<<"  "<<iterG<<"  "<< dGmax<<endl;
			}
            //if(iw==0)
                //cout<<"w="<<((ww/pi/T-1)/2.)<<"  "<<iterG<<"  "<< real(Del[0])<<"  "<<dGmax<<endl;//"  arg="<<arg(G[N-1]*Fi[N-1]/get_wm(N-1,ww))/pi<<endl;
            // search summ of self-cons eq. for each point of the grid
            for (int i=0; i<N; i++)
			{
                //ww = pi * T * (2. * iw + 1.);// +2. * pi * get_ksi(i) * get_ksi(i) * q[i] * G[i] * q[i] / 2.;
                wm = get_wm(i, ww) +2. * pi * get_ksi(i) * get_ksi(i) * G[i] * q[i] * q[i] / 2.;
                S1[iw][i] += 2. * get_tc(i) * get_type(i) * real(Fi1[i] / sqrt(wm * wm + Fi1[i] * conj(Fi[i])));// G[i] * Fi[i] / wm);// +conj(G[i]) * Fi1[i] / ww);
                //fout <<fixed<< (w/pi/T-1)/2 <<"\t"<< i <<"\t" << real(F[i].A[0]) << "\t" << imag(F[i].A[0]) << "\t" << real(G[i]) << "\t" << imag(G[i]) << "\t" <<real(Del[i])<<endl;// log10(abs(1.-G[i]*G[i]))<<  endl;
                //if(iw==0) fout1 <<fixed<< (w/pi/T-1)/2 <<"\t"<< i <<"\t" << real(F[i].A[1]) << "\t" << imag(F[i].A[1]) <<"\t" << real(F[i].A[2]) << "\t" << imag(F[i].A[2]) <<"\t" << real(F[i].A[3]) << "\t" << imag(F[i].A[3]) <<endl;// log10(abs(1.-G[i]*G[i]))<<  endl;

			}//*/
        }
    delete [] G;
    delete [] G1;
    delete [] Fi;
    delete [] Fi1;
    delete[] Fi_old;
    delete[] Fi1_old;

}
//////// SelfConsZero - Pair Potentail in bulk material ////////////

double SelfConsZero()
{
	double Del1, dDel, w;
    double S1, S2;

    Del1=1.76;  // Initial Delta is BCS value at T = 0
    dDel=1;
    while (abs(dDel)>epsDel)
        {
        S1=0;
        S2=log(T)/pi/(T)/2.;
        for(int iw=0; iw<w_obrez; iw++)
            {
            w=pi*T*(2.*iw+1.);
            S1 = S1 + Del1/sqrt(w*w+Del1*Del1);   // Use known Green Functions Fi in bulk superconductor
            S2 +=1./w;
            }

        dDel=Del1 - S1/S2;
        Del1= Del1 - dDel;
        }
    return Del1;
}

	//				for(int j=0;j<N;j++)
	//				{
	//					cout<<"Fi_old["<<j<<"]="<<Fi_old[j]<<" "<<"Fi["<<j<<"]="<<Fi[j]<<" "<<"Fi1["<<j<<"]="<<Fi1_old[j]<<" "<<"T="<<T<<'\n';
	//					cout<<"iter="<<iter<<'\n';
		//				cout<<"Fi["<<j<<"]="<<Fi[j]<<'\n';
			//			cout<<"Fi1["<<j<<"]="<<Fi1_old[j]<<" ";
			//		}
