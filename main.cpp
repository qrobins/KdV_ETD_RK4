//
//  main.cpp
//  KdVETDRK4
//
//  Created by Quentin Robinson on 12/19/16.
//  Copyright (c) 2016 Quentin.Robinson. All rights reserved.
//  
//  Scheme from Kassam & Trefethen
//

#include <cmath>
#include <iostream>
#include <complex>
#include <iterator>
#include <vector>
#include <fstream>
#include <cstdio>
#include <ctime>

using namespace std;

typedef complex<double> cx;

void KCfft(cx [], cx []);
void KCifft(cx [], cx []);
void flux(double [], double, double []);
void KurgTadm(double, double [], cx [], double, double[]);
void minmod(double [], double[], double []);

const int N = 16384;
const int m = log2(N);
const int nplots = 50;
const cx J = cx(0,1);

int main(int argc, const char * argv[])
{
    //std::clock_t start;
    //double duration;
    
    //start = clock();
    
    // spatial domain info
    double xmax = 800.0;   // domain length
    double dx = xmax/N;     // spatial step size
    double x[N];            // spatial domain
    // loop to assign values to spatial domain
    for (int i=0; i<N; ++i)
        x[i] = (i+1)*dx - xmax/2;
    
    // frequency domain info
    cx k[N];                // frequency domain
    // loops to assign values to frequency domain
    for (int i=0; i <= N/2; i++)
        k[i] = i*(2*M_PI/xmax);
    for (int i=-N/2+1; i <= -1; i++)
        k[i+N] = i*(2*M_PI/xmax);
    
    // time info
    double t = 0.0;       // initial time value
    double tmax = 100.0;   // final time value
    double dt = 0.001;   // time step size
    
    // KdV equation parameters
    double alpha = 0.66667;   // nonlinearity
    double eps = 1.0; //epsilonIn.0 * 0.001;     // dispersion
    double Fr = 2.0;      // Froude number

    // Forcing parameters
    double pm = 0.5;                // forcing amplitude
    double xi = 0.3;
    //double omega = omegaIn.0 * 0.01;             // forcing oscillation frequency
    //double osc_amp = pow(0.99,2);    // forcing oscillation amplitude
    // declare variables for surface elevation and forcing with their Fourier transforms
    cx u[N], u_hat[N], forcing[N], f_hat[N];
    cx u_lin[N], u_lin_hat[N];
    double u_burg[N];
    
    /*
    // exact solution initial condition
    for (int i=0; i<N; i++)
        u[i] = 4.0/3.0*pow(xi,2)*pow(1/cosh(xi*x[i]),2);
    // other initial condition for testing
    // alternate initial condition
    for (int i=0; i<N; i++){
        u[i] = .5*pow(1/cosh(2.*sqrt(3*alpha/(4*eps))*x[i]),2);
        u_burg[i] = real(u[i]);
        u_lin[i] = u[i];}
    */
     
    // fft of initial condition
    KCfft(u, u_hat); KCfft(u_lin, u_lin_hat);
    
    
    // loop to assign values to derivative of forcing function f
    
    for (int i=0; i<N; i++)
        forcing[i] = -2.0*pm*xi*pow(1/cosh(xi*x[i]),2)*tanh(xi*x[i]);
    /*
    // forcing function with known exact solution
    for (int i=0; i<N; i++)
        forcing[i] = -16.0/9.0*pow(xi,5)*pow(1/cosh(xi*x[i]),2)*tanh(xi*x[i]);
    */

    KCfft(forcing, f_hat); // fft of derivative of forcing function
    
    cout << "alpha = " << alpha << ", epsilon = " << eps << ", Froude = " << Fr << ", dt = " << dt;
    cout << ", pm = " << pm << ", xi = " << xi; //<< ", osc_amp = " << osc_amp << ", omega = " << omega;
    cout << ", f_x(x) = -2.0*pm*xi*pow(1/cosh(xi*x[i]),2)*tanh(xi*x[i])\n\n";
    
    cx L[N], E[N], E2[N];
    for (int i=0; i<N; i++)
    {
        L[i] = J*k[i]*(1-Fr) - eps/6*J*pow(k[i],3);
        E[i] = exp(dt*L[i]);
        E2[i] = exp(dt*L[i]/2.);
    }
    
    const int M = 32;         // no. of points for complex means
    cx r[M], LR[N][M], mean, f1[N], f2[N], f3[N];
    double Q[N];
    
    for (int i=0; i<M; i++)
    {
        r[i] = exp(J*M_PI*((double)(i+1.)-.5)/(double)M);    // roots of unity
        
        for (int j=0; j<N; j++)
            LR[j][i] = dt*conj(L[j]) + r[i];
    }
    
    for (int j=0; j<N; j++)
    {
        mean = 0;
        for (int ii=0; ii<M; ii++)
            mean += (exp(LR[j][ii]/2.0)-1.0)/LR[j][ii];
        
        mean = mean/(double)M;
        
        Q[j] = dt * real(mean);
        
        mean = 0;
        for (int ii=0; ii<M; ii++)
            mean += (-4.-LR[j][ii]+exp(LR[j][ii])*(4.-3.*LR[j][ii]+pow(LR[j][ii],2)))/pow(LR[j][ii],3)/(double)M;
        f1[j] = dt * real(mean);
        
        mean = 0;
        for (int ii=0; ii<M; ii++)
            mean += (2.+LR[j][ii]+exp(LR[j][ii])*(-2.+LR[j][ii]))/pow(LR[j][ii],3)/(double)M;
        f2[j] = dt * real(mean);
        
        mean = 0;
        for (int ii=0; ii<M; ii++)
            mean += (-4.-3.*LR[j][ii]-pow(LR[j][ii],2)+exp(LR[j][ii])*(4.-LR[j][ii]))/pow(LR[j][ii],3)/(double)M;
        f3[j] = dt * real(mean);
    }
    
    // data for y axis of waterfall plot
    double tdata[nplots+1];
    // declare variables to calculate and store result of full and linearized KdV
    cx uu[nplots+1][N], g[N], non[N], a[N], b[N], c[N], Na[N], Nb[N], Nc[N];
    cx uu_lin[nplots+1][N];
    // declare variables to calculate and store result of dispersionless case
    double uu_burg[nplots+1][N], eff1[N], eff2[N], eff3[N], eff4[N], place_hold[N];
    
    for (int i=0; i<N; i++)
        g[i] = 3.0/4.0*alpha*J*k[i];
    // store initial condition in results array
    for (int i=0; i<N; i++){
        uu[0][i] = u[i];
        uu_lin[0][i] = u_lin[i];
        uu_burg[0][i] = u_burg[i];
    }
    
    double filter = floor(N/3);
    if (fmod(filter,2)==0)
        filter = filter - 1;

    // Main time stepping loop
    for (int i=0; i<nplots; i++)
    {
        for (int ii=0; ii<round(tmax/dt/nplots); ii++)
        {
            /*
            // for time dependent forcing
            for (int j=0; j<N; j++)
                forcing[j] = -2.0*pm*xi*(1.0-osc_amp*sin(omega*t))*pow(1.0/cosh(xi*x[j]),2)*tanh(xi*x[j]);
            KCfft(forcing,f_hat);
            */
            t = t + dt;     // update time variable
            
            
            //Filter
            for (int j=N/2+1-(filter-1)/2-1; j<=N/2+(filter-1)/2; j++)
                u_hat[j] = 0;
            
            KCifft(u_hat,u);
            
            
            for (int j=0; j<N; j++)
                u[j] = pow(real(u[j]),2);
            KCfft(u, non);
            for (int j=0; j<N; j++)
                non[j] = g[j] * non[j] + f_hat[j];
            
            for (int j=0; j<N; j++)
                a[j] = E2[j] * u_hat[j] + Q[j] * non[j];
            KCifft(a,u);
            for (int j=0; j<N; j++)
                u[j] = pow(real(u[j]),2);
            KCfft(u,Na);
            for (int j=0; j<N; j++)
                Na[j] = g[j] * Na[j] + f_hat[j];
     
            for (int j=0; j<N; j++)
                b[j] = E2[j]*u_hat[j] + Q[j]*Na[j];
            KCifft(b,u);
            for (int j=0; j<N; j++)
                u[j] = pow(real(u[j]),2);
            KCfft(u,Nb);
            for (int j=0; j<N; j++)
                Nb[j] = g[j]*Nb[j] + f_hat[j];
     
            for (int j=0; j<N; j++)
                c[j] = E2[j]*a[j] + Q[j]*(2.0*Nb[j]-non[j]);
            KCifft(c,u);
            for (int j=0; j<N; j++)
                u[j] = pow(real(u[j]),2);
            KCfft(u,Nc);
            for (int j=0; j<N; j++)
                Nc[j] = g[j]*Nc[j] + f_hat[j];
     
            for (int j=0; j<N; j++)
                u_hat[j] = E[j]*u_hat[j] + non[j]*f1[j] + 2.0*(Na[j]+Nb[j])*f2[j] + Nc[j]*f3[j];
            
            KCifft(u_hat,u);
            
            
            
            // Solving the linear problem (i.e. g = 0)
            // filter
            for (int j=N/2+1-(filter-1)/2-1; j<=N/2+(filter-1)/2; j++)
                u_lin_hat[j] = 0;
            KCifft(u_lin_hat,u_lin);
            
            for (int j=0; j<N; j++){
                u_lin[j] = pow(real(u_lin[j]),2);
                non[j] = f_hat[j];
            }
            
            for (int j=0; j<N; j++)
                a[j] = E2[j] * u_lin_hat[j] + Q[j] * non[j];
            KCifft(a,u_lin);
            for (int j=0; j<N; j++){
                u_lin[j] = pow(real(u_lin[j]),2);
                Na[j] = f_hat[j];
            }
            
            for (int j=0; j<N; j++)
                b[j] = E2[j]*u_lin_hat[j] + Q[j]*Na[j];
            KCifft(b,u_lin);
            for (int j=0; j<N; j++){
                u_lin[j] = pow(real(u_lin[j]),2);
                Nb[j] = f_hat[j];
            }
            
            for (int j=0; j<N; j++)
                c[j] = E2[j]*a[j] + Q[j]*(2.0*Nb[j]-non[j]);
            KCifft(c,u_lin);
            for (int j=0; j<N; j++){
                u_lin[j] = pow(real(u_lin[j]),2);
                Nc[j] = f_hat[j];
            }
            
            for (int j=0; j<N; j++)
                u_lin_hat[j] = E[j]*u_lin_hat[j] + non[j]*f1[j] + 2.0*(Na[j]+Nb[j])*f2[j] + Nc[j]*f3[j];
            
            KCifft(u_lin_hat,u_lin);
            
            
            
            // Solving the dispersionless problem
            KurgTadm(dx,u_burg,forcing,Fr,eff1); 
            for (int j=0; j<N; j++)
                place_hold[j] = u_burg[j] + dt/2.0 * eff1[j];
            KurgTadm(dx,place_hold,forcing,Fr,eff2);
            for (int j=0; j<N; j++)
                place_hold[j] = u_burg[j] + dt/2.0 * eff2[j];
            KurgTadm(dx,place_hold,forcing,Fr,eff3);
            for (int j=0; j<N; j++)
                place_hold[j] = u_burg[j] + dt*eff3[j];
            KurgTadm(dx,place_hold,forcing,Fr,eff4);
            for (int j=0; j<N; j++)
                u_burg[j] = u_burg[j] + dt/6.0 * (eff1[j] + 2.0*eff2[j] + 2.0*eff3[j] + eff4[j]);
        }
        
        tdata[i+1] = t;
        for (int j=0; j<N; j++){
            uu[i+1][j] = u[j];
            uu_lin[i+1][j] = u_lin[j];
            uu_burg[i+1][j] = u_burg[j];
        }
    }
    

    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    //std::cout<<"printf: "<< duration <<'\n';
    /*
    // Output result
    cout << 0 << ' ';
    for (int j=0; j<N; j++)
        cout << x[j] << ' ';
    cout << '\n';
    for(int j=0; j<nplots+1; j++){
        cout << tdata[j] << ' ';
        for(int i = 0; i < N; i++)
            cout << real(uu[j][i]) << ' ' ;
        cout << '\n';
    }
    for(int j=0; j<nplots+1; j++){
        cout << tdata[j] << ' ';
        for(int i=0; i<N; i++)
            cout << real(uu_lin[j][i]) << ' ';
        cout << '\n';
    }
    for(int j=0; j<nplots+1; j++){
        cout << tdata[j] << ' ';
        for(int i=0; i<N; i++)
            cout << uu_burg[j][i] << ' ';
        cout << '\n';
    }
    */
    // Write result to file
    
    ofstream myfile ("surf_data.txt");
    if (myfile.is_open())
    {
        myfile << 0 << ' ';
        for (int j=0; j<N; j++)
            myfile << x[j] << ' ';
        myfile << '\n';
        for(int j=0; j<nplots+1; j++){
            myfile << tdata[j] << ' ';
            for(int i = 0; i < N; i++)
                myfile << real(uu[j][i]) << ' ' ;
            myfile << '\n';
        }
        for(int j=0; j<nplots+1; j++){
            myfile << tdata[j] << ' ';
            for(int i=0; i<N; i++)
                myfile << real(uu_lin[j][i]) << ' ';
            myfile << '\n';
        }
        for(int j=0; j<nplots+1; j++){
            myfile << tdata[j] << ' ';
            for(int i=0; i<N; i++)
                myfile << uu_burg[j][i] << ' ';
            myfile << '\n';
        }
        myfile.close();
    }
    else cout << "Unable to open file";
    
    
    return 0;
}





//Kincaid & Cheney fft
void KCfft(cx zeta[], cx zeta_hat[])
{
    cx u, v, w = exp(-2*M_PI*J/ (double)N);
    cx Z[N], D[N];
    
    for (int k=0; k<N; k++)
    {
        Z[k] = pow(w,k);
        zeta_hat[k] = zeta[k];
    }
    
    for (int n=0; n<m; n++)
    {
        for (int k=0; k<pow(2,m-n-1); k++)
        {
            for (int j=0; j<pow(2,n); j++)
            {
                u = zeta_hat[(int)pow(2,n)*k+j];
                v = Z[j*(int)pow(2,m-n-1)]*zeta_hat[(int)pow(2,n)*k+(int)pow(2,m-1)+j];
                D[(int)pow(2,n+1)*k+j] = (u+v)/2.;
                D[(int)pow(2,n+1)*k+j+(int)pow(2,n)] = (u-v)/2.;
            }
        }
        for (int j=0; j<N; j++)
            zeta_hat[j] = D[j];
    }
    
    for(int n=0; n<N; n++)
        zeta_hat[n] *= N;
}





//Inverse fft
void KCifft(cx zeta_hat[], cx zeta[])
{
    cx u, v, w = exp(2*M_PI*J/ (double)N);
    cx Z[N], D[N];
    
    for (int k=0; k<N; k++)
    {
        Z[k] = pow(w,k);
        zeta[k] = zeta_hat[k];
    }
    
    for (int n=0; n<m; n++)
    {
        for (int k=0; k<pow(2,m-n-1); k++)
        {
            for (int j=0; j<pow(2,n); j++)
            {
                u = zeta[(int)pow(2,n)*k+j];
                v = Z[j*(int)pow(2,m-n-1)]*zeta[(int)pow(2,n)*k+(int)pow(2,m-1)+j];
                D[(int)pow(2,n+1)*k+j] = (u+v)/2.;
                D[(int)pow(2,n+1)*k+j+(int)pow(2,n)] = (u-v)/2.;
            }
        }
        for (int j=0; j<N; j++)
            zeta[j] = D[j];
    }
}





// Kurganov Tadmor scheme for numerically capturing propagating shocks
void KurgTadm(double dx, double u[], cx force[], double fr, double kt[])
{
    
    double a[N], b[N], c[N], u_x[N], u_plus[N], u_minus[N], x;
    
    a[0] = (u[0] - u[N-1])/dx;
    b[N-1] = (u[0] - u[N-1])/dx;
    for(int i=1; i<N; i++){
        a[i] = (u[i] - u[i-1])/dx;
        b[i-1] = (u[i]-u[i-1])/dx;
    }
    
    minmod(a,b,u_x);
    
    for (int i=0; i<N; i++){
        u_plus[i] = u[i] - dx/2.0*u_x[i];
        u_minus[i] = u[i]+dx/2.0*u_x[i];
    }
    
    x = u_minus[N-1];
    for (int i=N-1; i>0; i--)
        u_minus[i] = u_minus[i-1];
    u_minus[0] = x;
    
    //spectral radius of the Jacobian of F(u)???
    for (int i=0; i<N; i++)
        a[i] = !(abs(u_minus[i])>abs(u_plus[i]))?abs(u_plus[i]):abs(u_minus[i]);
    
    b[N-1] = u_plus[0];
    for (int i=0; i<N-1; i++)
        b[i] = u_plus[i+1];
    
    flux(b,fr,u_x);
    
    b[N-1] = u_minus[0];
    for (int i=0; i<N-1; i++)
        b[i] = u_minus[i+1];
    
    flux(b,fr,c);
    
    for (int i=0; i<N; i++)
        u_x[i] = u_x[i] + c[i];
    
    flux(u_plus,fr,c);
    
    for (int i=0; i<N; i++)
        u_x[i] = u_x[i] - c[i];
    
    flux(u_minus,fr,c);
    
    for (int i=0; i<N; i++){
        u_x[i] = u_x[i] - c[i];
        u_x[i] = -u_x[i]/(2.0*dx);
    }
    
    b[N-1] = u_minus[0];
    c[N-1] = u_plus[0];
    b[N-1] = c[N-1] - b[N-1];
    for (int i=0; i<N-1; i++){
        b[i] = u_minus[i+1];
        c[i] = u_plus[i+1];
        b[i] = c[i] - b[i];
    }
    
    c[N-1] = a[0];
    b[N-1] = b[N-1]*c[N-1];
    for (int i=0; i<N-1; i++){
        c[i] = a[i+1];
        b[i] = b[i]*c[i];
    }
    
    for (int i=0; i<N; i++){
        a[i] = a[i]*(u_plus[i] - u_minus[i]);
        b[i] = (b[i] - a[i])/(2.0*dx);
        kt[i] = u_x[i] + b[i] + real(force[i]);
    }
}





// Minmod
void minmod(double a[], double b[], double u[])
{
    double c, d;
    
    for (int i=0; i<N; i++){
        if (a[i] == 0)
            c = 0;
        else
            c = a[i]/abs(a[i]);
        if (b[i] == 0)
            d = 0;
        else
            d = b[i]/abs(b[i]);
        u[i] = ((c)+(d))/2.0*(!(abs(b[i])<abs(a[i]))?abs(a[i]):abs(b[i]));
    }
}





// Flux
void flux(double f[], double fr, double y[])
{
    for (int i=0; i<N; i++)
        y[i] = (fr - 1)*f[i] - 3.0*pow(f[i],2);
}
