#include<iostream>
#include<fstream>
#include<sstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <vector>
#include <string>
using namespace std;

#define ZERO            1.e-10
#define DOANTICIPATE    true
#define FACPI		0

char dividends = 'A';

//simulation length
//#define	T               1000
//#define	Teq             3000
#define Tshock          20
//#define	tprint          1

//choices
#define RENORM          1           // 1-renormalize to inflation / 0 - no normalization

//demand
//#define G0              0.5         // fraction of savings

//revival
//#define phi             0.1         // revivial probability per unit of time
//#define taupi           0.2

//others
//#define delta           0.02        // dividends share

#define doshock false   //put false here

//GSL stuff
const gsl_rng_type *gslT ;
gsl_rng *gslr ;
gsl_rng *gslri ;

double min(double x, double y){ if(x<y) return x; else return y; }
double max(double x, double y){ if(x>y) return x; else return y; }
double red(double x){ return ( fabs(x) < ZERO ? 0.0 : x);}

int main(int argc, char *argv[ ] ){
    
    FILE    *out; //, *out2;
    char    add[10000], name[10000]; //  name2[10000];
    
    //variables
    double  bust, Pavg, Pold, u, S, Atot, firm_savings, debt_tot, Ytot, Wtot, e, Wavg, inflation, k, propensity, Dtot, rho;
    double  rhom, rhop, rp_avg, pi_avg, Gamma, u_avg, rm_avg;

    //simulation parameters
    int     N, T, Teq, tprint;
    double y0;
    
    //parameters
    double  gammap, gammaw, theta, Gamma0, eta0m, eta0p, beta, R, r, alpha, rho0, f, alpha_e, alpha_pi, Gammas;
    double G0, phi, taupi, delta;
    
    //others;
    double  Pmin, rp, rw, ren, Pnorm, arg, pay_roll, dY, p , tmp, budget, interests;
    double  Wmax, wage_norm, u_share, deftot;
    int     i, t, seed, seed_init, new_firm, new_len;
    double  profits;
    int    negative_count=0;
    double tau_meas , tau_tar;
    
    double pi_target, e_target;
    double wage_factor;
    
    
    if(argc!=31){
        printf("Incorrect input structure\n");
        exit(1);
    }
    
    sscanf(argv[1],  "%lf", &R); // Firm: Hiring/firing rate (2017-12) (default 2.0)
    sscanf(argv[2],  "%lf", &Gamma0); // Firm: loan rate effect on hire/fire (2017-13) (default 50)
    sscanf(argv[3],  "%lf", &Gammas); // Firm: loan rate effect on hire/fire (2017-13) (default 0.0)
    sscanf(argv[4],  "%lf", &r); // Firm: baseline gamma (2017-13) (default 1.0)
    sscanf(argv[5],  "%lf", &gammap); // Firm: adjustment ratio gamma_w / gamma_p (2015-18) (default: 0.01)
    sscanf(argv[6],  "%lf", &eta0m); // Firm: baseline firing propensity (2017-12) (default: 0.2)
    sscanf(argv[7],  "%lf", &wage_factor); // Factor to adjust wages to inflation expectations (default: 1.0)
    sscanf(argv[8],  "%lf", &rho0); // CB: baseline interest rate (2017-3) (default 0.02)
    sscanf(argv[9],  "%lf", &alpha_pi); // CB: reaction to inflation (2017-3) (default: 2.5)
    sscanf(argv[10], "%lf", &alpha_e); // CB: reaction to employment (2017-3) (default: 0.0)
    sscanf(argv[11], "%lf", &pi_target); // CB: inflation target (2017-3) (default 0.04)
    sscanf(argv[12], "%lf", &e_target); // CB: unemployment target (2017-3) (default 0.01)
    sscanf(argv[13], "%lf", &theta); // Bank: Default threshold (2015-16) (default 3.0)
    sscanf(argv[14], "%lf", &f); // Bank: bankruptcy effect on bank interest rates (2017-7) (default 0.5)
    sscanf(argv[15], "%lf", &alpha); // HH: real rate influence on consumption (2017-9) (default 0.0)
    sscanf(argv[16], "%lf", &G0); // HH: Baseline propensity to consume (2017-9) (default 0.5)
    sscanf(argv[17], "%lf", &delta); // HH: Dividend share (2017-9) (default 0.02)
    sscanf(argv[18], "%lf", &beta); // HH:  Intensity of choice (household demand) (default 0.5)
    sscanf(argv[19], "%lf", &tau_meas); // HH: expect pi weight of EWMA inflation (2018-7) (default 0.5)
    sscanf(argv[20], "%lf", &tau_tar); // HH: expect pi weight of CB target infl.(2018-7) (default 0.5)
    sscanf(argv[21], "%lf", &phi); // Economy: Revival probability per unit time (default 0.1)
    sscanf(argv[22], "%lf", &taupi); // Economy: Memory kernel for EWMA estimates (2017-4) (default 0.2)
    sscanf(argv[23], "%lf", &y0); // Economy: initial production (default 0.5)
    sscanf(argv[24], "%d",  &seed); // Simulation seed --> NOTE: there is no seed_init here
    sscanf(argv[25], "%d",  &seed_init); // Simulation seed --> NOTE: there is no seed_init here
    sscanf(argv[26], "%d",  &N); // Number of firms
    sscanf(argv[27], "%d",  &T); // Total runtime
    sscanf(argv[28], "%d",  &Teq); // Cutoff time
    sscanf(argv[29], "%d",  &tprint); // Log every nth period
    sscanf(argv[30], "%s",  add); // filename


    string s_add(add);
    string cons_chk("cons_shock");
    string theta_chk("theta_inf");
    
    //array
    vector  <double>  P(N), Y(N), D(N), A(N), W(N), PROFITS(N);
    vector  <int>     ALIVE(N), new_list(N);
    
    double R0 = R;
    double G = G0;
    double ao = 1.0;
    double M0 = ao*N;
    eta0p  = R*eta0m;
    gammaw = r*gammap;

    //printf("%.2e, %.2e, %.2e, %.2e , %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e , %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %d, %d, %d, %d, %d \n ", R, Gamma0, Gammas, r, gammap, eta0m, wage_factor, rho0, alpha_pi, alpha_e, pi_target, e_target, theta, f, alpha, G0, delta, beta, tau_meas, tau_tar, phi, taupi, y0, seed, N, T, Teq, tprint);
    /*
    printf("R = %.2e\nN = %d\ntheta = %.2e\ngp = %.2e\tgw = %.2e\nf = %.2e\nb = %.2e\nG = %.2e\nalpha=%f\n\n",R,N,theta,gammap,gammaw,f,beta,Gamma0,alpha);
    printf("rho0 = %.2e\tap = %.2e\tae = %.2e\n",rho0,alpha_pi,alpha_e);
    printf("pit = %.2e\tet = %.2e\n",pi_target,e_target);
    printf("taut = %.2e\ttaum = %.2e\n",tau_tar,tau_meas);
    printf("seed %d\n",seed);
    printf(" %.2e, %.2e, %.2e, %.2e , %.2e, %.2e, %.2e, %.2e, %.2e, %d, %d, %d, \n", Gammas, r, eta0m, wage_factor, G0, delta, phi, taupi, y0, T, Teq, tprint );
    */
    
    char params[10000];
    sprintf(params,"%.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e",rho0,alpha_pi,alpha_e,pi_target,e_target,theta,R,Gamma0,Gammas,alpha,tau_meas,tau_tar);
    
    //final averages
    double  avg_u, avg_bu, avg_k, avg_pi, avg_rho, avg_rhom, avg_rhop;
    double  sig_u, sig_bu, sig_k, sig_pi;
    double  max_u, min_u, max_pi, min_pi, max_rho, min_rho, max_rhom, min_rhom, max_rhop, min_rhop;
    
    int     collapse, avg_counter;
    
    double theta0 = theta;
    /* ************************************** INIT ************************************** */
    
    if(seed==-1)
        seed = 0; // time(NULL);
    
    //gsl_rng_env_setup() ;
    //gslri = gsl_rng_default ;
    gslri = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(gslri,seed_init);
    
    Pavg = 0.;
    Ytot = 0.;
    Wtot = 0.;
    Atot = 0.;
    Wmax = Wavg = 1.;
    Pold = 1.;
    inflation = pi_avg = 0.;
    rhom = rho = rm_avg = rho0;
    rhop = rp_avg = 0.;
    u_avg = 0.;
    
    for(i=0;i<N;i++){
        
        ALIVE[i] =  1;
        P[i] 	 =  1.  + 0.01*(2*gsl_rng_uniform(gslri)-1.);
        Y[i] 	 =  y0  + 0.01*(2*gsl_rng_uniform(gslri)-1.);
        D[i] 	 =  y0;
        W[i] 	 =  1.;
        PROFITS[i] = P[i]*min(D[i],Y[i]) - W[i]*Y[i];

        A[i] =  2.0*Y[i]*W[i]*gsl_rng_uniform(gslri);
        //printf("%lf, %lf \n", P[i], Y[i]);
        //printf("%lf", gsl_rng_uniform(gslr));
        Atot += A[i];
        Ytot += Y[i];
        Pavg += P[i]*Y[i];
        Wtot += Y[i]*W[i];
        
    }
    
    e = Ytot / N;
    u = 1. - e ;
    
    Pavg /= Ytot;
    Wavg = Wtot/Ytot;
    
    S = e*N;
    //S = 0.0;//modificato
    
    //fix total amount of money to N
    tmp = Atot + S;
    S = S*N/tmp*ao;
    Atot=0.;
    
    for(i=0;i<N;i++){
        A[i]  = A[i]*N/tmp*ao;
        Atot += A[i];
        
    }

    gslr = gsl_rng_alloc(gsl_rng_default) ;;
    gsl_rng_set(gslr, seed);
    //cout << "eccoci " << S << " " << Atot << endl;
    
    /* *********************************** MAIN CYCLE ************************************ */
    avg_counter = collapse = 0;
    avg_bu = avg_k = avg_pi = avg_u = 0.;
    avg_rho = avg_rhom = avg_rhop = 0.;
    sig_bu = sig_k = sig_pi = sig_u = 0.;
    
    max_rho = max_rhom = max_rhop = -1.0;
    min_rho = min_rhom = min_rhop = +1.0;
    max_pi = -1.0;
    min_pi = +1.0;
    
    
    max_u = 0.;
    min_u = 2.;
    
    bust=0.;
    
    if( alpha_pi < ZERO )
    {
        pi_target = 0.0;
        tau_tar   = 0.0;
    }
    
    sprintf(name,"%s.txt",add);
    out = (FILE*)fopen(name,"w");
    //sprintf(name2,"%s_firms.txt",add);
    //out2 = (FILE*)fopen(name2,"w");
    
    //cout << G0 << endl;
    
    for(t=0;t<T;t++){
        
        //printf("%d \n", t);
        //renormalize in unit of price
        
        
        if(RENORM==1){
            
            for(i=0;i<N;i++){
                P[i]/=Pavg;
                W[i]/=Pavg;
                A[i]/=Pavg;
                PROFITS[i] /= Pavg;
            }
            
            S    /= Pavg;
            Wavg /= Pavg;
            Wmax /= Pavg;
            Pold /= Pavg;
            M0 /= Pavg;
            Pavg  = 1.;
            
        }
        
        /* *********************************** UPDATE ************************************ */
        
        pi_avg = taupi*inflation + (1.-taupi)*pi_avg;
        rp_avg = taupi*rhop + (1.-taupi)*rp_avg;
        rm_avg = taupi*rhom + (1.-taupi)*rm_avg;
        u_avg = taupi*u + (1.-taupi)*u_avg;
        
        //update firms variables
        Wtot = 0.;
        Ytot = 0.;
        tmp  = 0.;
        Pmin = 1.e+300 ;
        
        if(beta>0.){
            wage_norm = 0.;
            for(i=0;i<N;i++)if(ALIVE[i]==1){
                arg = beta*(W[i]-Wmax)/Wavg ;
                if(arg>-100.) wage_norm += gsl_sf_exp(arg);
            }
        }
        
        new_len = 0 ;
        deftot = 0.;
        firm_savings = 0.;
        debt_tot = 0.;
        
        double pi_used = tau_tar * pi_target + tau_meas * pi_avg ;
        
        Gamma = max(Gamma0 * (rm_avg-pi_used),Gammas);
        for(i=0;i<N;i++){
            
            // living firms
            if(ALIVE[i]==1){
                
                pay_roll = Y[i]*W[i] ;
                
                // if not bankrupt update price / production / wages and compute interests / savings and debt
                
                if((A[i]>-theta*pay_roll)||(theta<0.)){
                    
                    
                    if(pay_roll>0.)
                        ren = Gamma * A[i] / pay_roll;
                    else ren = 0.;
                    
                    // this is not in the current code
                    //if (ren> 1.) ren=  1. ;
                    //if (ren<-1.) ren= -1. ;
                    
                    rp = gammap*gsl_rng_uniform(gslr);
                    rw = gammaw*gsl_rng_uniform(gslr);
                    
                    dY = D[i] - Y[i] ;
                    p  = P[i];
                    
                    if(beta>0.){
                        arg = beta*(W[i]-Wmax)/Wavg ;
                        u_share = 0.;
                        if(arg>-100.)u_share = u * N * (1.-bust) * gsl_sf_exp(arg) / wage_norm;
                    }
                    
                    else{
                        u_share = u;
                    }
                    
                    //excess demand
                    if(dY>0.){
                        
                        //increase production
                        double eta = eta0p*(1.+ren);
                        if(eta<0.0)eta=0.0;
                        if(eta>1.0)eta=1.0;
                        
                        Y[i] += min(eta*dY,u_share);
                        
                        //increase price
                        if(p<Pavg)
                            P[i] *= (1. + rp) ;
                        
                        //increase wage
                        if((PROFITS[i]>0.)&&(gammaw>ZERO)){
                            W[i] *= 1. + (1.0+ren) * rw * e;
                            W[i] = min(W[i],(P[i]*min(D[i],Y[i])+rhom*min(A[i],0.)+rhop*max(A[i],0.))/Y[i]);
                            W[i] = max(W[i],0.);
                        }
                        
                    }
                    
                    //excess production
                    else {
                        
                        //decrease production
                        double eta = eta0m*(1.-ren);
                        if(eta<0.0)eta=0.0;
                        if(eta>1.0)eta=1.0;
                        Y[i] += eta*dY ;
                        
                        //decrease price
                        if(p>Pavg)
                            P[i] *= (1. - rp) ;
                        
                        //decrease wage
                        if(PROFITS[i]<0.)
                        {
                            W[i] *= 1. - (1.0-ren)*rw*u;
                            W[i] = max(W[i],0.);
                            
                        }
                        
                    }
                    
                    if(DOANTICIPATE)
                    {
                        P[i] *= 1.0 + pi_used;
                        W[i] *= 1.0 + wage_factor * pi_used;
                    }
                    
                    Y[i] = max(Y[i],0.);
                    
                    Wtot += W[i]*Y[i];
                    tmp  += P[i]*Y[i];
                    Ytot += Y[i];
                    
                    firm_savings += max(A[i],0.);
                    debt_tot	 -= min(A[i],0.);
                    
                    Pmin = min(Pmin,P[i]);
                    
                    if((P[i]>1.0/ZERO)||(P[i]<ZERO)){
                        printf("price under/overflow... (1)\n");
                        if(P[i]>1.0/ZERO)collapse=3;
                        if(P[i]<ZERO)collapse=4;
                        t=T;
                    }
                }
                
                // if bankrupt shut down and compute default costs
                else {
                    deftot -= A[i];
                    Y[i] = 0.;
                    ALIVE[i] = 0;
                    A[i] = 0.;
                    new_list[new_len]=i;
                    new_len++;
                }
            }
            
            else{
                new_list[new_len]=i;
                new_len++;
            }
        }
        
        Pavg = tmp / Ytot ;
        Wavg = Wtot / Ytot ;
        e = Ytot / N ;
        u = 1. - e ;
        
        /* *********************************** INTERESTS ************************************ */
        
        // solve rounding errors
        if (fabs(S + firm_savings - deftot - debt_tot - M0) > 0.0)
        {
            double left = S + firm_savings - deftot - debt_tot - M0;
            //if too big, stop
            if (fabs(left) > ZERO*10 )
            {
                // printf("Accounting big problem! : tot_save = %.12e, tot_debt = %.12e, deftot = %.12e, M0 = %.12e, sum = %.12e\n", S+firm_savings, debt_tot, deftot, M0, left);
                // printf("Period of occurence: %d\n", t);
                exit(1);
            }
            //else remove from savings
            else
            {
                S -= left;
            }
            //check savings
            if (S < 0.0)
            {
                cout << "negative savings in interests" << endl;
                exit(1);
            }
        }

        rhom = rho;
        if(debt_tot>0.)rhom += (1.-f)*deftot / debt_tot ;
        
        //if(rhom>99) printf("%i, %lf \t", t, rhom);
        //rhom = min(rhom, 100);

        interests = rhom*debt_tot ;

        rhop = k = 0.;
        if( S + firm_savings > 0. ){
            rhop = (interests - deftot) / ( S + firm_savings );
            k = debt_tot / ( S + firm_savings) ;
        }

        S += rhop * S ;

        /* *********************************** CONSUMPTION ************************************ */
        
        propensity = G * (1.+ alpha*(pi_used-rp_avg) ) ;
        propensity = max(propensity,0.);
        propensity = min(propensity,1.);
        
//        if(t==Teq+1000)
//        {
//            S += 0.0001;
//            M0 += 0.0001;
//        }

        if(t > Teq &  t < Teq + Tshock & doshock)
        {
            //            for(i=0;i<N;i++){
            //                fprintf(out2,"%d %d %e %e %e %e %e\n",t,i,A[i],P[i],W[i],Y[i],D[i]);
            //            }
            if (s_add.find(cons_chk) != string::npos)
            {
                propensity *= 0.01;
                cout << t << " propensity " << propensity << " " << G << " " << alpha << endl;
            }
            if (s_add.find(theta_chk) != string::npos)
            {
                theta = theta0 * 100;
                cout << t << " theta " << theta << endl;
            }
            
        }
        else
        {
            theta = theta0;
        }
        
        budget = propensity * ( Wtot + max(S,0.) );
        
        Pnorm = 0.;
        for(i=0;i<N;i++)if(ALIVE[i]==1){
            arg = beta*(Pmin-P[i]) / Pavg ;
            if(arg>-100.) Pnorm += gsl_sf_exp(arg);
        }
        
        Dtot = 0.;
        profits = 0.;
        firm_savings = 0.;
        
        for(i=0;i<N;i++)if(ALIVE[i]==1){
            
            D[i] = 0.;
            
            arg = beta*(Pmin-P[i])/Pavg ;
            
            if(arg>-100.) D[i] = budget * gsl_sf_exp(arg) / Pnorm / P[i];
            
            PROFITS[i]  = P[i]*min(Y[i],D[i]) - Y[i]*W[i] + rhom*min(A[i],0.) + rhop*max(A[i],0.);
            
            S          -= P[i]*min(Y[i],D[i]) - Y[i]*W[i];
            A[i]       += PROFITS[i];
            
            if((A[i]>0.)&&(PROFITS[i]>0.)){
                if(dividends=='P'){
                    S    += delta*PROFITS[i];
                    A[i] -= delta*PROFITS[i];
                }
                if(dividends=='A'){
                    S    += delta*A[i];
                    A[i] -= delta*A[i];
                }
            }
            
            Dtot    += D[i];
            profits += PROFITS[i];
            
            firm_savings += max(A[i],0.);
            
        }
        
        /* ******************************* REVIVAL ******************************** */
        
        //revival
        deftot = 0.;
        
        for(i=0;i<new_len;i++)if(gsl_rng_uniform(gslr)<phi){
            
            new_firm = new_list[i];
            Y[new_firm] = max(u,0.)*gsl_rng_uniform(gslr);
            ALIVE[new_firm] = 1;
            P[new_firm] = Pavg;
            W[new_firm] = Wavg;
            
            A[new_firm] = W[new_firm]*Y[new_firm];
            deftot  += A[new_firm];
            
            firm_savings += A[new_firm];
            
            PROFITS[new_firm] = 0.;
            
        }
        
        /* ******************************* FINAL ******************************** */
        
        //new averages
        tmp  = 0.;
        Ytot = 0.;
        Wtot = 0.;
        bust = 0.;
        Wmax = 0.;
        Atot = 0.;
        
        debt_tot = 0.;
        
        for(i=0;i<N;i++){
            
            //final averages
            if(ALIVE[i]==1){
                
                if((firm_savings>0.)&&(A[i]>0.))A[i] -= deftot*A[i]/firm_savings;
                
                Wtot    += Y[i]*W[i];
                Ytot    += Y[i];
                tmp     += P[i]*Y[i];
                
                Wmax = max(W[i],Wmax);
                
                debt_tot -= min(A[i],0.);
                Atot     += A[i];
                
            }
            
            else bust += 1./N;
            
        }
        
        Pavg = tmp / Ytot;
        Wavg = Wtot / Ytot;
        
        inflation = (Pavg-Pold)/Pold;
        Pold = Pavg;
        
        e = Ytot / N  ;
        
        if((e-1 > ZERO)|| (e < -ZERO) || (S < ZERO) ){
            //printf("Error!! -> e = %.10e\tS=%.10e\n",e,S);
            collapse=2;
            t=T;
        }
        
        e = min(e,1.);
        e = max(e,0.);
        
        u =  1. - e;
        
        if(Ytot<ZERO){
            printf("Collapse\n");
            collapse = 1;
            t=T;
        }
        
        /************************************** CB ************************************/
        
        
        tmp = min((1.-u_avg)*1.025,e_target);
        
        rho = rho0 + FACPI * pi_target + alpha_pi * ( pi_avg - pi_target ) ;//+ alpha_e * gsl_sf_log((1.-u_avg)/tmp);
        
        /************************************** OUTPUT ************************************/
        
        
        if((t%tprint==0))
        {
            fprintf(out,"%d\t%e\t%e\t%e\t%e\t",t-Teq,u,bust,Pavg,Wavg);
            fprintf(out,"%.12e\t%.12e\t%.12e\t%.12e\t",S,Atot,firm_savings,debt_tot);
            fprintf(out,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",inflation,pi_avg,propensity,k,Dtot,rhom,rho,rhop,pi_used,tau_tar,tau_meas,R);
            fflush(out);
        }
        
        if(t>=Teq){
            
            if( rho < -ZERO )
                negative_count++;
            
            avg_u  += u;
            avg_pi += inflation;
            avg_rho += rho;
            avg_bu += bust;
            avg_k  += k;
            avg_rhom += rhom;
            avg_rhop += rhop;
            
            sig_u  += u*u;
            sig_pi += inflation*inflation;
            sig_bu += bust*bust;
            sig_k  += k*k;
            
            max_u = max(u,max_u);
            min_u = min(u,min_u);
            max_pi = max(inflation,max_pi);
            min_pi = min(inflation,min_pi);
            max_rho = max(rho,max_rho);
            min_rho = min(rho,min_rho);
            max_rhom =  max(rhom,max_rhom);
            min_rhom =  min(rhom,min_rhom);
            max_rhop =  max(rhop,max_rhop);
            min_rhop =  min(rhop,min_rhop);
            
            avg_counter++;
            
        }
        
    }
    
    if(avg_counter==0){
        max_u = 1.;
        min_u = 0.;
        avg_u = 1.;
        avg_pi = 0.;
        avg_bu = 1.;
        avg_k = 0.;
        sig_u = sig_bu = sig_k = sig_pi = 0.;
    }
    
    else{
        avg_k  /= avg_counter;
        avg_pi /= avg_counter;
        avg_rho /= avg_counter;
        avg_bu /= avg_counter;
        avg_u  /= avg_counter;
        avg_rhom /= avg_counter;
        avg_rhop  /= avg_counter;
        
        sig_k  /= avg_counter;
        sig_pi /= avg_counter;
        sig_bu /= avg_counter;
        sig_u  /= avg_counter;
        
        sig_k  -= avg_k*avg_k   ;
        sig_bu -= avg_bu*avg_bu ;
        sig_pi -= avg_pi*avg_pi ;
        sig_u  -= avg_u*avg_u   ;
        
    }
    
    fclose(out);
    //fclose(out2);
    
    double p_neg = double(negative_count)/avg_counter;
    
    char vars[10000];
    //sprintf(vars,"\nu: %.2e %.2e %.2e\npi: %.2e %.2e %.2e\noth: %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e",red(avg_u),red(max_u),red(min_u),red(avg_pi),red(max_pi),red(min_pi),
    //        red(avg_rho),red(max_rho),red(min_rho),red(avg_rhom),red(max_rhom),red(min_rhom),red(avg_rhop),red(max_rhop),red(min_rhop),red(p_neg));
    
    //printf("\n RESULTS:\n%s %s\n",params,vars);
    //printf("%.2e, %.2e, %.2e, %.2e , %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e , %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %.2e, %d, %d, %d, %d, %d \n ", R, Gamma0, Gammas, r, gammap, eta0m, wage_factor, rho0, alpha_pi, alpha_e, pi_target, e_target, theta, f, alpha, G0, delta, beta, tau_meas, tau_tar, phi, taupi, y0, seed, N, T, Teq, tprint);
    
}
