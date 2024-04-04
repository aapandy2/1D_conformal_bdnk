/* Simple solver for the conformal BDNK equations in slab symmetry.
 * 
 * Written by Alex Pandya, last updated 4/4/2024; numerical method
 * described in:
 * [1] https://arxiv.org/abs/2201.12317
 *
 * Choose options for the simulation in this file.
* */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <parameters.h>

/*number of time levels stored in array; keep this at 3 */
#define TL (3)

double dx;
double dt;

double x[N];
double xi[TL][N];
double ux[TL][N];

double eps[TL][N];

double xiD[TL][N];
double uxD[TL][N];

double xiP[TL][N];
double uxP[TL][N];

double Ttt[TL][N];
double Ttx[TL][N];

double VISC[TL][N];

double res[TL][N];

double flux_tx[N];
double flux_xx[N];

double Ttt_pv[N];
double Ttx_pv[N];

double T_tt0(double xi, double ux, double xicx, double uxcx, double xidot, 
             double uxdot);
double T_tx0(double xi, double ux, double xicx, double uxcx, double xidot, 
             double uxdot);
double T_xx0(double xi, double ux, double xicx, double uxcx, double xidot, 
             double uxdot);
double T_tt(double xi, double ux, double xicx, double uxcx, double xidot, 
            double uxdot);
double T_tx(double xi, double ux, double xicx, double uxcx, double xidot, 
            double uxdot);
double T_xx(double xi, double ux, double xicx, double uxcx, double xidot, 
            double uxdot);

int save_x_to_file();
int save_t_to_file();

/*should return qL at x=i+1/2*/
double WENO_reconst_qRx(double q[TL][N], int n, int i)
{
    i += 1;

    double v0 = 11./6.*q[n][i]   - 7./6.*q[n][i+1] + 1./3.*q[n][i+2];
    double v1 =  1./3.*q[n][i-1] + 5./6.*q[n][i]   - 1./6.*q[n][i+1];
    double v2 = -1./6.*q[n][i-2] + 5./6.*q[n][i-1] + 1./3.*q[n][i];

    double d0 = 1./10.;
    double d1 = 3./5.;
    double d2 = 3./10.;

    double beta0 = (1./4.  *pow(3.*q[n][i]-4.*q[n][i+1]+q[n][i+2],2.)
                   +13./12.*pow(   q[n][i]-2.*q[n][i+1]+q[n][i+2],2.)
                   );
    double beta1 = (1./4.  *pow(q[n][i+1]-q[n][i-1],2.)
                   +13./12.*pow(q[n][i-1]-2.*q[n][i]+q[n][i+1],2.)
                   );
    double beta2 = (1./4.  *pow(q[n][i-2]-4.*q[n][i-1]+3.*q[n][i],2.)
                   +13./12.*pow(q[n][i-2]-2.*q[n][i-1]+q[n][i],2.)
                   );
    double wt0 = d0/pow(epsW + beta0,2.);
    double wt1 = d1/pow(epsW + beta1,2.);
    double wt2 = d2/pow(epsW + beta2,2.);
    double w0  = wt0/(wt0+wt1+wt2);
    double w1  = wt1/(wt0+wt1+wt2);
    double w2  = wt2/(wt0+wt1+wt2);

    double qL = w0*v0 + w1*v1 + w2*v2;
    
    return qL;
}

/*should return qR at x=i+1/2 */
double WENO_reconst_qLx(double q[TL][N], int n, int i)
{
    double v0 =  1./3.*q[n][i]   +5./6.*q[n][i+1] -1./6. *q[n][i+2];
    double v1 = -1./6.*q[n][i-1] +5./6.*q[n][i]   +1./3. *q[n][i+1];
    double v2 =  1./3.*q[n][i-2] -7./6.*q[n][i-1] +11./6.*q[n][i];

    double d0 = 3./10.;
    double d1 = 3./5.;
    double d2 = 1./10.;

    double beta0 = (1./4.  *pow(3.*q[n][i]-4.*q[n][i+1]+q[n][i+2],2.)
                   +13./12.*pow(   q[n][i]-2.*q[n][i+1]+q[n][i+2],2.)
                   );
    double beta1 = (1./4.  *pow(q[n][i+1]-q[n][i-1],2.)
                   +13./12.*pow(q[n][i-1]-2.*q[n][i]+q[n][i+1],2.)
                   );
    double beta2 = (1./4.  *pow(q[n][i-2]-4.*q[n][i-1]+3.*q[n][i],2.)
                   +13./12.*pow(q[n][i-2]-2.*q[n][i-1]+q[n][i],2.)
                   );
    double wt0 = d0/pow(epsW + beta0,2.);
    double wt1 = d1/pow(epsW + beta1,2.);
    double wt2 = d2/pow(epsW + beta2,2.);
    double w0  = wt0/(wt0+wt1+wt2);
    double w1  = wt1/(wt0+wt1+wt2);
    double w2  = wt2/(wt0+wt1+wt2);

    double qL = w0*v0 + w1*v1 + w2*v2;

    return qL;
}

double reconst_qLx(double q[TL][N], int n, int i)
{
    return WENO_reconst_qLx(q,n,i);
}

double reconst_qRx(double q[TL][N], int n, int i)
{
    return WENO_reconst_qRx(q,n,i);
}

/*centered WENO FD; should be fourth order in smooth regions */
double wDx(double f[3][N], int n, int i)
{
    double eps_W = epsW;

    double beta1 = (1./4.  *pow(f[n][i-2]-4.*f[n][i-1]+3.*f[n][i],2.)
                   +13./12.*pow(f[n][i-2]-2.*f[n][i-1]+f[n][i],2.)
                   );
    double beta2 = (1./4.  *pow(f[n][i+1]-f[n][i-1],2.)
                   +13./12.*pow(f[n][i-1]-2.*f[n][i]+f[n][i+1],2.)
                   );
    double beta3 = (1./4.  *pow(3.*f[n][i]-4.*f[n][i+1]+f[n][i+2],2.)
                   +13./12.*pow(   f[n][i]-2.*f[n][i+1]+f[n][i+2],2.)
                   );

    double C1 = 1./6.;
    double C2 = 2./3.;
    double C3 = 1./6.;

    double alpha1 = C1/pow(eps_W + beta1, 2.);
    double alpha2 = C2/pow(eps_W + beta2, 2.);
    double alpha3 = C3/pow(eps_W + beta3, 2.);

    double w1 = alpha1/(alpha1+alpha2+alpha3);
    double w2 = alpha2/(alpha1+alpha2+alpha3);
    double w3 = alpha3/(alpha1+alpha2+alpha3);

    double Rp = (w1*(f[n][i-2]-4.*f[n][i-1]+3.*f[n][i])/(2.*dx)
                +w2*(f[n][i+1]-f[n][i-1])/(2.*dx)
                +w3*(-3.*f[n][i]+4.*f[n][i+1]-f[n][i+2])/(2.*dx)
                );

    return Rp;
}

double Dx(double arr[3][N], int n, int i)
{
    if(i > 1 && i < N-2)
    {
        return wDx(arr,n,i);
    }

    return (arr[n][i+1]-arr[n][i-1])/(2.*dx);
}


int set_initial_data()
{
    dx = (X_MAX-X_MIN)/(N-1.);
    dt = CFL * dx;

    /*define epsL, epsR, vL, vR for steady-state shocks */
    double epsL, epsR, vL, vR, vTemp, eps;
    for(int i=0; i<N; i++)
    {
        x[i] = X_MIN + i*dx;

#if ID_TYPE == GAUSSIAN
        /*Gaussian initial data.  Note the coefficient
         * and the constant added---without these, velocities
         * get too close to c and NANs appear*/ 
        eps      = GAUSSIAN_AMPLITUDE
                    *exp(-pow(x[i]-GAUSSIAN_MEAN,2.)
                            /GAUSSIAN_SPREAD) 
                    + GAUSSIAN_CONST;
        ux[0][i] = 0.;

#elif ID_TYPE == STEP
        if(x[i] <= 0.)
        {
            eps = STEP_EPS_L;
        }
        else
        {
            eps = STEP_EPS_R;
        }
        ux[0][i] = 0.;

#elif ID_TYPE == SMOOTH_SHOCK
        /*use PF jump conditions to get a shock in its rest frame */
        epsL      = SMOOTH_SHOCK_EPS_L;
        vL        = SMOOTH_SHOCK_V_L;
        epsR      = (epsL-9.*vL*vL*epsL)/(3.*(-1.+vL*vL));
        vR        = 1./(3.*vL);
        eps       = (epsR-epsL)/2.*(erf(x[i]/10.)+1.)+epsL;
        vTemp     = (vL-vR)/2.*(1.-erf(x[i]/10.))+vR;
        ux[0][i]  = vTemp/sqrt(1.-vTemp*vTemp);
#endif
       
        //definition of xi is ln(eps) 
        xi[0][i] = log(eps);

        /*set dissipative tensors to zero (perfect fluid initial data) */
        Ttt[0][i] = T_tt(xi[0][i], ux[0][i], 0., 0., 0., 0.);
        Ttx[0][i] = T_tx(xi[0][i], ux[0][i], 0., 0., 0., 0.);
    }

    for(int i=0; i<N; i++)
    {
        xiD[0][i] = 0.;
        uxD[0][i] = 0.;

        xiP[0][i] = Dx(xi,0,i);
        uxP[0][i] = Dx(ux,0,i);
    }

    //save x to a file called x.txt
    save_x_to_file();

    //save t to a file called t.txt
    save_t_to_file();

    return 0;
}

int set_gc(double arr[TL][N], int n)
{
#if BC == GHOST
    arr[n][0]    = arr[n][3];
    arr[n][1]    = arr[n][3];
    arr[n][2]    = arr[n][3];
    arr[n][N-3] = arr[n][N-4];
    arr[n][N-2] = arr[n][N-4];
    arr[n][N-1] = arr[n][N-4];
#elif BC == PERIODIC    
    arr[n][0]    = arr[n][N-6];
    arr[n][1]    = arr[n][N-5];
    arr[n][2]    = arr[n][N-4];
    arr[n][N-3] = arr[n][3];
    arr[n][N-2] = arr[n][4];
    arr[n][N-1] = arr[n][5];
#endif
    return 0;
}

int set_ghost_cells(int n)
{
    set_gc(xi,n);
    set_gc(ux,n);
    set_gc(xiP,n);
    set_gc(uxP,n);
    set_gc(xiD, n);
    set_gc(uxD, n);
    set_gc(Ttx,n);
    set_gc(Ttt,n);
    
    return 0;
}

double compute_A(double XI, double UX, double xiP, double uxP, double xiD, 
                 double uxD)
{
    return ((chi0*pow(exp(XI),0.75)*(4*UX*uxD + 4*sqrt(1 + pow(UX,2))*uxP +
    3*(1 + pow(UX,2))*xiD + 3*UX*sqrt(1 + pow(UX,2))*xiP))/(4.*sqrt(1 +
    pow(UX,2))));
}

double compute_Qx(double XI, double UX, double xiP, double uxP, double xiD, 
                  double uxD)
{
    return ((pow(exp(XI),0.75)*lambda0*(4*sqrt(1 + pow(UX,2))*uxD + 4*UX*uxP +
    UX*sqrt(1 + pow(UX,2))*xiD + (1 + pow(UX,2))*xiP))/4.);
}

double compute_m2sxx(double XI, double UX, double xiP, double uxP, double xiD, 
                     double uxD)
{
    return ((-4*pow(exp(XI),0.75)*eta0*sqrt(1 + pow(UX,2))*(UX*uxD + sqrt(1 +
    pow(UX,2))*uxP))/3.);
}

double T_tt(double xi, double ux, double xicx, double uxcx, double xidot, 
            double uxdot)
{
    double A     = compute_A    (xi, ux, xicx, uxcx, xidot, uxdot);
    double Qx    = compute_Qx   (xi, ux, xicx, uxcx, xidot, uxdot);
    double m2sxx = compute_m2sxx(xi, ux, xicx, uxcx, xidot, uxdot);

    double ut           = sqrt(1.+ux*ux);
    double Deltatt      = -1. + ut*ut;
    double Qt           = ux*Qx/ut;
    double m2etasigmatx = ux*m2sxx/ut;
    double m2etasigmatt = ux*m2etasigmatx/ut;
    double eps          = exp(xi);

    double ans = eps*(ut*ut+Deltatt/3.) + A*(ut*ut+Deltatt/3.) + 2.*Qt*ut 
                 + m2etasigmatt;

    return ans;
}

double T_tt0(double xi, double ux, double xicx, double uxcx, double xidot, 
             double uxdot)
{
    double ut           = sqrt(1.+ux*ux);
    double Deltatt      = -1. + ut*ut;
    double eps          = exp(xi);

    double ans = eps*(ut*ut+Deltatt/3.);

    return ans;
}

double T_tx(double xi, double ux, double xicx, double uxcx, double xidot, 
            double uxdot)
{
    double A     = compute_A    (xi, ux, xicx, uxcx, xidot, uxdot);
    double Qx    = compute_Qx   (xi, ux, xicx, uxcx, xidot, uxdot);
    double m2sxx = compute_m2sxx(xi, ux, xicx, uxcx, xidot, uxdot);

    double ut           = sqrt(1.+ux*ux);
    double Deltatx      = ut*ux;
    double Qt           = ux*Qx/ut;
    double m2etasigmatx = ux*m2sxx/ut;
    double m2etasigmatt = ux*m2etasigmatx/ut;
    double eps          = exp(xi);

    double ans = eps*(ut*ux+Deltatx/3.) + A*(ut*ux+Deltatx/3.) +Qt*ux + ut*Qx 
                 + m2etasigmatx;

    return ans;
}

double T_xx(double xi, double ux, double xicx, double uxcx, double xidot, 
            double uxdot)
{
    double A     = compute_A    (xi, ux, xicx, uxcx, xidot, uxdot);
    double Qx    = compute_Qx   (xi, ux, xicx, uxcx, xidot, uxdot);
    double m2sxx = compute_m2sxx(xi, ux, xicx, uxcx, xidot, uxdot);

    double ut           = sqrt(1.+ux*ux);
    double Deltaxx      = 1. + ux*ux;
    double Qt           = ux*Qx/ut;
    double m2etasigmatx = ux*m2sxx/ut;
    double m2etasigmatt = ux*m2etasigmatx/ut;
    double eps          = exp(xi);

    double ans = eps*(ux*ux+Deltaxx/3.) + A*(ux*ux+Deltaxx/3.)+2.*Qx*ux+m2sxx;

    return ans;
}

double T_tx0(double xi, double ux, double xicx, double uxcx, double xidot, 
             double uxdot)
{
    double ut           = sqrt(1.+ux*ux);
    double Deltatx      = ut*ux;
    double eps          = exp(xi);

    double ans = eps*(ut*ux+Deltatx/3.);

    return ans;
}

double T_xx0(double xi, double ux, double xicx, double uxcx, double xidot, 
             double uxdot)
{
    double ut           = sqrt(1.+ux*ux);
    double Deltaxx      = 1. + ux*ux;
    double eps          = exp(xi);

    double ans = eps*(ux*ux+Deltaxx/3.);

    return ans;
}

double compute_xiD(int n, int i, double XI, double UX, double xiP, double uxP, 
                   double T00, double T01)
{
    double ch = chi0/eta0;
    double l  = lambda0/eta0;

    double TttPF = ((9*exp(XI) + 8*exp(XI)*pow(UX,4) + 6*pow(UX,2)*(3*exp(XI) -
    2*pow(exp(XI),0.75)*eta0*uxP) + 3*pow(exp(XI),0.75)*eta0*pow(UX,3)*xiP)/(9
    + 6*pow(UX,2)));

    double TtxPF =  ((UX*sqrt(1 + pow(UX,2))*(8*exp(XI)*pow(UX,2) + 12*(exp(XI)
    - pow(exp(XI),0.75)*eta0*uxP) + 3*pow(exp(XI),0.75)*eta0*UX*xiP))/(9 +
    6*pow(UX,2)));

    double dtt, dtx;
    if(eta0 != 0)
    {
        dtt = (TttPF - T00)/eta0;
        dtx = (TtxPF - T01)/eta0;
    }
    else /*catch NANs */
    {
        dtt = 0.;
        dtx = 0.;
    }

    double eps;
    if(eta0*fabs(dtt) < TOL || eta0*fabs(dtx) < TOL)
    {
        dtt = 0.;
        dtx = 0.;
        eps = (-T00+sqrt(4.*pow(T00,2.)-3.*pow(T01,2.)));
        xi[n][i] = log(eps);
        ux[n][i] = 3.*T01/sqrt(pow(3.*T00+eps,2.)-pow(3.*T01,2.));

        //track that we're using PF solution
        VISC[0][i] = 0;
    }
    else
    {
        VISC[0][i] = 1;
    }

    return ((-2*(2*uxP + UX*(1 + pow(UX,2))*xiP))/(sqrt(1 + pow(UX,2))*(3 +
    2*pow(UX,2))) - (4*sqrt(1 + pow(UX,2))*(3*l + (-4 + 4*ch +
    6*l)*pow(UX,2))*dtt)/(pow(exp(XI),0.75)*(9*ch*l + 12*ch*(-1 + l)*pow(UX,2)
    + 4*(ch*(-3 + l) - l)*pow(UX,4))) + (4*(3*(ch + 2*l)*UX + (-4 + 4*ch +
    6*l)*pow(UX,3))*dtx)/(pow(exp(XI),0.75)*(9*ch*l + 12*ch*(-1 + l)*pow(UX,2)
    + 4*(ch*(-3 + l) - l)*pow(UX,4))));
}

double compute_uxD(int n, int i, double XI, double UX, double xiP, double uxP, 
                   double T00, double T01)
{
    double ch = chi0/eta0;
    double l  = lambda0/eta0;
    
    double TttPF = ((9*exp(XI) + 8*exp(XI)*pow(UX,4) + 6*pow(UX,2)*(3*exp(XI) - 
    2*pow(exp(XI),0.75)*eta0*uxP) + 3*pow(exp(XI),0.75)*eta0*pow(UX,3)*xiP)/(9
    + 6*pow(UX,2)));

    double TtxPF =  ((UX*sqrt(1 + pow(UX,2))*(8*exp(XI)*pow(UX,2) + 12*(exp(XI) 
    - pow(exp(XI),0.75)*eta0*uxP) + 3*pow(exp(XI),0.75)*eta0*UX*xiP))/(9 +
    6*pow(UX,2)));

    double dtt, dtx;
    if(eta0 != 0)
    {
        dtt = (TttPF - T00)/eta0;
        dtx = (TtxPF - T01)/eta0;
    }
    else /*catch NANs */
    {
        dtt = 0.;
        dtx = 0.;
    }

    double eps;
    if(eta0*fabs(dtt) < TOL || eta0*fabs(dtx) < TOL)
    {
        dtt = 0.;
        dtx = 0.;
        eps = (-T00+sqrt(4.*pow(T00,2.)-3.*pow(T01,2.)));
        xi[n][i] = log(eps);
        ux[n][i] = 3.*T01/sqrt(pow(3.*T00+eps,2.)-pow(3.*T01,2.));
    }

    return (-((sqrt(1 + pow(UX,2))*(8*UX*uxP + 3*xiP))/(12 + 8*pow(UX,2))) +
    (3*UX*sqrt(1 + pow(UX,2))*(4*ch + l + 2*(2*ch +
    l)*pow(UX,2))*dtt)/(pow(exp(XI),0.75)*(9*ch*l + 12*ch*(-1 + l)*pow(UX,2) +
    4*(ch*(-3 + l) - l)*pow(UX,4))) - (3*(1 + pow(UX,2))*(3*ch + 2*(2*ch +
    l)*pow(UX,2))*dtx)/(pow(exp(XI),0.75)*(9*ch*l + 12*ch*(-1 + l)*pow(UX,2) +
    4*(ch*(-3 + l) - l)*pow(UX,4))));
}

/*conservative centered pi flux at i+1/2 */
double flux_x(int n, int i, int COMP)
{
    double xiL = reconst_qLx(xi,n,i);
    double xiR = reconst_qRx(xi,n,i);
    double uxL = reconst_qLx(ux,n,i);
    double uxR = reconst_qRx(ux,n,i);

    double xiDL = reconst_qLx(xiD,n,i);
    double xiDR = reconst_qRx(xiD,n,i);
    double uxDL = reconst_qLx(uxD,n,i);
    double uxDR = reconst_qRx(uxD,n,i);

    /*note: eps', ux' are computed as arrays, and
            then reconstructed in pi flux */
    double xicxL = reconst_qLx(xiP,n,i);
    double uxcxL = reconst_qLx(uxP,n,i);
    double xicxR = reconst_qRx(xiP,n,i);
    double uxcxR = reconst_qRx(uxP,n,i);

    double iphL, iphR, jump_iph;
    if(COMP == 0)
    {
        iphL = T_tx(xiL, uxL, xicxL, uxcxL, xiDL, uxDL);
        iphR = T_tx(xiR, uxR, xicxR, uxcxR, xiDR, uxDR);
        jump_iph = T_tt(xiR, uxR, xicxR, uxcxR, xiDR, uxDR) 
                   - T_tt(xiL, uxL, xicxL, uxcxL, xiDL, uxDL);
    }
    else
    {
        iphL = T_xx(xiL, uxL, xicxL, uxcxL, xiDL, uxDL); 
        iphR = T_xx(xiR, uxR, xicxR, uxcxR, xiDR, uxDR);
        jump_iph = T_tx(xiR, uxR, xicxR, uxcxR, xiDR, uxDR) 
                   - T_tx(xiL, uxL, xicxL, uxcxL, xiDL, uxDL);
    }

    /*define KT speed (speed of light) and jump at i+1/2 */
    double a        = 1.;
    double flux_iph = 0.5*(iphL + iphR - a*jump_iph);

    return flux_iph;
}

double Ttx_cx(int n, int i)
{
    double flux_iph = flux_x(n, i,   0);
    double flux_imh = flux_x(n, i-1, 0);

    double ans = (flux_iph-flux_imh)/dx;

    return ans;
}

double Txx_cx(int n, int i)
{
    double flux_iph = flux_x(n, i,   1);
    double flux_imh = flux_x(n, i-1, 1);

    double ans = (flux_iph-flux_imh)/dx;

    return ans;
}

int Heun_solve_system(int n)
{
    set_ghost_cells(n);
    double xicx, uxcx;
    for(int i=3; i<N-3; i++)
    {
        /*evolve Ttt, Ttx, Tty, Jt (predictor step) */
        Ttt[n+1][i] = (Ttt[n][i] + dt*( -Ttx_cx(n,i)
                                      ));
        Ttx[n+1][i] = (Ttx[n][i] + dt*( -Txx_cx(n,i)
                                      ));

#if PRIM == BDNK_PRIM
        xi[n+1][i] = xi[n][i] + dt*xiD[n][i];
        ux[n+1][i] = ux[n][i] + dt*uxD[n][i];
#elif PRIM == PF_PRIM
        eps[n+1][i] = (-Ttt[n+1][i]
                      +sqrt(4.*pow(Ttt[n+1][i],2.)-3.*pow(Ttx[n+1][i],2.))
                      );
        ux[n+1][i]  = 3.*Ttx[n+1][i]/sqrt(pow(3.*Ttt[n+1][i]+eps[n+1][i],2.)
                                         -pow(3.*Ttx[n+1][i],2.));
#endif
    }

    set_ghost_cells(n+1);
    for(int i=3; i<N-3; i++)
    {
        //compute spatial derivative arrays
        xiP[n+1][i] = Dx(xi,n+1,i);
        uxP[n+1][i] = Dx(ux,n+1,i);

        xicx = Dx(xi,n+1,i);
        uxcx = Dx(ux,n+1,i);

        xiD[n+1][i] = compute_xiD(n+1, i, xi[n+1][i], ux[n+1][i], xicx, uxcx, 
                                  Ttt[n+1][i], Ttx[n+1][i]);
        uxD[n+1][i] = compute_uxD(n+1, i, xi[n+1][i], ux[n+1][i], xicx, uxcx,
                                  Ttt[n+1][i], Ttx[n+1][i]);
    }

    set_ghost_cells(n+1);
    for(int i=3; i<N-3; i++)
    {
        /*evolve Ttt, Ttx, Tty, Jt (corrector step) */
        Ttt[n+2][i] = (Ttt[n][i] + dt/2.*(-Ttx_cx(n,i)
                                          -Ttx_cx(n+1,i)
                                         ));
        Ttx[n+2][i] = (Ttx[n][i] + dt/2.*(-Txx_cx(n,i)
                                          -Txx_cx(n+1,i)
                                         ));

#if PRIM == BDNK_PRIM
        xi[n+2][i] = xi[n][i] + dt/2.*(xiD[n][i]+xiD[n+1][i]);
        ux[n+2][i] = ux[n][i] + dt/2.*(uxD[n][i]+uxD[n+1][i]);
#elif PRIM == PF_PRIM
        eps[n+2][i] = (-Ttt[n+2][i]
                      +sqrt(4.*pow(Ttt[n+2][i],2.)-3.*pow(Ttx[n+2][i],2.))
                      );
        ux[n+2][i]  = 3.*Ttx[n+2][i]/sqrt(pow(3.*Ttt[n+2][i]+eps[n+2][i],2.)
                                         -pow(3.*Ttx[n+2][i],2.));
#endif
    }

    set_ghost_cells(n+2);
    for(int i=3; i<N-3; i++)
    {
        xicx = Dx(xi,n+2,i);
        uxcx  = Dx(ux,n+2,i);

        //compute spatial derivative arrays
        xiP[n+2][i] = Dx(xi,n+2,i);
        uxP[n+2][i] = Dx(ux,n+2,i);

        xiD[n+2][i] = compute_xiD(n+2, i, xi[n+2][i], ux[n+2][i], xicx, uxcx,
                                  Ttt[n+2][i], Ttx[n+2][i]);
        uxD[n+2][i] = compute_uxD(n+2, i, xi[n+2][i], ux[n+2][i], xicx, uxcx,
                                  Ttt[n+2][i], Ttx[n+2][i]);
    }
    set_ghost_cells(n+2);

    return 0;
}

int copy_forward(int n)
{
    /*copy n+1 timelevel (stored in TL-1 array position) 
      into n level (0 array position)*/
    for(int i=0; i<N; i++)
    {
        xi[0][i] = xi[TL-1][i];
        ux[0][i] = ux[TL-1][i];

        xiD[0][i] = xiD[TL-1][i];
        uxD[0][i] = uxD[TL-1][i];

        xiP[0][i] = xiP[TL-1][i];
        uxP[0][i] = uxP[TL-1][i];

        Ttt[0][i] = Ttt[TL-1][i];
        Ttx[0][i] = Ttx[TL-1][i];
    }
    return 0;
}

int compute_eps(int n)
{
    for(int i=0; i<N; i++)
    {
        eps[n][i] = exp(xi[n][i]);
    }
    return 0;
}

int output_arr(double arr[TL][N], char* filename_str, int n)
{
    int GRID_SAMPLE_RATE_X = 1;

    char file_loc[100];
    strcpy(file_loc, DIREC);

    /*output arr to file*/
    strcat(file_loc, filename_str);
    FILE *farr = fopen(file_loc, "a");

    for(int i=0; i<N; i+=GRID_SAMPLE_RATE_X)
    {
        fprintf(farr, "%e  ", arr[n][i]);
    }
    fprintf(farr, "\n");
    fclose(farr);
    strcpy(file_loc, DIREC);

    return 0;
}

int save_x_to_file()
{
    int GRID_SAMPLE_RATE_X = 1;

    char file_loc[100];
    strcpy(file_loc, DIREC);

    /*output arr to file*/
    strcat(file_loc, "/x.txt");
    FILE *farr = fopen(file_loc, "a");

    for(int i=0; i<N; i+=GRID_SAMPLE_RATE_X)
    {
        fprintf(farr, "%e  ", x[i]);
    }
    fprintf(farr, "\n");
    fclose(farr);
    strcpy(file_loc, DIREC);

    return 0;
}

int save_t_to_file()
{
    char file_loc[100];
    strcpy(file_loc, DIREC);

    /*output arr to file*/
    strcat(file_loc, "/t.txt");
    FILE *farr = fopen(file_loc, "a");

    for(int n=0; n<MAX_TIMESTEP; n++)
    {
        fprintf(farr, "%e\n", n*dt);
    }
    fprintf(farr, "\n");
    fclose(farr);
    strcpy(file_loc, DIREC);

    return 0;
}

int print_to_file(int n)
{

    /*output variables */
    output_arr(xi,  "/xi.txt", n);
    output_arr(ux,   "/ux.txt",  n);
    output_arr(Ttt, "/Ttt.txt", n);
    output_arr(Ttx, "/Ttx.txt", n);
    output_arr(xiD,  "/xiD.txt", n);
    output_arr(uxD,   "/uxD.txt",  n);
    output_arr(VISC, "/VISC.txt", n);

    /*output residuals */
    compute_eps(n);
    output_arr(eps, "/eps.txt", n);

    return 0;
}

int main()
{
    set_initial_data();

    for(int n=0; n<MAX_TIMESTEP; n++)
    {
        Heun_solve_system(0);

        if(n % TS_STEP == 0)
        {
            print_to_file(0);
            printf("\nTimestep: %d", n);
        }

        copy_forward(0);
    }

    printf("\n");
    return 0;
}
