//
//  main.cpp
//  Basic implementation of comparison between SSO & PSO
//  Created by En-Rong Tsai on 2019/3/14.
//  Copyright © 2019 En-Rong Tsai. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define MaxNsol  300
#define MaxNvar  150

#define Nrun   30
#define Ngen   1000
#define Nsol   100
#define Nvar   100
//1 1 1
//1 2 2
//2 1 2
//2 2 1 orthogonality L4(2^3)
#define Cg    0.4 //1 1 2 2
#define Cp    0.8 //1 2 1 2
#define Cw    0.9 //1 2 2 1
//其中 1 = low level ; 2 = high level
#define Xub    100
#define Xlb   -100

clock_t         start                                      ;
unsigned int    run, gen                                   ;
// ---------------------------------------------------------
void            PARA_init(void)                            ;

float           fRAND_MAX                                  ;
float           fRAND_MAX1                                 ;
float           Xub_Xlb                                    ;

unsigned int    iCg                                        ;
unsigned int    iCp                                        ;
unsigned int    iCw                                        ;
// ---------------------------------------------------------
double          FIT_cal(float*XX)                          ;
// ---------------------------------------------------------
void            SSO_init(void)                             ;
void            SSO_update(void)                           ;

unsigned int    gBest, genBest                             ;
float           X [MaxNsol][MaxNvar]                       ;
float           P [MaxNsol][MaxNvar]                       ;
float           F [MaxNsol]                                ;
float           pF[MaxNsol]                                ;
//----------------------------------------------------------
void            PSO_init(void)                             ;
void            PSO_update(void)                           ;

float           V[MaxNsol][MaxNvar]                        ;
//----------------------------------------------------------
void            OUTPUT(void)                               ;
float           run_time                                   ;
//----------------------------------------------------------

int main(void)
{
    srand((unsigned) time (NULL))  ;
    
    PARA_init();
    for (run=0; run<Nrun; run++)
    {
        start=clock();
        printf("SSO ");
        SSO_init();
        for (gen=1; gen<Ngen; gen++)   SSO_update();
        OUTPUT();
        
        start=clock();
        printf("PSO ");
        PSO_init();
        for (gen=1; gen<Ngen; gen++)   PSO_update();
        OUTPUT();
        printf("\n")      ;
    }
}
/***************************************************************************/
void PARA_init(void)
{
    fRAND_MAX  = (float) RAND_MAX;
    fRAND_MAX1 = fRAND_MAX + 1   ;
    Xub_Xlb    = (Xub-Xlb)/fRAND_MAX1; //printf("Xub_Xlb=%f\n", Xub_Xlb);
    
    iCg = Cg * RAND_MAX;
    iCp = Cp * RAND_MAX;
    iCw = Cw * RAND_MAX;
}
/***************************************************************************/
void SSO_update(void)
{
    register int    sol               ;
    register int    var               ;
    unsigned int    rnd               ;
    
    for (sol=0; sol<Nsol; sol++)
    {
        for (var=0; var<Nvar; var++)
        {
            rnd = (double)rand();
            /*
             rnd = (double)rand()/(float)(RAND_MAX+1);
             rnd = (double)rand()/fRAND_MAX;
             rnd()<Cg
             ==> (double)rand()/fRAND_MAX<Cg
             ==> (double)rand()<Cg*fRAND_MAX
             ==> iCg = Cg*fRAND_MAX
             */
            if      (rnd<iCg) X[sol][var] = P[gBest][var];
            else if (rnd<iCp) X[sol][var] = P[sol]  [var];
            else if (rnd>iCw) X[sol][var] = Xub_Xlb * (double)rand() + Xlb;
            //               else if (rnd>iCw) X[sol][var] = (Xub-Xlb) * (double)rand()/fRAND_MAX1 + Xlb;
        }
        
        F[sol]=FIT_cal(X[sol]);
        if (F[sol]<pF[sol])
        {
            pF[sol]=F[sol];
            memcpy(P[sol], X[sol], Nvar * sizeof(float));
            if (F[sol]<pF[gBest])  { gBest=sol; genBest=gen; }
        }
    }
}
/***************************************************************************/
void PSO_update(void)
{
    register int    sol               ;
    register int    var               ;
    float           rnd1, rnd2, w     ;
    
    for (sol=0; sol<Nsol; sol++)
    {
        rnd1 = (double)rand()/(float)(RAND_MAX+1);
        rnd2 = (double)rand()/(float)(RAND_MAX+1);
        w    = 0.95*(1-(float)gen/(float)Ngen);
        for (var=0; var<Nvar; var++)
        {
            V[sol][var] = w * V[sol][var]
            + 2*(P[gBest][var]-X[sol][var]) * rnd1
            + 2*(P[sol][var]-X[sol][var])   * rnd2 ;
            
            if      (V[sol][var]> 1)  V[sol][var] =  1 ;
            else if (V[sol][var]<-1)  V[sol][var] = -1 ;

            X[sol][var] += V[sol][var];
            if      (X[sol][var]<Xlb)   X[sol][var] = Xlb ;
            else if (X[sol][var]>Xub)   X[sol][var] = Xub ;
        }
        
        F[sol]=FIT_cal(X[sol]);
        if (F[sol]<pF[sol])
        {
            pF[sol]=F[sol];
            memcpy(P[sol], X[sol], Nvar * sizeof(float));
            if (F[sol]<pF[gBest])  { gBest=sol; genBest=gen; }
        }
    }
}
/***************************************************************************/
void SSO_init(void)
{
    register int    sol, var   ;
    
    for (gBest=sol=0; sol<Nsol; sol++)
    {
        for (var=0; var<Nvar; var++)
            P[sol][var] = X[sol][var] = Xub_Xlb * (double)rand() + Xlb;
        //               P[sol][var] = X[sol][var] = X(ub-Xlb) * (double)rand()/(float)(RAND_MAX+1) + Xlb;
        
        pF[sol]=F[sol]=FIT_cal(X[sol]);
        if (F[sol]<F[gBest]) gBest=sol;
    }
}
/***************************************************************************/
void PSO_init(void)
{
    unsigned int    sol, var                  ;
    
    SSO_init();
    
    for (sol=0; sol<Nsol; sol++)
        for (var=0; var<Nvar; var++)
            V[sol][var]=(double)rand()/(float)(RAND_MAX+1) - 2 ;
}
/***************************************************************************/
double FIT_cal(float*XX)
{
    register int  var  ;
    double        SUM  ;
    
    for (SUM=var=0; var<Nvar; var++) SUM += XX[var]*XX[var] ;
    return SUM;
}
/***************************************************************************/
void OUTPUT(void)
{
    register int  var  ;
    
    printf("%d %f %d ", run, (float)(clock()-start)/CLOCKS_PER_SEC, genBest) ;
    //     for (var=0; var<Nvar; var++) printf("%f ", P[gBest][var]);
    printf("%f\n", pF[gBest]);
}

