/*
 *  diffIRT.c
 *
 *
 *  Created by Dylan Molenaar at the University of Amsterdam.
 *  Copyright 2013 Dylan Molenaar.
 *
 *  This file is part of the diffIRT package for R.
 *
 *  The diffIRT package is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>

void LLdiff(double *y, double *x, int *N, int *nit, double *ai, double *vi, double *sd2A, double *sd2V, double *ter, int *model, double *A, double *W, int *nq, double *epsAlt, double *ou, double *LL);
void DERdiff(double *y, double *x, int *N, int *nit, int *constrain, int *npar, double *ai, double *vi, double *sd2A, double *sd2V, double *ter, int *model, double *A, double *W, int *nq, double *epsAlt, double *delta, double *ou, double *LL, double *deriv);
void LLfacA(double *y, double *x, int *N, int *nit, double *ai, double *vi, double *sd2A, double *sd2V, double *ter, double *Ap, int *model, double *A, double *W, int *nq, double *epsAlt, double *ou, double *LL);
void LLfacV(double *y, double *x, int *N, int *nit, double *ai, double *vi, double *sd2A, double *sd2V, double *ter, double *Vp, int *model, double *A, double *W, int *nq, double *epsAlt, double *ou, double *LL);
void mdiffQQ(double *y, int *N, double *ai, double *vi, double *sd2A, double *sd2V, double *ter, int *model, double *A, double *W,int *nq, double *epsAlt, double *D);
void makeT(int *nit,double *design, int *dim, int *r, double *T);

void ddiff(double *y, double *x, double *ap, double *ai, double *vp, double *vi,double *ter, int *model, double *epsAlt, double *ou){

  double aal=exp(*ap-*ai);
  double yy=(*y-*ter)/(aal*aal);
  double kl, ks, exl=0;

  if(*model==1) exl=(*vp-*vi);      // D-diffusion model
  if(*model==2) exl=exp(*vp-*vi);   // Q-diffusion model

  if((M_PI*yy*(*epsAlt))<1){
    kl=sqrt(-2*log(M_PI*yy*(*epsAlt))/(M_PI*M_PI*yy));
    if(kl < 1/(M_PI*sqrt(yy)) ) kl=1/(M_PI*sqrt(yy));
  } else kl=1/(M_PI*sqrt(yy));
  if((2*sqrt(2*M_PI*yy)*(*epsAlt))<1){
    ks=2+sqrt(-2*yy*log(2*(*epsAlt)*sqrt(2*M_PI*yy)));
    if(ks < (sqrt(yy)+1)) ks=(sqrt(yy)+1);
  } else ks=2;
  double p=0;
  if(ks<kl){
    for(int k=(-floor((ks-1)/2)); k<=ceil((ks-1)/2); k++) p=p+(.5+2*k)*exp(-(.5+2*k)*(.5+2*k)/(2*yy));
    p=p/sqrt(2*M_PI*yy*yy*yy);
  } else{
    int KK=ceil(kl);
    for(int k=1; k<=KK;k++) p=p+(k*exp(-k*k*M_PI*M_PI*yy/2)*sin(k*M_PI*.5));
    p=p*M_PI;
  }
  p=p*1/(aal*aal)*exp(exl*aal*(*x-.5)-.5*exl*exl*(*y-*ter));
  *ou=p;
  return;
}

void ddiffFacV(double *y, double *x, double *ap, double *ai, double *vp, double *vi,double *ter, int *model, double *epsAlt, double *ou){

  double aal=exp(*ap-*ai);
  double yy=(*y-*ter)/(aal*aal);
  double kl, ks, exl=0;

  if(*model==1) exl=(*vp-*vi);      // D-diffusion model
  if(*model==2) exl=(*vp/exp(*vi));   // Q-diffusion model

  if((M_PI*yy*(*epsAlt))<1){
    kl=sqrt(-2*log(M_PI*yy*(*epsAlt))/(M_PI*M_PI*yy));
    if(kl < 1/(M_PI*sqrt(yy)) ) kl=1/(M_PI*sqrt(yy));
  } else kl=1/(M_PI*sqrt(yy));
  if((2*sqrt(2*M_PI*yy)*(*epsAlt))<1){
    ks=2+sqrt(-2*yy*log(2*(*epsAlt)*sqrt(2*M_PI*yy)));
    if(ks < (sqrt(yy)+1)) ks=(sqrt(yy)+1);
  } else ks=2;
  double p=0;
  if(ks<kl){
    for(int k=(-floor((ks-1)/2)); k<=ceil((ks-1)/2); k++) p=p+(.5+2*k)*exp(-(.5+2*k)*(.5+2*k)/(2*yy));
    p=p/sqrt(2*M_PI*yy*yy*yy);
  } else{
    int KK=ceil(kl);
    for(int k=1; k<=KK;k++) p=p+(k*exp(-k*k*M_PI*M_PI*yy/2)*sin(k*M_PI*.5));
    p=p*M_PI;
  }
  p=p*1/(aal*aal)*exp(exl*aal*(*x-.5)-.5*exl*exl*(*y-*ter));
  *ou=p;
  return;
}

void LLdiff(double *y, double *x, int *N, int *nit, double *ai, double *vi, double *sd2A, double *sd2V, double *ter, int *model, double *A, double *W, int *nq, double *epsAlt, double *ou, double *LL){
  double xter, tmp1, tmp2, AtransAp, AtransVp;

  for(int i=0; i<(*N); i++){
    tmp2=0;
      for(int n1=0; n1<(*nq); n1++) for(int n2=0; n2<(*nq); n2++){{
        tmp1=1;
        AtransAp=sqrt(exp(*sd2A))*A[n1];                    // transform the nodes
        AtransVp=sqrt(exp(*sd2V))*A[n2];
        for(int j=0; j<(*nit); j++){
            if(y[i+j*(*N)]!=-999){                //to omit missing values. watchout that a subject cannot have only missings!
              xter=exp(ter[j]);
              ddiff(&y[i+j*(*N)],&x[i+j*(*N)],&AtransAp,&ai[j],&AtransVp,&vi[j],&xter,model,epsAlt,ou);
              tmp1=tmp1*(*ou);                    //calculate likelihood of response vector of each subject
            }
        }
        tmp2=tmp2+tmp1*W[n1]*W[n2];             //calculate marginal likelihood for each subject
      }}
    if(tmp2==0||!R_FINITE(tmp2)) LL[i]=1e-300;
    else LL[i]=tmp2;                            //store marginal likelihood for each subject
  }
  *ou=0;
  for(int i=0; i<(*N); i++) (*ou)+=log(LL[i]);  //take log of the likelihood for each subjects, and add them togheter
  *ou=-2*(*ou);                                 // -2 x log(LL)
  return;
}

void DERdiff(double *y, double *x, int *N, int *nit, int *constrain, int *npar, double *ai, double *vi, double *sd2A, double *sd2V, double *ter, int *model, double *A, double *W, int *nq, double *epsAlt, double *delta, double *ou, double *LL, double *deriv){

  LLdiff(y, x, N, nit, ai, vi, sd2A, sd2V, ter, model, A, W, nq, epsAlt, ou, LL);
  double tmp1=*ou;
  double tmp2=0;


  for(int j=1; j<(*npar+1);j++){      //j starts at 1 !!
    for(int i=0; i<*nit;i++){
      if(constrain[i]==j) ai[i]=ai[i]+*delta;
     if(constrain[i+*nit]==j) vi[i]=vi[i]+*delta;
     if(constrain[i+2*(*nit)]==j) ter[i]=ter[i]+*delta;
    }
    if(constrain[3*(*nit)]==j)   *sd2A=*sd2A+*delta;
    if(constrain[3*(*nit)+1]==j)  *sd2V=*sd2V+*delta;

    LLdiff(y, x, N, nit, ai, vi, sd2A, sd2V, ter, model, A, W, nq, epsAlt, ou, LL);
    tmp2=*ou;
    deriv[j-1]=(tmp2-tmp1)/(*delta);

    for(int i=0; i<*nit;i++){
      if(constrain[i]==j) ai[i]=ai[i]-*delta;
      if(constrain[i+*nit]==j) vi[i]=vi[i]-*delta;
      if(constrain[i+2*(*nit)]==j) ter[i]=ter[i]-*delta;
     }
     if(constrain[3*(*nit)]==j)   *sd2A=*sd2A-*delta;
     if(constrain[3*(*nit)+1]==j)  *sd2V=*sd2V-*delta;
    }
 return;
}


void LLfacA(double *y, double *x, int *N, int *nit, double *ai, double *vi, double *sd2A, double *sd2V, double *ter, double *Ap, int *model, double *A, double *W, int *nq, double *epsAlt, double *ou, double *LL){
  double xter, tmp1, tmp2, AtransVp;

  for(int i=0; i<(*N); i++){
    tmp2=0;
      for(int n1=0; n1<(*nq); n1++){
        tmp1=1;
        AtransVp=sqrt(exp(*sd2V))*A[n1];                    // transform the nodes
        for(int j=0; j<(*nit); j++){
            if(y[i+j*(*N)]!=-999){                //to omit missing values. watchout that a subject cannot have only missings!
              xter=exp(ter[j]);
              ddiff(&y[i+j*(*N)],&x[i+j*(*N)],&Ap[i],&ai[j],&AtransVp,&vi[j],&xter,model,epsAlt,ou);
              tmp1=tmp1*(*ou);                    //calculate likelihood of response vector of each subject
            }
        }
        tmp2=tmp2+tmp1*W[n1];             //calculate marginal likelihood for each subject
      }
    if(tmp2==0||!R_FINITE(tmp2)) LL[i]=1e-300;
    else LL[i]=tmp2*1/sqrt(2*M_PI*(exp(*sd2A)))*exp(-.5*Ap[i]/sqrt(exp(*sd2A)));                            //store marginal likelihood for each subject
  }
  *ou=0;
  for(int i=0; i<(*N); i++) (*ou)+=log(LL[i]);  //take log of the likelihood for each subjects, and add them togheter
  *ou=-2*(*ou);                                 // -2 x log(LL)
  return;
}

void LLfacV(double *y, double *x, int *N, int *nit, double *ai, double *vi, double *sd2A, double *sd2V, double *ter, double *Vp, int *model, double *A, double *W, int *nq, double *epsAlt, double *ou, double *LL){
  double xter, tmp1, tmp2, AtransAp;

  for(int i=0; i<(*N); i++){
    tmp2=0;
      for(int n1=0; n1<(*nq); n1++){
        tmp1=1;
        AtransAp=sqrt(exp(*sd2A))*A[n1];                    // transform the nodes
        for(int j=0; j<(*nit); j++){
            if(y[i+j*(*N)]!=-999){                //to omit missing values. watchout that a subject cannot have only missings!
              xter=exp(ter[j]);
              ddiffFacV(&y[i+j*(*N)],&x[i+j*(*N)],&AtransAp,&ai[j],&Vp[i],&vi[j],&xter,model,epsAlt,ou);
              tmp1=tmp1*(*ou);                    //calculate likelihood of response vector of each subject
            }
        }
        tmp2=tmp2+tmp1*W[n1];             //calculate marginal likelihood for each subject
      }
    if(tmp2==0||!R_FINITE(tmp2)) LL[i]=1e-300;
    else LL[i]=tmp2*1/sqrt(2*M_PI*(exp(*sd2V)))*exp(-.5*Vp[i]/sqrt(exp(*sd2V)));                            //store marginal likelihood for each subject
  }
  *ou=0;
  for(int i=0; i<(*N); i++) (*ou)+=log(LL[i]);  //take log of the likelihood for each subjects, and add them togheter
  *ou=-2*(*ou);                                 // -2 x log(LL)
  return;
}

void ddiffQQ(double *y, double *ap, double *ai, double *vp, double *vi,double *ter, int *model, double *epsAlt, double *ou){
  double aal=exp(*ap-*ai);
  double yy=(*y-*ter)/(aal*aal);
  double kl, ks, exl=0;

  if(*model==1) exl=(*vp-*vi);      // D-diffusion model
  if(*model==2) exl=exp(*vp-*vi);   // Q-diffusion model

    if((M_PI*yy*(*epsAlt))<1){
     kl=sqrt(-2*log(M_PI*yy*(*epsAlt))/(M_PI*M_PI*yy));
      if(kl < 1/(M_PI*sqrt(yy)) ) kl=1/(M_PI*sqrt(yy));
    } else kl=1/(M_PI*sqrt(yy));
    if((2*sqrt(2*M_PI*yy)*(*epsAlt))<1){
      ks=2+sqrt(-2*yy*log(2*(*epsAlt)*sqrt(2*M_PI*yy)));
      if(ks < (sqrt(yy)+1)) ks=(sqrt(yy)+1);
      } else ks=2;
double p=0;
if(ks<kl){
  for(int k=(-floor((ks-1)/2)); k<=ceil((ks-1)/2); k++) p=p+(.5+2*k)*exp(-(.5+2*k)*(.5+2*k)/(2*yy));
  p=p/sqrt(2*M_PI*yy*yy*yy);
} else{
  int KK=ceil(kl);
  for(int k=1; k<=KK;k++) p=p+(k*exp(-k*k*M_PI*M_PI*yy/2)*sin(k*M_PI*.5));
  p=p*M_PI;
}
int x=0;
double p0=p*1/(aal*aal)*exp(exl*aal*(x-.5)-.5*exl*exl*(*y-*ter));
x=1;
double p1=p*1/(aal*aal)*exp(exl*aal*(x-.5)-.5*exl*exl*(*y-*ter));
*ou=p0+p1;
return;
}

void mdiffQQ(double *y, int *N, double *ai, double *vi, double *sd2A, double *sd2V, double *ter, int *model, double *A, double *W,int *nq, double *epsAlt, double *D){

double xter, tmp, A1trans, A2trans;
double ou=0;

for(int i=0; i<(*N); i++){
tmp=0;
  for(int n1=0; n1<(*nq); n1++) for(int n2=0; n2<(*nq); n2++){{

      A1trans=sqrt(exp(*sd2A))*A[n1];     // transform the nodes
      A2trans=sqrt(exp(*sd2V))*A[n2];
        if(y[i]!=-999)         //to omit missing values. watchout that a subject cannot have only missings!
            xter=exp(*ter);
            ddiffQQ(&y[i],&A1trans,ai,&A2trans,vi,&xter,model,epsAlt,&ou); //calculate densities for each observation
        tmp=tmp+ou*W[n1]*W[n2]; //calculate marginal density for each observation (i.e., ap and vp integrated out)
  }}
   if(tmp==0||!R_FINITE(tmp)) D[i]=-999;
   else D[i]=tmp;                  //store marginal density for each response
}
  return;
}

void equal(double *big, double *small, int *dim,int *ou){
//length small (lS) should be smaller or equal to LB;
int d=0;
 int lS=small[0];
 int lB=big[0];
 for(int s=1; s<=(lS); s++){
 for(int b=1; b<=(lB); b++){
    if(big[*dim*b]==small[*dim*s]) d=d+1;
 }}

 if(d==(lS)) *ou=1;
 else *ou=0;
 }

void makeT(int *nit,double *design, int *dim, int *r, double *T){

int ou;

for(int i=0; i<*r;i++){
for(int j=0; j<*dim;j++){

if(i<=j){
   equal(&design[j],&design[i],dim,&ou);
 T[i+j*(*r)]=ou;

 }
 }
 }
}



