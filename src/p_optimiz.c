/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "p_optimiz.h"

/*********************************************************/

void Opt_Dist_F(phydbl *dist, phydbl *F, model *mod)
{
  phydbl ax,bx,cx;
  if(*dist < BL_MIN) *dist = BL_MIN;

  ax = *dist;
  bx = 1.5*(*dist);

  ax = 10.*(*dist);
  bx =     (*dist);
  cx = .10*(*dist);

  Dist_F_Brent(ax,bx,cx,1.E-10,1000,dist,F,mod);
  
  return;
}

/*********************************************************/

phydbl Dist_F_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
	phydbl *param, phydbl *F, model *mod)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL, init_lnL, curr_lnL;

//optimize distance, not likelihood
  phydbl old_param, cur_param;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x = w = v = bx;
  old_lnL = UNLIKELY;
  fw = fv = fx = -Lk_Dist(F,fabs(bx),mod);
  curr_lnL = init_lnL = -fw;

//optimize distance
  old_param = cur_param = fabs(bx);

  for(iter=1;iter<=BRENT_ITMAX;iter++)
    {
      xm=0.5*(a+b);

      tol2=2.0*(tol1=tol*fabs(x)+BRENT_ZEPS);

//optimize distance
      //if ((iter > 1) && fabs (old_param-cur_param) < mod->s_opt->min_diff_lk_local)
      if ((iter > 1) && fabs (old_param-cur_param) < 1.E-06)
	{
	  *param = x;
	  curr_lnL = Lk_Dist(F,*param,mod);
	  return -curr_lnL;
	}
      
      if(fabs(e) > tol1) 
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=fabs(q);
	  etemp=e;
	  e=d;
	  if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    {
	      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	      //printf("Golden section step\n");
	    }
	  else
	    {
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-x);
	      //printf("Parabolic step\n");
	    }
        }
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  //printf("Golden section step (default)\n");
	}
      
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      if(u<BL_MIN) u = BL_MIN;
      (*param) = fabs(u);
      old_lnL = curr_lnL;
      fu = -Lk_Dist(F,fabs(u),mod);
      curr_lnL = -fu;
      
      if(fu <= fx) 
	{
	  if(iter > n_iter_max) return -fu;

	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	} 
      else
	{
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) 
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
	  else if (fu <= fv || v == x || v == w) {
            v=u;
            fv=fu;
	  }
	}
//optimize distance
      old_param = cur_param;
      cur_param = *param;
    }
  Exit("Too many iterations in BRENT.");
  return(-1);
}

/*********************************************************/
