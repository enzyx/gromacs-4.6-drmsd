/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "main.h"
#include "vec.h"
#include "txtdump.h"

void calc_force(int natom,rvec f[],rvec fff[])
{
  int  i,j,m;
  int  jindex[] = { 0, 5, 10};
  rvec dx,df;
  real msf1,msf2;
  
  for(j=0; (j<2); j++) {
    clear_rvec(fff[j]);
    for(i=jindex[j]; (i<jindex[j+1]); i++) {
      for(m=0; (m<DIM); m++) {
	fff[j][m] += f[i][m];
      }
    }
  }
  
  msf1 = iprod(fff[0],fff[0]);
  msf2 = iprod(fff[1],fff[1]);
  if (debug) {
    pr_rvecs(debug,0,"force",f,natom);
    
    fprintf(debug,"FMOL:  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f\n",
	    fff[0][XX],fff[0][YY],fff[0][ZZ],fff[1][XX],fff[1][YY],fff[1][ZZ]);
    fprintf(debug,"RMSF:  %10.3e  %10.3e\n",msf1,msf2);
  }
}


void calc_f_dev(int natoms,real charge[],rvec x[],rvec f[],
		t_idef *idef,real *xiH,real *xiS)
{
  enum { wwwO, wwwH1, wwwH2, wwwS, wwwNR };
  int  lj_index[wwwNR] = { 0,    4,     4,     8 };
  real rmsf,dFH,dFS,dFO,dr_14;
  real q[wwwNR],c6[wwwNR],c12[wwwNR],c12ratio;
  rvec fff[2],dx;
  int  i,j,aj;
  
  for(i=0; (i<wwwNR); i++) {
    q[i]   = charge[i];
    c12[i] = idef->iparams[lj_index[i]].lj.c12;
    c6[i]  = idef->iparams[lj_index[i]].lj.c6;
  }
  
  calc_force(natoms,f,fff);
  rmsf = norm(fff[0]);
  dFS = 0;
  dFH = 0;
  dFO = 0;
  for(i=0; (i<4); i++) {
    for(j=4; (j<8); j++) {
      if (c12[i] != 0) {
	rvec_sub(x[i],x[j],dx);
	aj       = j % 4;
	dr_14    = pow(iprod(dx,dx),-7);
	c12ratio = -12*sqrt(c12[aj]/c12[i])*dr_14*iprod(fff[0],dx);
	
	switch (i) {
	case wwwH1:
	case wwwH2:
	  dFH += c12ratio;
	  break;
	case wwwS:
	  dFS += c12ratio;
	  break;
	case wwwO:
	  dFS += c12ratio;
	  break;
	}
      }
    }
  }

  if (debug) {    
    fprintf(debug,"FFF: dFS=%10.3e,  dFH=%10.3e,  dFO=%10.3e, rmsf=%10.3e\n",
	    dFS,dFH,dFO,rmsf);
  }
  if (dFH == 0)
    *xiH = 1;
  else
    *xiH=rmsf/(10*c12[wwwH1]*dFH);
    
  if (dFS == 0)
    *xiS = 1;
  else
    *xiS=rmsf/(10*c12[wwwS]*dFS);
}
