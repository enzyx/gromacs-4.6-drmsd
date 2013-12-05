/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by Rainer Reisenauer
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
#include <math.h>
#include <string.h>

#include "sysstuff.h"
#include "gstat.h"
#include "macros.h"
#include "main.h"
#include "mshift.h"
#include "mtop_util.h"
#include "pbc.h"
#include "smalloc.h"
#include "statutil.h"
#include "tpxio.h"
#include "typedefs.h"
#include "vec.h"
#include "xvgr.h"
#include "gmx_ana.h"
#include "copyrite.h"

#include "drmsdpot.h"

typedef struct t_drmsd_data
{
    real  t, drmsd, vpot;
    struct t_drmsd_data *next;
} t_drmsd_data;

void init_drmsd_data_element(t_drmsd_data *ddat){
	ddat->next = NULL;
	ddat->t = 1;
	ddat->drmsd = 1;
	ddat->vpot = 1;
}

void dump_drmsd_xvg(const char *filename, const output_env_t oenv, t_drmsd_data *ddat_head, real drmsd_ref)
{
#define NFILE asize(fnm)
	FILE *out = NULL;
	out = xvgropen(filename, "Read variable", "Time (ps)", "dRMSD (nm)", oenv); /* open output file */

	t_drmsd_data *ddat = ddat_head;
	while (ddat->next !=NULL){
		fprintf(out, "%.0f \t %7f \t %7f \t%7f\n", ddat->t, ddat->drmsd, drmsd_ref, ddat->vpot);
		ddat=ddat->next;
	}

	ffclose(out);
}

int gmx_drmsd(int argc, char *argv[])
{
    const char     *desc[] = {
        "[TT]g_drmsd[tt] computes distances and potentials for the distanceRMSD.",
        "If given a .gro coordinate and index file it calculates the reference ",
        "distances and writes them to a topology"
    };
    t_pargs         pa[]      = {

    };

    FILE           *out = NULL;
    t_tpxheader     header;
    t_inputrec      ir;		/*input record, mdp options*/
    gmx_mtop_t      mtop;	/* topology */
    rvec           *xtop;	/*vector with entry for every atom */
    gmx_localtop_t *top;	/* The node-local topology struct */
    t_fcdata        fcd;
    t_commrec      *cr;
    t_graph        *g;
    int             ntopatoms, natoms;
    t_trxstatus    *status;
    t_drmsd_data   *ddat_head=NULL, *ddat, *exdd, *exhead;
    real            t;
    rvec           *x;
    matrix          box;
    t_pbc           pbc, *pbc_null;
    output_env_t    oenv;
    gmx_rmpbc_t     gpbc = NULL;

    /* for potential calculation */
    real vpot;
    rvec f[2*DIM*fcd.drmsdp.npairs], fshift[DIM*CENTRAL];
    real lambda;
    real dvdlambda = 0;
    t_mdatoms *md;
    int *global_atom_index = 0;

    t_filenm        fnm[] = {
    	{ efTPX, "-s", NULL, ffREAD },
    	{ efTRX, "-f", NULL, ffREAD },
    	{ efXVG, "-o",  "drmsd", ffWRITE }
    };

#define NFILE asize(fnm)

    cr  = init_par(&argc, &argv); /* communication between multiple nodes */
    CopyRight(stderr, argv[0]);

    /* parsing arguments PCA are common options like -b, -e, nicelevel
     * NFILE is the number of files, fnm is an array of filenames
     * oenv is the output environment*/
    parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

    read_tpxheader(ftp2fn(efTPX, NFILE, fnm), &header, FALSE, NULL, NULL);
    snew(xtop, header.natoms); /* allocate memory for xtop vector */
    read_tpx(ftp2fn(efTPX, NFILE, fnm), &ir, box, &ntopatoms, xtop, NULL, NULL, &mtop); // read topology from tpr file

    lambda = ir.fepvals->all_lambda[ir.fepvals->init_fep_state][efptBONDED];

    //fprintf(stderr,"\n\n bonded_lambda = %f\n\n",lambda);


    //fprintf(stderr,"drmsd reference is %f\n",ir.drmsd_ref);
    //fprintf(stderr,"drmsd reference for state B is %f\n",ir.drmsd_refB);
    //fprintf(stderr,"drmsd fc is %f\n",ir.drmsd_fc);


    top = gmx_mtop_generate_local_top(&mtop, &ir); /* generate local topology */

    /* for periodic boundary conditions */
    g        = NULL;
    pbc_null = NULL;
    if (ir.ePBC != epbcNONE)
    {
        if (ir.bPeriodicMols)
        {
            pbc_null = &pbc;
        }
        else
        {
            g = mk_graph(stderr, &top->idef, 0, mtop.natoms, FALSE, FALSE);
        }
    }

    init_drmsd_pot(stderr, &mtop, &ir, cr, 0, &fcd, 0);

    /* initialise first ddat structure */
	snew(ddat_head,1);
	init_drmsd_data_element(ddat_head);
	ddat = ddat_head;

    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box); /* read number of atoms in system */

    // Values from drmsd calculation
    //fprintf(stderr, "fcd.drmsdp.rmsd     = %7f\n", fcd.drmsdp.rmsd);
    //fprintf(stderr, "fcd.drmsdp.rmsd_ref = %7f\n", fcd.drmsdp.rmsd_ref);
    //fprintf(stderr, "fcd.drmsdp.npairs   = %7d\n", fcd.drmsdp.npairs);
    //fprintf(stderr, "fcd.drmsdp.fc       = %7f\n", fcd.drmsdp.fc);
    //fprintf(stderr, "fcd.drmsdp.dt       = %7f\n", fcd.drmsdp.dt[0]);


    if (ir.ePBC != epbcNONE)
    {
        gpbc = gmx_rmpbc_init(&top->idef, ir.ePBC, natoms, box);
    }

    do
    {
        //fprintf(stderr,"starting while %d\n",i);
    	/* fprintf(stderr,"got time %f\n",t); */
    	/* fprintf(stderr,"trxstatus->nxframe = %f\n",status->NATOMS); */
    	/* fprintf(out,  "%.0f  %7d\n", t, status->xframe->x); */


        if (ir.ePBC != epbcNONE)
        {
            if (ir.bPeriodicMols)
            {
                set_pbc(&pbc, ir.ePBC, box);
            }
            else
            {
                gmx_rmpbc(gpbc, natoms, box, x);
            }
        }

        //fprintf(stderr,"calculating drmsd %d\n",i);
        calc_drmsd(top->idef.il[F_DRMSDP].nr, top->idef.il[F_DRMSDP].iatoms,
        		top->idef.iparams, (const rvec*) x, pbc_null, &fcd, lambda);


        /* //output of memory locations
        fprintf(stderr,"0x%8x next \t 0x%8x t \t%f \t 0x%8x drmsd \t%f \t 0x%8x vpot \t%f\n",\
        		ddat->next,ddat->t,ddat->t,ddat->drmsd,ddat->drmsd,ddat->vpot,ddat->vpot);
        		*/

        snew(ddat->next,1);
        ddat = ddat->next;
        init_drmsd_data_element(ddat);

        ddat->drmsd = fcd.drmsdp.rmsd;
        ddat->t = t;
        ddat->vpot = if_drmsd_pot(fcd.drmsdp.npairs, top->idef.il[F_DRMSDP].iatoms,
        		top->idef.iparams, (const rvec*) x, f, fshift, pbc_null, g,
        		lambda, &dvdlambda, NULL, &fcd, global_atom_index);

        //output of entire list so far every iteration
        /*t_drmsd_data *read_ddat = ddat_head->next;
        read_ddat = ddat_head->next;
    	while (read_ddat->next !=NULL){
    		//fprintf(stderr, "%.0f \t %7f \t %7f \t%7f -- ddat_head\n", ddat_head->t, ddat_head->drmsd, fcd.drmsdp.rmsd_refB, ddat_head->vpot);
    		fprintf(stderr, "%.0f \t %7f \t %7f \t%7f\n", read_ddat->t, read_ddat->drmsd, fcd.drmsdp.rmsd_refB, read_ddat->vpot);
    		read_ddat=read_ddat->next;
    	}*/

        //fprintf(stderr,"%f \t %f \t %f\n",t,ddat->drmsd,ddat->vpot);
    }
    while (read_next_x(oenv, status, &t, natoms, x, box));

    //dump output to dmrsd.xvg
    dump_drmsd_xvg(opt2fn("-o", NFILE, fnm), oenv, ddat_head->next, fcd.drmsdp.rmsd_ref);

    /* output to screen
    ddat = ddat_head;
    do
    {
        fprintf(stderr,"ddat->t = %f; \t ddat->drmsd = %f; \t ddat->vpot= %f \n",
        		ddat->t,ddat->drmsd,ddat->vpot);
        ddat = ddat->next;
    } while(ddat->next!=NULL);
    */



    close_trj(status);
    if (ir.ePBC != epbcNONE)
    {
        gmx_rmpbc_done(gpbc);
    }

    ffclose(out);

    gmx_finalize_par();


    fprintf(stderr,"\nit.works\n");

    thanx(stderr);

	return 0;
}
