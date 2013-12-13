/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by Rainer Reisenauer, Manuel Luitz
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
#include "names.h"

#include "drmsdpot.h"
#include "gmx_wham.h"

typedef struct t_drmsd_data
{
    real  t, drmsd, vpot;
    struct t_drmsd_data *next;
} t_drmsd_data;

typedef struct t_drmsd_head
{
	real drmsd_ref;
	real f_const;
	real lambda;
	struct t_drmsd_data *next;

}t_drmsd_head;

void init_drmsd_data_elemnt(t_drmsd_data *ddat)
{
	ddat->next  = NULL;
	ddat->t     = 0;
	ddat->drmsd = 0;
	ddat->vpot  = 0;
}

void init_drmsd_head(t_drmsd_head *dhead)
{
	dhead->next = NULL;
	dhead->drmsd_ref = 0;
	dhead->f_const = 0;
	dhead->lambda = 0;
}

void read_fid(char* fid, char* fn){
	/*find file identifier after the last _
	 * and clear file type (.xtc) of filename */

	if(strchr(fn, '_')){
		char *fid_tmp;
		snew(fid_tmp,STRLEN);
		strcpy(fid_tmp,strrchr(fn, '_')+1);
		sprintf(fid   ,"_%.*s",(int)(strlen(fid_tmp)-4),fid_tmp);
	}
	else {
		fid = "";
	}
}

void read_fileset(int set_number, char *ftprnm,
				  char *ftrxnm, output_env_t oenv,
				  t_drmsd_head *ddat_head)
{
    t_inputrec      ir;		            /* input record */
    gmx_mtop_t      mtop;	            /* topology */
    gmx_localtop_t *top;	            /* The node-local topology struct */
    t_fcdata        fcd;
    real            t;
    rvec           *x;
    matrix          box;
    t_pbc          *pbc;
    gmx_rmpbc_t     gpbc = NULL;
    t_trxstatus    *status;
    int             ntopatoms, natoms;
    int             ePBC;
    t_drmsd_data   *ddat;

    /* for potential calculation */
    real vpot, drmsd_ref, lambda, L1, dvdlambda = 0.0;
    rvec  fshift[DIM*CENTRAL];
    rvec *f = NULL;
    int *global_atom_index = 0;


    /* read topology from tpr file */
    //ePBC = read_tpx(ftp2fn(efTPX, NFILE, fnm), &ir, box, &ntopatoms, NULL, NULL, NULL, &mtop);
    ePBC = read_tpx(ftprnm, &ir, box, &ntopatoms, NULL, NULL, NULL, &mtop);


    /* Read the lambda value from inputrec */
    if ( ir.efep != efepNO ){
        lambda = ir.fepvals->all_lambda[efptBONDED][ir.fepvals->init_fep_state];
    }
    else
    {
        lambda = 0.0;
    }
    L1 = 1.0 - lambda;

    /* generate local topology */
    top = gmx_mtop_generate_local_top(&mtop, &ir);

    /* Initialize the drmsd potential data structure in fcd
     * and setting bIsREMD to false
     */
    init_drmsd_pot(stderr, &mtop, &ir, NULL, 0, &fcd, FALSE);

    /* Calculate the lambda dependent reference drmsd */
    drmsd_ref = L1 * fcd.drmsdp.rmsd_ref + lambda * fcd.drmsdp.rmsd_refB;

    /* We need a dummy force vector for the ifunc call */
    snew(f, mtop.natoms);

    /* save drmsd run parameters */
    ddat_head->f_const = fcd.drmsdp.fc;
    ddat_head->lambda = lambda;
    ddat_head->drmsd_ref = drmsd_ref;

    /* Initialize first ddat structure */
    snew(ddat, 1);
    init_drmsd_data_elemnt(ddat);
    ddat_head->next = ddat;

    /* read number of atoms in system */
    natoms = read_first_x(oenv, &status, ftrxnm, &t, &x, box);

    /* for periodic boundary conditions */
    if (ePBC != epbcNONE)
    {
        snew(pbc, 1);
    }
    else
    {
        pbc = NULL;
    }
    gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms, box);

    do
    {
        /* Initialization for correct distance calculations */
        if (pbc)
        {
            set_pbc(pbc, ePBC, box);
            /* make molecules whole again */
            gmx_rmpbc(gpbc, natoms, box, x);
        }

        calc_drmsd(top->idef.il[F_DRMSDP].nr, top->idef.il[F_DRMSDP].iatoms,
        		top->idef.iparams, (const rvec*) x, pbc, &fcd, lambda);

        ddat->drmsd = fcd.drmsdp.rmsd;
        ddat->t     = t;
        ddat->vpot  = if_drmsd_pot(top->idef.il[F_DRMSDP].nr, top->idef.il[F_DRMSDP].iatoms,
        		top->idef.iparams, (const rvec*) x, f, fshift, pbc, NULL,
        		lambda, &dvdlambda, NULL, &fcd, global_atom_index);

        /* Extend the linked list */
        snew(ddat->next, 1);
        ddat = ddat->next;
        init_drmsd_data_elemnt(ddat);
    }
    while (read_next_x(oenv, status, &t, natoms, x, box));

    if (ePBC != epbcNONE)
    {
        gmx_rmpbc_done(gpbc);
    }

    close_trj(status);
}

void dump_drmsd_xvg(const char *filename, char *f_id, const output_env_t oenv,
        t_drmsd_head *d_head)
{
#define NFILE asize(fnm)
    FILE *out = NULL;
    char buf[STRLEN];
    sprintf(buf,"%.*s%s.xvg",(int)(strlen(filename)-4),filename,f_id);

    /* open output file and write header lines */
    out = xvgropen(buf, "distance RMSD", "Time (ps)", "dRMSD (nm)", oenv);
    fprintf(out,"@ subtitle \"lambda = %.3f, reference dRMSD = %.4f nm, force constant = %6.f kJ/mol/nm^2\"\n",
    		d_head->lambda, d_head->drmsd_ref, d_head->f_const);
    fprintf(out,"@ legend on\n");
    fprintf(out,"@ legend box on\n");
    fprintf(out,"@ s0 legend \"distance RMSD\"\n");
    fprintf(out,"@ s1 legend \"dRMSD potential\"\n");

    t_drmsd_data *ddat = d_head->next;
    while (ddat->next != NULL)
    {
        fprintf(out, "%12.7f %12.7f %12.7f\n", ddat->t, ddat->drmsd,
                ddat->vpot);
        ddat = ddat->next;
    }

    ffclose(out);
}

void do_wham(const output_env_t oenv, int nwin, t_drmsd_head **dd_head){
    t_UmbrellaHeader         header;
    t_UmbrellaWindow        *window;
}

void free_ddata(t_drmsd_head **ddat_head, int nfilesets){
	int i;
	t_drmsd_data *ddat;
	t_drmsd_data *ddat_del;
	for(i = 0; i < nfilesets; i++){
		ddat = ddat_head[i]->next;
		while(ddat->next != NULL){
			ddat_del = ddat;
			ddat = ddat-> next;
			sfree(ddat_del);
		}
		ddat_del = ddat;
		sfree(ddat_del);
	}
	sfree(ddat_head);
}

int gmx_drmsd(int argc, char *argv[])
{
    const char     *desc[] = {
        "[TT]g_drmsd[tt] computes distances and potentials for the distanceRMSD.",
        "If given a .gro coordinate and index file it calculates the reference ",
        "distances and writes them to a topology [not implemented yet]",
        "-o defines a filename, default is drmsd.xvg",
        "	if multiple files are given, the filename is extended by the",
        "	file identifier of the trajectory after the last underscore and",
        "	excluding the file extension."
    };
    t_pargs         pa[]      = {
    };
    FILE           *out = NULL;
    int 			ntrxfile = 0, ntprfile = 0;
    int 			i;
    char       	  **ftrxnms, **ftprnms;
    char 		  **file_ids;
    gmx_bool		bf_given_index;
    t_drmsd_head  **ddat_head;
    output_env_t    oenv;

    t_filenm        fnm[] = {
    	{ efTPX, "-s", NULL, ffOPTRDMULT },
    	{ efTRX, "-f", NULL,  ffOPTRDMULT},
    	{ efXVG, "-o",  "drmsd", ffWRITE }
    };

#define NFILE asize(fnm)
    CopyRight(stderr, argv[0]);

    /* parsing arguments PCA are common options like -b, -e, nicelevel
     * NFILE is the number of files, fnm is an array of filenames
     * oenv is the output environment
     */
    parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

    /* get names of trx files */
    if (opt2bSet("-f", NFILE, fnm))
    {
        ntrxfile = opt2fns(&ftrxnms, "-f", NFILE, fnm);
    }
    /* get names of tpr files */
    if (opt2bSet("-s", NFILE, fnm))
    {
        ntprfile = opt2fns(&ftprnms, "-s", NFILE, fnm);
    }

    if(ntprfile != ntrxfile){
    	fprintf(stderr,"Error: number of tpr files does not match number of trajectories.\n");
    	return 0;
    }

    /* initialize array of data heads */
    snew(ddat_head,ntrxfile);
    /*initiallize array of file identifiers */
    snew(file_ids,ntrxfile);

    /* read trajectories */
    for(i = 0; i < ntrxfile; i++){

		snew(ddat_head[i],1);
		init_drmsd_head(ddat_head[i]);

		snew(file_ids[i],STRLEN);
		read_fid(file_ids[i],ftrxnms[i]);
		fprintf(stderr,"returned = %s\n",file_ids[i]);

		fprintf(stderr,"file_ids[i] got = %s\n",file_ids[i]);
		read_fileset(i, ftprnms[i], ftrxnms[i], oenv, ddat_head[i]);
    }

    /* dump output to dmrsd.xvg_i */
    for(i = 0; i < ntrxfile; i++){
    	dump_drmsd_xvg(opt2fn("-o", NFILE, fnm), file_ids[i], oenv, ddat_head[i]);
    }

    //do_wham(oenv, ntrxfile, ddat_head);

    gmx_finalize_par();

    free_ddata(ddat_head, ntrxfile);

    thanx(stderr);

	return 0;
}
