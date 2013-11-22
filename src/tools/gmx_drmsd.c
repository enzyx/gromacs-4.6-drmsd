#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <string.h>

#include "sysstuff.h"
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

#include "drmsdpot.h"

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
    real            t;
    rvec           *x;
    matrix          box;
    t_pbc           pbc, *pbc_null;
    FILE           *fplog;
    output_env_t    oenv;



    t_filenm        fnm[] = {
    	{ efTPX, NULL, NULL, ffREAD },
    	{ efTRX, "-f", NULL, ffREAD },
    	{ efLOG, "-l",  "drmsd", ffWRITE },
    	{ efXVG, "-o",  "drmsd", ffWRITE }
    };

#define NFILE asize(fnm)

    cr  = init_par(&argc, &argv); /* communication between multiple nodes */
    /* CopyRight(stderr, argv[0]); */

    /* parsing arguments PCA are common options like -b, -e, nicelevel
     * NFILE is the number of files, fnm is an array of filenames
     * oenv is the output environment*/
    parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

    gmx_log_open(ftp2fn(efLOG, NFILE, fnm), cr, FALSE, 0, &fplog); /* initialize log */

    read_tpxheader(ftp2fn(efTPX, NFILE, fnm), &header, FALSE, NULL, NULL);
    snew(xtop, header.natoms); /* allocate memory for xtop vector */
    read_tpx(ftp2fn(efTPX, NFILE, fnm), &ir, box, &ntopatoms, xtop, NULL, NULL, &mtop); /* read topology from tpr file */

    fprintf(stderr,"drmsd reference is %f\n",ir.drmsd_ref);
    fprintf(stderr,"drmsd fc is %f\n",ir.drmsd_fc);


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
            g = mk_graph(fplog, &top->idef, 0, mtop.natoms, FALSE, FALSE);
        }
    }

    out = xvgropen(opt2fn("-o", NFILE, fnm), "Read variable", "Time (ps)", "[?]", oenv); /* open output file */

    init_drmsd_pot(fplog, &mtop, &ir, cr, 0, &fcd, 0);

    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box); /* read number of atoms in system */

    history_t      *hist;
    calc_drmsd_pot( cr->ms, top->idef.il[F_DRMSDP].nr, top->idef.il[F_DRMSDP].iatoms, top->idef.iparams,
    		(const rvec*) x, pbc_null, &fcd, hist);

    fprintf(stderr, "\n fcd.drmsdp.rmsd     = %7f\n", fcd.drmsdp.rmsd);
    fprintf(stderr, "\n fcd.drmsdp.rmsd_ref = %7f\n", fcd.drmsdp.rmsd_ref);
    fprintf(stderr, "\n fcd.drmsdp.npairs   = %7f\n", fcd.drmsdp.npairs);
    fprintf(stderr, "\n fcd.drmsdp.fc       = %7f\n", fcd.drmsdp.fc);
    fprintf(stderr, "\n fcd.drmsdp.dt       = %7f\n", fcd.drmsdp.dt);


    /* read out atoms of binding
    int j;
    int isize = 5;
    for (j = 1; (j < isize); j+=3)
        {
    		fprintf(stderr,"\n value = %d",top->idef.il[F_DRMSDP].iatoms[j]);
        }
	*/

    while (read_next_x(oenv, status, &t, natoms, x, box)){
    	/* fprintf(stderr,"got time %f\n",t); */
    	/* fprintf(stderr,"trxstatus->nxframe = %f\n",status->NATOMS); */
    	/* fprintf(out,  "%.0f  %7d\n", t, status->xframe->x); */
    }
    close_trj(status);

    /* fprintf(stderr,  "%.0f  %7d\n", t, status->BOX[0]); */

    ffclose(out);

    gmx_finalize_par();

    gmx_log_close(fplog);

    fprintf(stderr,"\nit.works\n");

	return 0;
}
