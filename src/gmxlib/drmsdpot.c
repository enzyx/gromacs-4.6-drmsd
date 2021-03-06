/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by Manuel Luitz, Rainer Bomblies
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

#include "sysstuff.h"
#include "vec.h"
#include "typedefs.h"
#include "main.h"
#include "mtop_util.h"
#include "drmsdpot.h"
#include "smalloc.h"
#include "copyrite.h"
#include "pbc.h"
#include "xvgr.h"
#include "gmxfio.h"

#include "names.h"

void print_drmsd_header(FILE *fp, const output_env_t oenv, real lambda,
        real drmsd_ref, real f_const)
{
    /* open output file and write header lines */
    xvgr_header(fp, "distance RMSD", "Time (ps)", "dRMSD (nm)", exvggtXNY,
            oenv);

    fprintf(fp, "@ subtitle \"lambda = %.3f, reference dRMSD = %.4f nm, "
            "force constant = %6.f kJ/mol/nm^2\"\n", lambda, drmsd_ref,
            f_const);
    fprintf(fp, "@ legend on\n");
    fprintf(fp, "@ legend box on\n");
    fprintf(fp, "@ s0 legend \"distance RMSD\"\n");
    fprintf(fp, "@ s1 legend \"dRMSD potential\"\n");
}

void print_drmsd_data(FILE *fp, real time, real drmsd, real vpot){
    fprintf(fp, "%12.4f %12.7f %12.7f\n", time, drmsd, vpot);
}

/* Initialize I/O */
extern FILE *open_drmsd_out(const char *fn, const t_inputrec *ir,
        const output_env_t oenv, gmx_bool append)
{
    FILE *fp;
    real lambda, L1, drmsd_ref;

    /* Read the lambda value from inputrec */
    if ( ir->efep != efepNO ){
        lambda = ir->fepvals->all_lambda[efptBONDED][ir->fepvals->init_fep_state];
    }
    else
    {
        lambda = 0.0;
    }
    L1 = 1.0 - lambda;

    /* Calculate the lambda dependend reference drmsd */
    drmsd_ref = L1 * ir->drmsd_ref + lambda * ir->drmsd_refB;

    if (append)
    {
        fp = gmx_fio_fopen(fn, "a+");
    }
    else
    {
        fp = gmx_fio_fopen(fn, "w+");
        print_drmsd_header(fp, oenv, lambda, drmsd_ref, ir->drmsd_fc);
    }
    /* Royal flush! */
    fflush(fp);
    return fp;
}

void init_drmsd_pot(FILE *fplog, const gmx_mtop_t *mtop, t_inputrec *ir,
        const t_commrec *cr, gmx_bool bPartDecomp, t_fcdata *fcd,
        gmx_bool bIsREMD)
{
    int nmol, drmsdmols;
    t_ilist *il;
    gmx_mtop_ilistloop_t iloop;
    t_drmsdpotdata *drmsddata;

    /* Get pointer to drmsdpotdata structure */
    drmsddata = &(fcd->drmsdp);

    /* Distance rmsd potential should be only applied to one molecule at a time */
    iloop = gmx_mtop_ilistloop_init(mtop);
    drmsdmols = 0;
    while (gmx_mtop_ilistloop_next(iloop, &il, &nmol))
    {
        if (il[F_DRMSDP].nr > 0)
        {
            drmsdmols += nmol;
        }
    }
    if (drmsdmols > 1)
    {
        gmx_fatal(FARGS,
                "Distance RMSD potential can only be applied to one molecule\n");
    }

    /* Count the total number of distance rmsd interactions in the system */
    drmsddata->npairs = gmx_mtop_ftype_count(mtop, F_DRMSDP);
    if (drmsddata->npairs == 0)
    {
        return;
    }

    if (fplog)
    {
        fprintf(fplog, "Initializing the distance RMSD potential\n");
    }

    /* Feeding constants from inputrec to local structure */
    drmsddata->fc        = ir->drmsd_fc;
    drmsddata->rmsd_ref  = ir->drmsd_ref;
    drmsddata->rmsd_refB = ir->drmsd_refB;

    /* Setup the array of distances */
    snew(drmsddata->dt, drmsddata->npairs);

    /* Bonded interactions between atoms beyond the longest cut-off distance
     * cause problems to the domain decomposition. Tell the user to use
     * particle decomposition instead (Which might be slower) */
    if (cr && PAR(cr) && !bPartDecomp)
    {
        const char *notestr =
                "NOTE: atoms involved in distance RMSD potential restraints should be "
                        "within the longest cut-off distance, if this is not the "
                        "case mdrun generates a fatal error, in that case use "
                        "particle decomposition (mdrun option -pd)";

        if (MASTER(cr))
        {
            fprintf(stderr, "\n%s\n\n", notestr);
        }
        if (fplog)
        {
            fprintf(fplog, "%s\n", notestr);
        }
    }

    /* Tell the user about drmsd potential setup */
    if (fplog)
    {
        fprintf(fplog,
                "There are %d atom pairs involved in the distance rmsd potential\n",
                drmsddata->npairs);
        please_cite(fplog, "Hansdampf2013");
    }
}

void calc_drmsd(int nfa, const t_iatom forceatoms[],
        const t_iparams ip[], const rvec x[], const t_pbc *pbc, t_fcdata *fcd, real lambda)
{
    atom_id ai, aj;
    int fa, type, N;
    rvec dx;
    real d, dmdref, drmsd, L1, dvdlsum;
    real dref, drefB, drmsd_ref, drmsd_refB;
    t_drmsdpotdata *drmsdpotdata;

    drmsdpotdata = &(fcd->drmsdp);

    N       = drmsdpotdata->npairs;
    drmsd   = 0.0;
    dvdlsum = 0.0;
    fa      = 0;

    L1         = 1.0 - lambda;
    drmsd_ref  = drmsdpotdata->rmsd_ref;
    drmsd_refB = drmsdpotdata->rmsd_refB;

    /* This check should not be necessary but to prevent things */
    if (N != (int)(nfa/DIM)){
        gmx_fatal(FARGS,
                "Number of distance RMSD pairs differs from expected!\n"
                "This is most likely a parallelization issue.\n"
                "Try using OpenMP threads instead of MPI or particle decomposition.\n");
    }

    while (fa < nfa)
    {
        type  = forceatoms[fa];        /* instantaneous atom pair */
        ai    = forceatoms[fa + 1];
        aj    = forceatoms[fa + 2];
        dref  = ip[type].drmsdp.dref;  /* The reference distance of atompair ai, aj */
        drefB = ip[type].drmsdp.drefB; /* The reference distance ai, aj in state B */

        /* Get shortest distance ai, aj with respect to periodic boundary conditions */
        if (pbc)
        {
            pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        }
        else
        {
            rvec_sub(x[ai], x[aj], dx);
        }

        d        = sqrt(iprod(dx, dx));
        dmdref   = sqr(d - L1 * dref - lambda * drefB);
        dvdlsum += (d - L1 * dref - lambda * drefB) * (dref - drefB);
        drmsd   += dmdref;

        /* We could save d to drmsdpotdata here but it is not used further */
        //drmsdpotdata->dt[fa/3] = d;
        fa    += DIM;
    }

    /* Calculate the final distance based RMSD */
    drmsd = sqrt(1. / N * drmsd);

    /* Calculate dVdl here because it contains a sum over all drmsd atom pairs */
    drmsdpotdata->dvdl  = drmsdpotdata->fc * (drmsd - L1 * drmsd_ref - lambda * drmsd_refB )
                          * ( dvdlsum / ( N * drmsd ) + drmsd_ref - drmsd_refB );

    /* Calculate drmsd potential here */
    drmsdpotdata->vpot = 0.5 * drmsdpotdata->fc
                         * sqr(drmsd - L1 * drmsd_ref - lambda * drmsd_refB);

    drmsdpotdata->rmsd  = drmsd;
}

real if_drmsd_pot(int npairs, const t_iatom forceatoms[], const t_iparams ip[],
        const rvec x[], rvec f[], rvec fshift[], const t_pbc *pbc,
        const t_graph *g, real lambda, real *dvdlambda, const t_mdatoms *md,
        t_fcdata *fcd, int *global_atom_index)
{
    atom_id ai, aj;
    rvec dx;
    ivec dt;
    int fa, type, ki = CENTRAL, m, N;
    real fc, vtot, vpair, dvdlpair, drmsd_ref, drmsd_refB, drmsd, d, dref, drefB, f_scal, fij, ffc;
    real L1;
    t_drmsdpotdata *drmsdpotdata;

    drmsdpotdata = &(fcd->drmsdp);

    /* Fetch data from drmsdpotdata */
    fc         = drmsdpotdata->fc;
    drmsd_ref  = drmsdpotdata->rmsd_ref;
    drmsd_refB = drmsdpotdata->rmsd_refB;
    drmsd      = drmsdpotdata->rmsd;
    N          = drmsdpotdata->npairs;

    L1 = 1.0 - lambda;

    /* Calculate the potential energy */
    //vpair    = 0.5 * fc * sqr(drmsd - L1 * drmsd_ref - lambda * drmsd_refB);
    vpair    = (drmsdpotdata->vpot) / N;
    dvdlpair = drmsdpotdata->dvdl / N;
    vtot     = 0;

    /* Calculate force prefactor outside the loop */
    ffc = -fc / N * (drmsd - L1 * drmsd_ref - lambda * drmsd_refB) / drmsd;

    /* Loop over drmsd pairs */
    for (fa = 0; fa < npairs; fa += DIM)
    {
        type  = forceatoms[fa];
        ai    = forceatoms[fa + 1];
        aj    = forceatoms[fa + 2];
        dref  = ip[type].drmsdp.dref;
        drefB = ip[type].drmsdp.drefB;

        if (pbc)
        {
            ki = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        }
        else
        {
            rvec_sub(x[ai], x[aj], dx);
        }
        d      = sqrt(iprod(dx, dx));
        f_scal = ffc * (d - L1 * dref - lambda * drefB) / d;

        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            ki = IVEC2IS(dt);
        }

        for (m = 0; m < DIM; m++)
        {
            fij                 = f_scal * dx[m];
            f[ai][m]           += fij;
            f[aj][m]           -= fij;
            fshift[ki][m]      += fij;
            fshift[CENTRAL][m] -= fij;
        }

        /* Add the energy per atom pair to the total energy */
        vtot       += vpair;
        *dvdlambda += dvdlpair;
    }

    return vtot;
}
