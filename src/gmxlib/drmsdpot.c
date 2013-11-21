/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by Manuel Luitz
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

#include "names.h"

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
    drmsddata->fc = ir->drmsd_fc;
    drmsddata->rmsd_ref = ir->drmsd_ref;

    /* Setup the array of distances */
    snew(drmsddata->dt, drmsddata->npairs);

    /* Bonded interactions between atoms beyond the longest cut-off distance
     * cause problems to the domain decomposition. Tell the user to use
     * particle decomposition instead (Which might be slower) */
    if (cr && PAR(cr) && !bPartDecomp)
    {
        const char *notestr =
                "NOTE: atoms involved in distance rmsd potential restraints should be "
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

        if (ir->nstdrmsdpout != 0)
        {
            if (fplog)
            {
                fprintf(fplog,
                        "\nWARNING: Can not write distance rmsd restraint data to "
                                "energy file with domain decomposition\n\n");
            }
            if (MASTER(cr))
            {
                fprintf(stderr,
                        "\nWARNING: Can not write distance rmsd restraint data to "
                                "energy file with domain decomposition\n");
            }
            ir->nstdrmsdpout = 0;
        }
    }

    /* Tell the user about drmsd potential setup */
    fprintf(stderr,
            "There are %d atom pairs involved in the distance rmsd potential\n",
            drmsddata->npairs);
    please_cite(fplog, "Hansdampf2013");
}

void calc_drmsd_pot(const gmx_multisim_t *ms, int nfa, const t_iatom forceatoms[],
        const t_iparams ip[], const rvec x[], const t_pbc *pbc, t_fcdata *fcd,
        history_t *hist)
{
    atom_id ai, aj;
    int fa, type;
    rvec dx;
    real d, dmdref, dRMSD;
    real dref;
    t_drmsdpotdata *drmsdpotdata;

    drmsdpotdata = &(fcd->drmsdp);

    dRMSD = 0.0;
    fa = 0;

    while (fa < nfa)
    {
        type = forceatoms[fa]; /* instantaneous atom pair */
        ai = forceatoms[fa + 1];
        aj = forceatoms[fa + 2];
        dref = ip[type].drmsdp.dref; /* The reference distance of atompair ai, aj */

        /* Get shortest distance ai,aj with respect to periodic boundary conditions */
        if (pbc)
        {
            pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        }
        else
        {
            rvec_sub(x[ai], x[aj], dx);
        }

        d = sqrt(iprod(dx, dx));
        dmdref = sqr(d - dref);
        dRMSD += dmdref;

        /* We can save d to drmsdpotdata here but it is not used further */
        //drmsdpotdata->dt[fa/3] = d;
        fa += 3;
    }

    dRMSD = sqrt(1. / (drmsdpotdata->npairs) * dRMSD);

    drmsdpotdata->rmsd = dRMSD;
}

real ta_drmsd_pot(int npairs, const t_iatom forceatoms[], const t_iparams ip[],
        const rvec x[], rvec f[], rvec fshift[], const t_pbc *pbc,
        const t_graph *g, real lambda, real *dvdlambda, const t_mdatoms *md,
        t_fcdata *fcd, int *global_atom_index)
{
    atom_id ai, aj;
    rvec dx;
    ivec dt;
    int fa, type, ki = CENTRAL, m;
    real fc, vtot, vpair, rmsd_ref, rmsd, d, dref, f_scal, fij;
    t_drmsdpotdata *drmsdpotdata;

    drmsdpotdata = &(fcd->drmsdp);

    //TODO write to file

    fc = drmsdpotdata->fc;
    rmsd_ref = drmsdpotdata->rmsd_ref;
    rmsd = drmsdpotdata->rmsd;

    /* Calculate the potential energy */
    vpair = 0.5 * fc * sqr(rmsd - rmsd_ref);
    vtot = 0;

    /* Loop over drmsd pairs */
    for (fa = 0; fa < npairs; fa += 3)
    {
        type = forceatoms[fa];
        ai = forceatoms[fa + 1];
        aj = forceatoms[fa + 2];
        dref = ip[type].drmsdp.dref;

        if (pbc)
        {
            ki = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        }
        else
        {
            rvec_sub(x[ai], x[aj], dx);
        }
        d = sqrt(iprod(dx, dx));
        f_scal = -fc / (npairs / 3.) * (rmsd - rmsd_ref) / (rmsd) * (d - dref)
                / d;

        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), dt);
            ki = IVEC2IS(dt);
        }

        for (m = 0; m < DIM; m++)
        {
            fij = f_scal * dx[m];

            f[ai][m] += fij;
            f[aj][m] -= fij;
            fshift[ki][m] += fij;
            fshift[CENTRAL][m] -= fij;
        }

        /* Add the energy per atom pair to the total energy */
        vtot += vpair;
    }

    return vtot;
}
