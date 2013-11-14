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

#include "names.h"


void init_drmsd_pot(FILE *fplog, const gmx_mtop_t *mtop,
                 t_inputrec *ir, const t_commrec *cr, gmx_bool bPartDecomp,
                 t_fcdata *fcd, gmx_bool bIsREMD)
{
    t_drmsdpotdata  *dd;

    /* Get pointer to drmsdpotdata structure */
    dd = &(fcd->drmsdp);

    /* START DEBUGING */
    fprintf(stderr, "Initializing the distance rmsd parameters\n");
    fprintf(stderr, "drmsd-pot: %s\n", EBOOL(ir->bDrmsdPot));
    fprintf(stderr, "drmsd-ref: %f\n", ir->drmsd_ref);
    fprintf(stderr, "drmsd-fc: %f\n", ir->drmsd_fc);
    fprintf(stderr, "nstdrmsdpout: %d\n", ir->nstdrmsdpout);
    fprintf(stderr, "We have %d drmsd pairs\n", gmx_mtop_ftype_count(mtop, F_DRMSDP));

    t_iparams *ip;
    t_ilist *il;
    int nmol, i;
    gmx_mtop_ilistloop_t iloop;
    iloop     = gmx_mtop_ilistloop_init(mtop);
    while(gmx_mtop_ilistloop_next(iloop, &il, &nmol))
        {
            fprintf(stderr, "Number of molecules of this type: %d\n", nmol);
            fprintf(stderr, "Number of drmsd interactions in this moleculetype: %d\n", il[F_DRMSDP].nr);
            for (i = 0; i < il[F_DRMSDP].nr; i += 3)
            {
                fprintf(stderr, "Reference distance drmsd: %f\n", mtop->ffparams.iparams[il[F_DRMSDP].iatoms[i]].drmsdp.dref);
            }
        }
    /* END */

    /* Count the total number of distance rmsd interactions in the system */
    if (gmx_mtop_ftype_count(mtop, F_DRMSDP) == 0)
        {
            dd->nres = 0;
            return;
        }

    if (fplog)
    {
        fprintf(fplog, "Initializing the distance RMSD potential\n");
    }
}

void calc_drmsd_pot(const gmx_multisim_t *ms,
                     int nfa, const t_iatom forceatoms[], const t_iparams ip[],
                     const rvec x[], const t_pbc *pbc,
                     t_fcdata *fcd, history_t *hist)
{

}

real ta_drmsd_pot(int nfa, const t_iatom forceatoms[], const t_iparams ip[],
               const rvec x[], rvec f[], rvec fshift[],
               const t_pbc *pbc, const t_graph *g,
               real lambda, real *dvdlambda,
               const t_mdatoms *md, t_fcdata *fcd,
               int *global_atom_index)
{
	return 0.0;
}
