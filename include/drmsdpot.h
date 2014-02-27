/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by Manuel Luitz, Rainer Bomlies
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
#ifndef _distance_rmsd_pot_h
#define _distance_rmsd_pot_h
#include "typedefs.h"

/*
 * Implementation of the distance rmsd restraint potential function
 */

/*
 * Create header for drmsd output file
 */
void print_drmsd_header(FILE *fp, const output_env_t oenv, real lambda,
        real drmsd_ref, real f_const);

/*
 * Print the drmsd data to output file
 */
void print_drmsd_data(FILE *fp, real time, real drmsd, real vpot);

/*
 * Open a drmsd output file
 */
FILE *open_drmsd_out(const char *fn, const t_inputrec *ir,
                 const output_env_t oenv, gmx_bool append);


/*
 * Setup
 *
 * bIsREMD	Is this part of a replica exchange simulation
 */
void init_drmsd_pot(FILE *fplog, const gmx_mtop_t *mtop,
                 t_inputrec *ir, const t_commrec *cr, gmx_bool bPartDecomp,
                 t_fcdata *fcd, gmx_bool bIsREMD);

/*
 * Calculate the drmsd value of the current  step.
 * This has to be done before domain decomposition
 */
void calc_drmsd(int nfa, const t_iatom forceatoms[], const t_iparams ip[],
                     const rvec x[], const t_pbc *pbc,
                     t_fcdata *fcd, real lambda);
/*
 * Interaction function of the distance RMSD potential
 * This function may run on different threads and gets
 * only a subset of atom coordinates
 */
t_ifunc if_drmsd_pot;

#endif
