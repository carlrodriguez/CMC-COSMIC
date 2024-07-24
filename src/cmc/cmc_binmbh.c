/* -*- linux-c -*- */
/* vi: set filetype=c.doxygen: */

#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include "cmc.h"
#include "cmc_vars.h"

/**
* @brief calculate the units used
*
* @param obj[2] ?
* @param bs_units ?
*/
void bmbh_calcunits(fb_obj_t *obj[2], fb_units_t *bs_units)
{
	bs_units->v = sqrt(FB_CONST_G*(obj[0]->m + obj[1]->m)/(obj[0]->m * obj[1]->m) * \
			(obj[1]->obj[0]->m * obj[1]->obj[1]->m / obj[1]->a));
	bs_units->l = obj[1]->a;
	bs_units->t = bs_units->l / bs_units->v;
	bs_units->m = bs_units->l * fb_sqr(bs_units->v) / FB_CONST_G;
	bs_units->E = bs_units->m * fb_sqr(bs_units->v);
}

/**
* @brief ?
*
* @param t ?
* @param ksin ?
* @param kbin ?
* @param W ?
* @param bmax ?
* @param hier ?
* @param rng ?
*
* @return ?
*/
fb_ret_t binmbh(double *t, long kbin, double v[3], double dist, fb_hier_t *hier, gsl_rng *rng, double time)
{
	int j;
	long jbin;
    int num_bh=0;
	double m0, m1, a1, e1, m10, m11;
	fb_units_t fb_units;
	fb_input_t input;
	fb_ret_t retval;
	
	/* a useful definition */
	jbin = star[kbin].binind;

	/* set parameters */
	input.ks = 0;
	input.Dflag = 0;
	input.dt = 0.0;
	input.tcpustop = 60.0;
	input.absacc = 1.0e-9;
	input.relacc = 1.0e-9;
	input.ncount = 500;
	input.tidaltol = 1.0e-5;
	input.speedtol = 1;

	input.PN1 = 0;
	input.PN2 = 0;
	input.PN25 = 0;
	input.PN3 = 0;
	input.PN35 = 0;
    /* Only used for stellar-mass BHs*/
    input.BH_REFF = BH_RADIUS_MULTIPLYER;
    input.BHNS_TDE_FLAG = BHNS_TDE;
	input.firstlogentry[0] = '\0';
	input.fexp = 1.0;
	fb_debug = 0;

	/* initialize a few things for integrator */
	*t = 0.0;
	hier->nstar = 3;
	fb_init_hier(hier);

	/* create binary */
	hier->hier[hier->hi[2]+0].obj[0] = &(hier->hier[hier->hi[1]+1]);
	hier->hier[hier->hi[2]+0].obj[1] = &(hier->hier[hier->hi[1]+2]);
	hier->hier[hier->hi[2]+0].t = *t;


	/* give the objects some properties */
	for (j=0; j<hier->nstar; j++) {
		hier->hier[hier->hi[1]+j].ncoll = 1;
		sprintf(hier->hier[hier->hi[1]+j].idstring, "%d", j);
		hier->hier[hier->hi[1]+j].n = 1;
		hier->hier[hier->hi[1]+j].obj[0] = NULL;
		hier->hier[hier->hi[1]+j].obj[1] = NULL;
		hier->hier[hier->hi[1]+j].Lint[0] = 0.0;
		hier->hier[hier->hi[1]+j].Lint[1] = 0.0;
		hier->hier[hier->hi[1]+j].Lint[2] = 0.0;
	}
	
	hier->hier[hier->hi[1]+0].id[0] = 0; // MBH is always star 0 
	hier->hier[hier->hi[1]+1].id[0] = binary[jbin].id1;
	hier->hier[hier->hi[1]+2].id[0] = binary[jbin].id2;

    /* Note: we may want TDEs without necessarily having star-star collisions*/
	if (SS_COLLISION || BH_LOSS_CONE) {
		hier->hier[hier->hi[1]+0].R = cenma.m*units.mstar / FB_CONST_MSUN * 9.8664506e-9 * FB_CONST_AU * BH_RADIUS_MULTIPLYER; // Convert mass to BH_radius
		hier->hier[hier->hi[1]+1].R = binary[jbin].rad1 * units.l;
		hier->hier[hier->hi[1]+2].R = binary[jbin].rad2 * units.l;
        if(binary[jbin].bse_kw[0] == 14) hier->hier[hier->hi[1]+1].R *= BH_RADIUS_MULTIPLYER; 
        if(binary[jbin].bse_kw[1] == 14) hier->hier[hier->hi[1]+2].R *= BH_RADIUS_MULTIPLYER; 
	} else {
		hier->hier[hier->hi[1]+0].R = 0.0;
		hier->hier[hier->hi[1]+1].R = 0.0;
		hier->hier[hier->hi[1]+2].R = 0.0;
	}

	hier->hier[hier->hi[1]+0].m = cenma.m * units.mstar;
	hier->hier[hier->hi[1]+1].m = binary[jbin].m1 * units.mstar;
	hier->hier[hier->hi[1]+2].m = binary[jbin].m2 * units.mstar;

	hier->hier[hier->hi[1]+0].k_type = 14;
	hier->hier[hier->hi[1]+1].k_type = binary[jbin].bse_kw[0];
	hier->hier[hier->hi[1]+2].k_type = binary[jbin].bse_kw[1];

	hier->hier[hier->hi[1]+0].chi = 0; 
	hier->hier[hier->hi[1]+1].chi = binary[jbin].bse_bhspin[0];
	hier->hier[hier->hi[1]+2].chi = binary[jbin].bse_bhspin[1];

	hier->hier[hier->hi[2]+0].m = hier->hier[hier->hi[1]+1].m + hier->hier[hier->hi[1]+2].m;

	hier->hier[hier->hi[1]+0].Eint = cenma.E * units.E;  
	hier->hier[hier->hi[1]+1].Eint = binary[jbin].Eint1 * units.E;
	hier->hier[hier->hi[1]+2].Eint = binary[jbin].Eint2 * units.E;

	hier->hier[hier->hi[2]+0].a = binary[jbin].a * units.l;
	hier->hier[hier->hi[2]+0].e = binary[jbin].e;

	hier->obj[0] = &(hier->hier[hier->hi[1]+0]);
	hier->obj[1] = &(hier->hier[hier->hi[2]+0]);
	hier->obj[2] = NULL;

	/* logging */
	parafprintf(binintfile, "********************************************************************************\n");
	parafprintf(binintfile, "type=BMBH t=%.9g\n", TotalTime);
	parafprintf(binintfile, "params[CMC Code Units]: distance=%g v_x=%g v_y=%g v_z=%g\n", dist, v[0], v[1], v[2]);
	/* set units to 1 since we're already in CGS */
	fb_units.v = fb_units.l = fb_units.t = fb_units.m = fb_units.E = 1.0;
	parafprintf(binintfile, "input: ");
	binint_log_obj(hier->obj[0], fb_units);
	parafprintf(binintfile, "input: ");
	binint_log_obj(hier->obj[1], fb_units);

	/* set the positions and velocities */
	hier->obj[0]->x[0] = 0.0; 
	hier->obj[0]->x[1] = 0.0; 
	hier->obj[0]->x[2] = 0.0;
	
	hier->obj[0]->v[0] = 0.0;
	hier->obj[0]->v[1] = 0.0; 
	hier->obj[0]->v[2] = 0.0;

	hier->obj[1]->x[0] = 0.0;
	hier->obj[1]->x[1] = 0.0;
	hier->obj[1]->x[2] = dist*units.l;
	
	hier->obj[1]->v[0] = v[0]*units.l/units.t;
	hier->obj[1]->v[1] = v[1]*units.l/units.t;
	hier->obj[1]->v[2] = v[2]*units.l/units.t;

	/* get the units and normalize */
	bmbh_calcunits(hier->obj, &fb_units);
	fb_normalize(hier, fb_units);

    /* This needs to be set here after fb_units is defined*/
	input.tstop = time*units.t/fb_units.t; 
	
	/* trickle down the binary properties, then back up */
	fb_randorient(&(hier->hier[hier->hi[2]+0]), rng, curr_st);
	fb_downsync(&(hier->hier[hier->hi[2]+0]), *t);
	fb_upsync(&(hier->hier[hier->hi[2]+0]), *t);
	
	/* call fewbody! */
	retval = fewbody(input, fb_units, hier, t, rng, curr_st);

	/* and return */
	return(retval);
}
