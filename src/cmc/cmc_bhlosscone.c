/* -*- linux-c -*- */
/* vi: set filetype=c.doxygen: */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include "cmc.h"
#include "cmc_vars.h"

/**
* @brief ?
*
* @param fname file name
*/
void create_rwalk_file(char *fname) {

    MPI_File mpi_rwalk_file;
    char mpi_rwalk_file_buf[10000];
    char mpi_rwalk_file_wrbuf[10000000];
    long long mpi_rwalk_file_len=0, mpi_rwalk_file_ofst_total=0;
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &mpi_rwalk_file);
	 if(tcount==1)
		 MPI_File_set_size(mpi_rwalk_file, 0);

  pararootfprintf(rwalk_file, "\n");
  pararootfprintf(rwalk_file, 
          "# 1:index, 2:Time, 3:r, 4:Trel, 5:dt, 6:l2_scale, 7:n_steps, 8:beta 9:n_local, 10:W, 11:P_orb, 12:n_orb\n");
  mpi_para_file_write(mpi_rwalk_file_wrbuf, &mpi_rwalk_file_len, &mpi_rwalk_file_ofst_total, &mpi_rwalk_file);
  MPI_File_close(&mpi_rwalk_file);
}

/**
* @brief ?
*
* @param fname ?
* @param index star index
* @param Trel ?
* @param dt ?
* @param l2_scale ?
* @param n_steps ?
* @param beta ?
* @param n_local ?
* @param W ?
* @param P_orb ?
* @param n_orb ?
*/
void write_rwalk_data(char *fname, long index, double Trel, double dt, 
    double l2_scale, double n_steps, double beta, double n_local, double W, 
    double P_orb, double n_orb) {

	double r = star_r[index];
    MPI_File mpi_rwalk_file;
    char mpi_rwalk_file_buf[10000];
    char mpi_rwalk_file_wrbuf[10000000];
    long long mpi_rwalk_file_len=0, mpi_rwalk_file_ofst_total=0;
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &mpi_rwalk_file);
	 if(tcount==1)
		 MPI_File_set_size(mpi_rwalk_file, 0);

  	parafprintf(rwalk_file, "%li %g %g %g %g %g %g %g %g %g %g %g\n", 
      index, TotalTime, r, Trel, dt, sqrt(l2_scale), n_steps, beta, n_local, W, P_orb, n_orb);

  	mpi_para_file_write(mpi_rwalk_file_wrbuf, &mpi_rwalk_file_len, &mpi_rwalk_file_ofst_total, &mpi_rwalk_file);
  	MPI_File_close(&mpi_rwalk_file);
}

/**
* @brief This is the random walk procedure as outlined by Freitag & Benz (2002). Change of notation: beta here is theta in the paper.
*
* @param index star index
* @param v[4] ?
* @param vcm[4] ?
* @param beta Deflection angle due to two-body relaxation in a given timestep.
* @param dt MC timestep
*/
void bh_rand_walk(long index, double v[4], double vcm[4], double beta, double dt, gsl_rng *rng)
{ 
	double w[3], n_orb, P_orb, deltabeta_orb, L2, Rdisr, Jlc, vlc;
	double deltamax, deltasafe, delta, dbeta;
	double w_mag, l2_scale;
	int i;
    int g_index;
	g_index = get_global_idx(index);

	char fname[80];
	long is_in_ids;
	double Trel, n_local, M2ave; 
	double W, n_steps= 1.;
    double t_conversion;
	
	is_in_ids= 0;
	sprintf(fname, "%s.rwalk_steps.dat", outprefix);
	n_local= calc_n_local(g_index, AVEKERNEL, clus.N_MAX);
	W = 4.0 * sigma_array.sigma[index] / sqrt(3.0*PI);
	M2ave= calc_average_mass_sqr(g_index, clus.N_MAX);
	Trel= (PI/32.)*cub(W)/ ( ((double) clus.N_STAR) * n_local * (4.0* M2ave) );
	if (index==1 && tcount%SNAPSHOT_DELTACOUNT==0 && SNAPSHOTTING && WRITE_RWALK_INFO) {
		is_in_ids=1;
		create_rwalk_file(fname);
	};
	/* simulate loss cone physics for central mass */
	//MPI: Parallelized, but might have mistakes since I am not clear as to what some functions are doing.
	
	/*First, calculate the orbital period */
	P_orb = calc_P_orb(index);	

	/*Then, calculate deltabeta_orb using eq. 30 in Freitag & Benz (2002). */
	n_orb = dt * ((double) clus.N_STAR)/log(GAMMA * ((double) clus.N_STAR)) / P_orb; 
	
	l2_scale= 1.;

	/* scale down L2 if the time step is larger than BH_LC_FDT*Trel */
	/* This is inconsistent, as for stars with dt< BH_LC_FDT*Trel the probability
	 * of hitting the loss cone becomes smaller, compared to the case with 
	 * dt=BH_LC_FDT*Trel
	 */
	/* if (BH_LC_FDT>0. && dt> BH_LC_FDT*Trel) { */
	if (BH_LC_FDT>0.) {
		n_steps= dt/BH_LC_FDT/Trel;
		l2_scale= 1./n_steps;
	};
	
	/*As a reminder, deltebeta_orb is deltatheta_orb in Freitag & Benz (2002).*/
	deltabeta_orb = 1.0/sqrt(n_orb) * sqrt(l2_scale)*beta;
	
	/*Define L2 as the total quadratic deflection angle during a timestep*/
	L2 = l2_scale*fb_sqr(beta);

    /*coefficient to get time from quadratic angle deflection below*/
    t_conversion = dt/L2;
	
	if (BH_R_DISRUPT_NB>0.) {/*The default value of BH_R_DISRUPT_NB is 0*/
		Rdisr= BH_R_DISRUPT_NB;
	} else if (STELLAR_EVOLUTION){
		double Rss;
		//dprintf("cenma.m= %g, star[%li].m= %g\n", cenma.m, index, star[index].m);
		if (star[index].binind > 0){
			double a = binary[star[index].binind].a;
			double e = binary[star[index].binind].e;
			Rdisr= pow(2.*cenma.m/star_m[g_index], 1./3.)*a*(1+e);
			//TODO: quick and cheap check for apoMBH of binary
		} else {
			Rdisr= pow(2.*cenma.m/star_m[g_index], 1./3.)*star[index].rad;
		}
		Rss=4.24e-06*cenma.m/SOLAR_MASS_DYN*RSUN/units.l; // Mass of MBH 
		Rdisr= MAX(Rdisr, Rss);
	} else {
		Rdisr= pow(2.*cenma.m/star_m[g_index], 1./3.)*star[index].rad;
	}

	/*Calculate the angular momentum at the LC using eq. 27 from Freitag & Benz (2002)*/
	Jlc= sqrt(2.*cenma.m*madhoc*Rdisr);
	vlc= Jlc/star_r[g_index];

	/*Calculate the particle's velocity in the encounter CM frame */
	for (i=0; i<3; i++) {
		w[i]= v[i+1]- vcm[i+1];
	}
	w_mag= sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
	delta= 0.0;
    int in_loss_cone=0;
    double time_for_binary;

	while (L2 > 0.0) { /* If L2 <= 0, the random walk is over*/
		L2 -= fb_sqr(delta); /*L2 is updated after the random walk*/

        /*If the tangential speed of the particle is less than vlc,the star has entered the loss cone and is disrupted */
        in_loss_cone = (sqrt(fb_sqr(w[0]+vcm[1])+fb_sqr(w[1]+vcm[2])) <= vlc);

		if (in_loss_cone){
			if(star[index].binind > 0) { //Binary

                /*Figure out how long to run fewbody (how long is this step?)*/
                time_for_binary = t_conversion*fb_sqr(delta);

                /*Then call fewbody*/
                fb_ret_t retval;
                double t=0;
                fb_hier_t hier;
                hier.nstarinit = 3;
                fb_malloc_hier(&hier);
                retval = binmbh(&t, index, w, star_r[g_index], &hier, rng, time_for_binary);

                /* free Fewbody memory */
                fb_free_hier(hier); // TODO: Sloppy to malloc and free this every iteration; only do it once

                /* Finally write to file */
                if(WRITE_BH_LOSSCONE_INFO)
                    parafprintf(bhlossconefile, "%g 1 %g %ld %ld %g %g %g %g %g %g %ld %ld %g %g %g %g %g %g %g %g \n", TotalTime, star[index].r, binary[star[index].binind].id1, binary[star[index].binind].id2, binary[star[index].binind].m1 * units.mstar / MSUN, binary[star[index].binind].m2 * units.mstar / MSUN,  binary[star[index].binind].rad1 * units.l / RSUN, binary[star[index].binind].rad2 * units.l / RSUN, binary[star[index].binind].bse_radc[0] * units.l / RSUN , binary[star[index].binind].bse_radc[1] * units.l / RSUN, binary[star[index].binind].bse_kw[0], binary[star[index].binind].bse_kw[1], binary[star[index].binind].a * units.l / AU, binary[star[index].binind].e, star[index].r_peri * units.l / AU, v[1], v[2], v[3], star[index].E, star[index].J);

			} else{ //Single
				cenma.m_new += star_m[g_index]; 
				//TODO: SMBH: this is just for bookkeeping (to track the energy deleted by destroying stars).  HOWEVER, the energy of the cluster also 
				//changes by virtue of the fact that you're increasing the SMBH mass.  It's possible we're double counting here.  
				cenma.E_new +=  (2.0*star_phi[g_index] + star[index].vr * star[index].vr + star[index].vt * star[index].vt) / 
				2.0 * star_m[g_index] * madhoc;

                /* Destroy the star and complete the random walk */
                destroy_obj(index);
                L2 = 0.0; 

                /* Finally write to file */
                if (WRITE_BH_LOSSCONE_INFO) 
                    parafprintf(bhlossconefile, "%g 0 %g %ld -100 %g -100 %g -100 %g -100 %ld -100 -100 -100 %g %g %g %g %g %g\n", TotalTime, star[index].r, star[index].id, star[index].m * units.mstar / MSUN, star[index].rad  * units.l / RSUN, star[index].se_rc * units.l / RSUN, star[index].se_k, star[index].r_peri * units.l / AU, v[1], v[2], v[3], star[index].E, star[index].J );
			}
		}

        /*The amplitude of the random walk step is calculated using eq. 31 of Freitag & Benz (2002).*/
        deltamax= 0.1*FB_CONST_PI;
        deltasafe= CSAFE*(sqrt(fb_sqr(w[0]+vcm[1])+fb_sqr(vcm[2]+w[1]))-vlc)/w_mag;
        delta = MAX(deltabeta_orb, MIN(deltamax, MIN(deltasafe, sqrt(L2)))); 

        /* If we have a binary in the loss cone then take steps on order of p_orb */
        if(in_loss_cone) delta = deltabeta_orb;

        /*Set the direction of the random walk step by drawing a random angle dbeta*/
        dbeta = 2.0 * PI * rng_t113_dbl_new(curr_st); 

        do_random_step(w, dbeta, delta); 
	} 

	if (tcount%SNAPSHOT_DELTACOUNT==0 && SNAPSHOTTING && WRITE_RWALK_INFO) {
		write_rwalk_data(fname, g_index, Trel, dt, l2_scale, n_steps, beta,
				n_local, W, P_orb, n_orb);
	}

	/*Free up the star structure since we have already r_peri and r_apo in the bhlosscone file */
	star[index].r_peri = 0.0;
	star[index].r_apo = 0.0;

};

/* here the notation of Freitag & Benz (2002) is used */
/**
* @brief ?
*
* @param w Star's velocity vector 
* @param beta Random angle 
* @param delta Amplitude of the random walk step
*/
void do_random_step(double *w, double beta, double delta) {
   double theta, phi, w_mag, new_w_dir[3];
   double a; /* this is a variable to store an intermediate result*/ 

   /* calculate direction of w*/
   w_mag= sqrt(fb_sqr(w[0])+ fb_sqr(w[1])+ fb_sqr(w[2]));
   theta= acos(w[2]/w_mag);
   phi= atan2(w[1], w[0]);
   
   /* rotate new vector (w_mag, beta, delta) into the direction of w */
   a= cos(theta)* sin(beta)* sin(delta)+ sin(theta)* cos(delta);
   new_w_dir[0]= sin(phi)* cos(beta) * sin(delta) + cos(phi)* a;
   new_w_dir[1]= -cos(phi)* cos(beta)* sin(delta) + sin(phi)* a;
   new_w_dir[2]= -sin(theta)* sin(beta)* sin(delta) + cos(theta)* cos(delta);

   w[0]= w_mag* new_w_dir[0];
   w[1]= w_mag* new_w_dir[1];
   w[2]= w_mag* new_w_dir[2];
};

/**
* @brief calculate star's radial orbital period
*
* @param index star index
*
* @return star's radial orbital period
*/
double calc_P_orb(long index)
{
	double E, J, Porb, error, Porbapproxmin, Porbapproxmax, Porbapprox;
	orbit_rs_t orbit_rs;
	calc_p_orb_params_t params;
	gsl_function F;
	struct Interval star_interval;
	int g_index;
	g_index = get_global_idx(index);


    /* default values for star_interval */
	star_interval.min= 1;
	star_interval.max= clus.N_MAX+1;

	E = star[index].E + MPI_PHI_S(star_r[g_index], g_index);
	J = star[index].J;

	/*TODO: we're doing this twice?!  Put in a flag to only compute
	 * new orbital parameters here if we're using the loss cone*/
	/*Find out new orbit of star, primarily peri and apo- center distances*/
	orbit_rs = calc_orbit_new(index, E, J);

	/*Assigning the values of r_peri and r_apo to the star structure so they can be accessed for the bhlosscone file later*/
	star[index].r_peri = orbit_rs.rp;
	star[index].r_apo = orbit_rs.ra;

	/*A cheap way to compute the radial orbital period: just 
	 * geometrically average the equivalent orbital period of a
	 * circular orbit at pericenter and apocenter.  This works exactly
	 * in the Keplarian case, but not for a general spherical potential*/
	Porbapproxmin = 2.0 * PI * orbit_rs.rp / (J / orbit_rs.rp);
	Porbapproxmax = 2.0 * PI * orbit_rs.ra / (J / orbit_rs.ra);
	Porbapprox = sqrt(Porbapproxmin * Porbapproxmax);

	if (orbit_rs.ra-orbit_rs.rp< CIRC_PERIOD_THRESHOLD) {
		dprintf("Orbit is considered circular for period calculation.\n");
		dprintf("ra-rp= %g and is less than the threshold %g\n", 
			orbit_rs.ra-orbit_rs.rp, CIRC_PERIOD_THRESHOLD);
		orbit_rs.circular_flag = 1;
	}

	if (orbit_rs.circular_flag == 1) {
		/* We're returning the azimuthal period here, which is not the same as the
		   radial period for the general cluster potential.  This shouldn't make
		   any difference for the BH loss cone stuff, since the orbit is circular. */
		return(2.0 * PI * star_r[g_index] / star[index].vt);
	} else {

		params.E = E;
		params.J = J;
		params.index = g_index;
		params.kmax= orbit_rs.kmax+1; 
		params.kmin= orbit_rs.kmin;
		params.rp = orbit_rs.rp;
		params.ra = orbit_rs.ra;
		F.params = &params;
	 
        /* test if the interval of rmax is the same as [kmax,kmax+1] */
		if (!orbit_rs.circular_flag) {
			if (function_Q(g_index, params.kmax, E, J)>0 || function_Q(g_index, params.kmax-1,E, J)<0) {
			  dprintf("r and phi interval do not match: id= %li, r_kmax= %li\n",
			    star[index].id, params.kmax-1);
			  dprintf("star index is %li\n", index);
			  dprintf("f_Q[r_kmax]= %g; f_Q[r_kmax+1]= %g\n", function_Q(g_index, params.kmax-1, E, J),
			    function_Q(g_index, params.kmax, E, J));
			  dprintf("phi_kmax= %li, phi_kmax+1= %li\n", orbit_rs.kmax, orbit_rs.kmax+1);
			  dprintf("f_Q[phi_kmax]= %g; f_Q[phi_kmax+1]= %g\n", 
			    function_Q(g_index, orbit_rs.kmax, E, J),
			    function_Q(g_index, orbit_rs.kmax+1, E, J));
			  dprintf("(r[r_kmax]-rmax=%g, r[r_kmax+1]-rmax)= %g\n", 
			    star[params.kmax-1].r-orbit_rs.ra,star[params.kmax].r-orbit_rs.ra);
			};
                        if (calc_vr(params.rp, index, E, J)< 0.) {
                          dprintf("Harrrg: vr(rmin)< 0.! Damn it! Index: %li, Id: %li\n", index, 
                            star[index].id);
                        };
                        if (calc_vr(params.ra, index, E, J)< 0.) {
                          dprintf("Harrrg: vr(rmax)< 0.! Damn it! Index: %li, Id: %li\n", index, 
                            star[index].id);
                        };
                        if (params.kmax!=orbit_rs.kmax+1) 
                          dprintf("kmax in orbit_rs and params differ! kmax_o= %li, kmax_p=%li, Index: %li, Id: %li\n", 
                            orbit_rs.kmax, params.kmax, index, star[index].id);
                        if ((params.kmin!=orbit_rs.kmin)&& (params.kmin>1)) 
                          dprintf("kmin in orbit_rs and params differ! kmin_o= %li, kmin_p=%li, Index: %li, Id: %li\n", 
                            orbit_rs.kmin, params.kmin, index, star[index].id);
		};

		// Cabrera 230731: Added error handling to default to Porbapprox
		gsl_error_handler_t *old_handler;
		old_handler = gsl_set_error_handler_off();

        //MPI: These seem to be never used, so they are not parallelized as of now. Also unclear how these functions will be used, so dont understand if index transformation is reqd or not. Might later have to double check if the parallelization is correct.
		if (0) { /* use standard potential function with Stefan's speedup trick here */
			F.function = &calc_p_orb_f;
			status = gsl_integration_qags(&F, orbit_rs.rp, orbit_rs.ra, 0, 1.0e-3, 1000, workspace_lc_porb_integral, &Porb, &error);
		}

		if (0) { /* use fast potential function here (not much of a speedup over Stefan's technique in practice) */
			F.function = &calc_p_orb_f2;
			status = gsl_integration_qags(&F, orbit_rs.rp, orbit_rs.ra, 0, 1.0e-3, 1000, workspace_lc_porb_integral, &Porb, &error);
		}

		if (1) { /* use Gauss-Chebyshev for factor of ~few speedup over standard method */
			F.function = &calc_p_orb_gc;
			status = gsl_integration_qaws(&F, orbit_rs.rp, orbit_rs.ra, table_lc_porb_integral,
				       	1.0e-3, 1.0e-3, 1000, workspace_lc_porb_integral, &Porb, &error);

		}

		// Print error if there's a problem with the integration, and return Porbapprox
		if (status) {
			eprintf("gsl_integration_qa[g,w]s failed (gsl_errno=%d, index=%d, g_index=%d, orbit_rs.rp=%e, orbit_rs.ra=%e); check cmc_bhlosscone.c for [g,w] routine; returning Porbapprox=%e\n", status, index, g_index, orbit_rs.rp, orbit_rs.ra, Porbapprox);
			Porb = Porbapprox;
		}

		// Reset to previous error handler
		gsl_set_error_handler(old_handler);
		
		return(Porb);
	}
}

/**
* @brief integrand for calc_P_orb
*
* @param x Variable of integration
* @param params Parameters of the orbit
*
* @return ?
*/
double calc_p_orb_f(double x, void *params) {
	calc_p_orb_params_t myparams = *(calc_p_orb_params_t *) params;
	double radicand;

	radicand = 2.0 * myparams.E - fb_sqr(myparams.J/x) - 2.0 * (potential(x) + MPI_PHI_S(x, myparams.index));

	if (radicand < 0.0) {
		dprintf("radicand=%g<0; setting to zero; index=%ld\n", radicand, myparams.index);
		radicand = 0.0;
	}
	
	return(2.0 / sqrt(radicand));
}

/**
* @brief integrand for calc_P_orb
*
* @param x Variable of integration
* @param params Parameters of the orbit
*
* @return ?
*/
double calc_p_orb_f2(double x, void *params) {
	calc_p_orb_params_t myparams = *(calc_p_orb_params_t *) params;
	double radicand;

	radicand = 2.0 * myparams.E - fb_sqr(myparams.J/x) - 2.0 * (fastpotential(x, myparams.kmin, myparams.kmax) + MPI_PHI_S(x, myparams.index));

	if (radicand < 0.0) {
		dprintf("radicand=%g<0; setting to zero; index=%ld\n", radicand, myparams.index);
		radicand = 0.0;
	}
	
	return(2.0 / sqrt(radicand));
}

/**
* @brief integrand for calc_P_orb, using Gauss-Chebyshev for regularizing the integrand near the endpoints
*
* @param x Variable of integration
* @param params Parameters of the orbit
*
* @return ?
*/
double calc_p_orb_gc(double x, void *params) {
	calc_p_orb_params_t myparams = *(calc_p_orb_params_t *) params;
	double radicand;
	double phik, phik1, phi0, phi1, rk, rk1, rminus, rplus;
	double E, J, rp, ra;
	double result;
	long index, kmin, kmax;

	E = myparams.E;
	J = myparams.J;
	index = myparams.index;
	kmin = myparams.kmin;
	kmax = myparams.kmax;
	rp = myparams.rp;
	ra = myparams.ra;

	/*kmin and kmax are the global indeces that enclose the position of a particle. 
	Our star is located between star_r[kmin] and star_r[kmax]*/

	if (x <= star_r[kmin+1]) { 
		/* return integrand regularized at r=rp */
		phik = star_phi[kmin] + MPI_PHI_S(star_r[kmin], index);
		phik1 = star_phi[kmin+1] + MPI_PHI_S(star_r[kmin+1], index);
		rk = star_r[kmin];
		rk1 = star_r[kmin+1];
		
		/* The expression for the potential at any point r along the orbit is analogous to eq. 33 in Rodriguez et al. (2021) */
		/* Note:  phi(r) = phi0 + phi1/r */
		phi0 = phik + (phik1 - phik)/(1.0-rk/rk1); /*Interpolated potential between r[kmin] and r[kmin+1]*/
		phi1 = (phik - phik1)/(1.0/rk - 1.0/rk1); 

		if (E-phi0==0.) {
		  dprintf("E is phi0 near rp! Damn it!");
		};

		/*The roots for the energy equation (eq. 26 in Rodriguez et al. (2021)) using the definition for the potential at a point r:  phi(r) = phi0 + phi1/r */
		rminus = (phi1 - sqrt(fb_sqr(phi1)+2.0*fb_sqr(J)*(E-phi0))) / (2.0*(E-phi0));
		if (kmax == kmin + 1) {			
			/* then rminus = ra, so must cancel (ra-x)/(rminus-x) term analytically */
			result= 2.0*x*sqrt(1.0/((2.0*phi0-2.0*E)));
			if (gsl_isinf(result)) {
			  dprintf("result is infinite near rp! Damn it! kmax==kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			};
			if (gsl_isnan(result)) {
			  dprintf("result is NaN near rp! Damn it! kmax==kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			};
			return(2.0*x*sqrt(1.0/((2.0*phi0-2.0*E))));
		} else {
		        result= 2.0*x*sqrt((ra-x)/((2.0*phi0-2.0*E)*(rminus-x)));
			if (gsl_isinf(result)) {
			  dprintf("result is infinite near rp! Damn it! kmax!=kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			};
			if (gsl_isnan(result)) {
			  dprintf("result is NaN near rp! Damn it! kmax!=kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			};
			return(2.0*x*sqrt((ra-x)/((2.0*phi0-2.0*E)*(rminus-x))));
		}
	} else if (x >= star_r[kmax-1]) { 
		/* return integrand regularized at r=ra*/
		phik = star_phi[kmax-1] + MPI_PHI_S(star_r[kmax-1], index);
		phik1 = star_phi[kmax] + MPI_PHI_S(star_r[kmax], index);
		rk = star_r[kmax-1];
		rk1 = star_r[kmax];
		phi0 = phik + (phik1 - phik)/(1.0-rk/rk1);
		phi1 = (phik - phik1)/(1.0/rk - 1.0/rk1);

		rplus = (phi1 + sqrt(fb_sqr(phi1)+2.0*fb_sqr(J)*(E-phi0))) / (2.0*(E-phi0));
		if (kmax == kmin + 1) {
			result=2.0*x*sqrt(1.0/((2.0*phi0-2.0*E)));
			/* then rplus = rp, so must cancel (x-rp)/(x-rplus) term analytically */
			if (gsl_isinf(result)) {
			  dprintf("result is infinite near ra! Damn it! kmax==kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			  dprintf("x=%g, ra= %g, rp= %g\n", x, ra, rp);
			};
			if (gsl_isnan(result)) {
			  dprintf("result is NaN near ra! Damn it! kmax==kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			  dprintf("x=%g, ra= %g, rp= %g\n", x, ra, rp);
			};
			return(2.0*x*sqrt(1.0/((2.0*phi0-2.0*E))));
		} else {
			result= 2.0*x*sqrt((x-rp)/((2.0*phi0-2.0*E)*(x-rplus)));
			if (gsl_isinf(result)) {
			  dprintf("result is infinite near ra! Damn it! kmax!=kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			  dprintf("x=%g, ra= %g, rp= %g\n", x, ra, rp);
			};
			if (gsl_isnan(result)) {
			  dprintf("result is NaN near ra! Damn it! kmax==kmin+1\n");
			  dprintf("kmax=%li, kmin=%li, index=%li, rk=%g, rk1=%g, phi0=%g, phi1=%g\n",
			    kmax, kmin, index, rk, rk1, phi0, phi1);
			  dprintf("phik= %g, phik1= %g, E= %g, J=%g, id=%li\n", phik, phik1, E, J, star[index].id);
			  dprintf("x=%g, ra= %g, rp= %g\n", x, ra, rp);
			};
			return(2.0*x*sqrt((x-rp)/((2.0*phi0-2.0*E)*(x-rplus))));
		}
	} else {
		radicand = 2.0 * (E - (potential(x) + MPI_PHI_S(x, index)))- fb_sqr(J/x);
		if (radicand < 0.0) {
			dprintf("radicand=%g<0; setting to zero; index=%ld\n", radicand, index);
			dprintf("kmin= %li, kmax= %li, rp=%g, ra=%g, Id: %li\n",
			  kmin, kmax, rp, ra, star[index].id);
			radicand = 0.0;
		};
		result= 2.0 * sqrt((x-rp)*(ra-x)/radicand);
		if (gsl_isinf(result)) {
		  dprintf("result is infinite! Damn it!\n");
		  dprintf("kmax=%li, kmin=%li, index=%li\n",
		    kmax, kmin, index);
		  dprintf("E= %g, J=%g, id=%li\n", E, J, star[index].id);
  	          dprintf("x=%g, ra= %g, rp= %g\n", x, ra, rp);
                  dprintf("x-rp= %g, ra-x= %g, rp-r[kmin+1]= %g\n", x-rp, ra-x, rp-star_r[kmin+1]);
                  dprintf("ra-r[kmax-1]= %g\n", ra-star_r[kmax-1]);
		};
		if (gsl_isnan(result)) {
		  dprintf("result is NaN! Damn it!\n");
		  dprintf("kmax=%li, kmin=%li, index=%li\n",
		    kmax, kmin, index);
		  dprintf("E= %g, J=%g, id=%li\n", E, J, star[index].id);
  	          dprintf("x=%g, ra= %g, rp= %g, radicand= %g\n", x, ra, rp, radicand);
                  dprintf("x-rp= %g, ra-x= %g, rp-r[kmin+1]= %g\n", x-rp, ra-x, rp-star_r[kmin+1]);
                  dprintf("ra-r[kmax-1]= %g\n", ra-star_r[kmax-1]);
		};
		return(2.0 * sqrt((x-rp)*(ra-x)/radicand));
	}
}

/**
* @brief ?
*
* @param r ?
*
* @return ?
*/
struct Interval get_r_interval(double r) {
  long kmax, kmin, i;
  struct Interval star_interval;

  if (SEARCH_GRID) {
   star_interval= search_grid_get_interval(r_grid, r);
   kmax= star_interval.max;
   kmin= star_interval.min;
  } else {
   kmax= clus.N_MAX+1;
   kmin= 1;
  };
  if (kmin==kmax-1) {
   i= kmin;
  } else {
   i =  FindZero_r(kmin, kmax, r);
  };

  star_interval.max= i+1;
  star_interval.min= i;
  return (star_interval);
}

int analyze_fewbody_output(fb_hier_t *hier, fb_ret_t *retval, long index){
    
    int mbhid=0,binid=0,sinid=0;
    /* One object -- either double TDE or the binary is unchanged 
     * (and is technically in a triple with the MBH)*/
        // Remember nobj is number of objects directly below this object, 
        // n is the total number of stars under this (all things under the hierarchy)
	star_t tempstar1, tempstar2;
    int i;
    double vs[20];
    long knew, knew1, knew2;
    double r_imbh_frame[3], v_imbh_frame[3], rhat[3], v_t[3];
    double rmag, v_rmag;
	char string1[1024], string2[1024];

	fb_units_t cmc_units, printing_units;
    cmc_units.v = sqrt((star_m[0]+star_m[get_global_idx(index)])/(star_m[0]*star_m[get_global_idx(index)]) * 
               (binary[star[index].binind].m1 * binary[star[index].binind].m2 / binary[star[index].binind].a) * madhoc);
    cmc_units.l = binary[star[index].binind].a;
    cmc_units.t = cmc_units.l / cmc_units.v;
    cmc_units.m = cmc_units.l * sqr(cmc_units.v);
    cmc_units.E = cmc_units.m * sqr(cmc_units.v);

	printing_units.v = cmc_units.v * units.l / units.t;
	printing_units.l = cmc_units.l * units.l;
	printing_units.t = cmc_units.t * units.t;
	printing_units.m = cmc_units.m * units.m;
	printing_units.E = cmc_units.E * units.E;

    /* First check if the integration actually worked; return -1 if error */
	if ( !( (fabs(retval->DeltaEfrac) < 1.0e-3 || fabs(retval->DeltaE) < 1.0e-3) && 
		 (fabs(retval->DeltaLfrac) < 1.0e-3 || fabs(retval->DeltaL) < 1.0e-3) ) && 
         (!((fabs(retval->DeltaE_GWfrac) > 1.0e-3 || fabs(retval->DeltaE_GW > 1.0e-3)
            ) && retval->PN_ON == 1))) /* did we have a significant energy error that wasn't from gravitational waves? */
    {
		parafprintf(binintfile, "outcome: energy and/or angular momentum error\n");
		print_interaction_error();
        return -1;
	} else if ( isnan(retval->DeltaE) || isnan(retval->DeltaL) ) {
		parafprintf(binintfile, "outcome: NaN returned by fewbody\n");
		print_interaction_error();
        return -1;
	} else if (retval->retval == 0) {
		/* bad outcome; ignore for now */
		parafprintf(binintfile, "outcome: stopped\n");
		print_interaction_error();
        return -1;
	} else if (hier->obj[0]->n == 4) {
		/* outcome is a quadruple */
		parafprintf(binintfile, "outcome: error\n");
		print_interaction_error();
        return -1;
	} else 
		parafprintf(binintfile, "outcome: %s (%s)\n", fb_sprint_hier(*hier, string1), fb_sprint_hier_hr(*hier, string2));

    /* Five Cases -- 
     *  --Binary unchanged  return 0
     *  --Binary disrupted  return 1
     *  --One TDE           return 2
     *  --Two TDEs          return 3
     *  --Binary merger     return 4
     *  
     *  Given this, it's easier to look for specific cases rather than generically process the full output */

    // ONE BIG TODO:
    // we're not differentiating between binary merger -> TDE versus both TDEs...

    double MBH_TDE_ACCRETION = 0.5;//
    // TODO: make this a parameter and set it *everywhere* (including inside fewbody)

    if (hier->nobj == 1){ /* One top-level object */
        if (hier->obj[0]->n == 1){ /* Only one star; either binary merger that's TDEd or double TDE*/
            /*Add mass of binary to the MBH*/
            cenma.m_new += MBH_TDE_ACCRETION*star_m[get_global_idx(index)]; 

            /*Energy too*/
            cenma.E_new +=  (2.0*star_phi[get_global_idx(index)] + star[index].vr * star[index].vr + star[index].vt * star[index].vt) / 
            2.0 * star_m[get_global_idx(index)] * madhoc + star[index].Eint; /*I don't think the binary should have any internal energy, but better safe than sorry*/

            /*Don't forget binding energy*/
            cenma.E_new -= binary[star[index].binind].m1 * binary[star[index].binind].m2 * sqr(madhoc) 
                / (2.0 * binary[star[index].binind].a) 
                - binary[star[index].binind].Eint1 - binary[star[index].binind].Eint2;

            /*Carl: TODO: double check minus sign there.
             * Also, double check whether we need MBH_TDE_ACCRETION for cenma energy as well*/

            /* Destroy the binary and complete the random walk */
            destroy_obj(index);
            return 3;

        } else if (hier->obj[0]->n == 2){ /*Or maybe the binary merged or one TDE*/
            if ((hier->obj[0]->obj[0]->id[0] == 0) && (hier->obj[0]->obj[0]->ncoll == 1)){
                binid = 1;
            } else if ((hier->obj[0]->obj[1]->id[0] == 0) && (hier->obj[0]->obj[1]->ncoll == 1)){
                mbhid = 1;
            }
            if ((mbhid == 1) || (binid == 1)){ /*IMBH is unmerged, must be a binary merger*/

                /* Create a new star for the binary merger*/
				knew = create_star(index, 0);

                /*Extract the position/velocity wrt the IMBH (assumed cluster center) from fewbody*/
                for(i=0; i<3; i++){
                    r_imbh_frame[i] = hier->obj[0]->obj[binid]->x[i] - hier->obj[0]->obj[mbhid]->x[i];
                    v_imbh_frame[i] = hier->obj[0]->obj[binid]->v[i] - hier->obj[0]->obj[mbhid]->v[i];
                }

                rmag = fb_mod(r_imbh_frame);
                for(i=0; i<3; i++) rhat[i] /= rmag;

                /*Set radial position*/
                star_r[get_global_idx(knew)] = rmag*cmc_units.l; 

                /*Set velocities as well*/
                v_rmag = fb_dot(rhat,v_imbh_frame);
                for(i=0; i<3; i++) v_t[i] = v_imbh_frame[i] - rhat[i]*v_rmag;

                star[knew].vr = v_rmag*cmc_units.v;
                star[knew].vt = fb_mod(v_t)*cmc_units.v;

                /*Set mass; this gets overwritten below*/
                star_m[get_global_idx(knew)] = hier->obj[0]->obj[binid]->m * cmc_units.m/madhoc;

                /*set potential*/
                star_phi[get_global_idx(knew)] = potential(star_r[get_global_idx(knew)]);

                /* Calculate new energies by recomputing E = PE + KE using new velocity */
                set_star_EJ(knew);

                /* set rnew, vrnew, vtnew */
                set_star_news(knew);
                
                /* I don't even know if this is necessary */
                set_star_olds(knew);

                /* mark stars as interacted */
                star[knew].interacted = 1;
                
                /* Finally, actually merge the two stars using COSMIC */
                cp_SEvars_to_star(index, 0, &tempstar1);
                cp_m_to_star(index, 0, &tempstar1);
                cp_SEvars_to_star(index, 1, &tempstar2);
                cp_m_to_star(index, 1, &tempstar2);
                merge_two_stars(&tempstar1, &tempstar2, &(star[knew]), vs, curr_st);
                star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
                vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);

                /* Destroy the original binary */
                destroy_obj(index);

                return 4;
            } else { /*If not, then single TDE*/ 
                /*find the unTDE'd star*/
                if(hier->obj[0]->obj[0]->ncoll > 1)
                    sinid = 1;
                else 
                    mbhid = 1;

				knew = create_star(index, 0);
                star[knew].id = hier->obj[0]->obj[sinid]->id[0];

                /*Extract the position/velocity wrt the IMBH (assumed cluster center) from fewbody*/
                for(i=0; i<3; i++){
                    r_imbh_frame[i] = hier->obj[0]->obj[sinid]->x[i] - hier->obj[0]->obj[mbhid]->x[i];
                    v_imbh_frame[i] = hier->obj[0]->obj[sinid]->v[i] - hier->obj[0]->obj[mbhid]->v[i];
                }

                rmag = fb_mod(r_imbh_frame);
                for(i=0; i<3; i++) rhat[i] /= rmag;

                /*Set radial position*/
                star_r[get_global_idx(knew)] = rmag*cmc_units.l; 

                /*Set velocities as well*/
                v_rmag = fb_dot(rhat,v_imbh_frame);
                for(i=0; i<3; i++) v_t[i] = v_imbh_frame[i] - rhat[i]*v_rmag;

                star[knew].vr = v_rmag*cmc_units.v;
                star[knew].vt = fb_mod(v_t)*cmc_units.v;

                /*Set mass*/
                star_m[get_global_idx(knew)] = hier->obj[0]->obj[binid]->m * cmc_units.m/madhoc;

                /*set potential*/
                star_phi[get_global_idx(knew)] = potential(star_r[get_global_idx(knew)]);

                /* Calculate new energies by recomputing E = PE + KE using new velocity */
                set_star_EJ(knew);

                /* set rnew, vrnew, vtnew */
                set_star_news(knew);
                
                /* I don't even know if this is necessary */
                set_star_olds(knew);

                /* mark stars as interacted */
                star[knew].interacted = 1;

                /* Copy stellar evolution parameters from binary member */
                if(hier->obj[0]->obj[sinid]->id[0] == binary[star[index].binind].id1)
                    binid = 0;
                else
                    binid = 1;
                cp_SEvars_to_newstar(index, binid, knew);

                /*Add mass of other star to the MBH*/
                cenma.m_new += hier->obj[0]->obj[mbhid]->m * cmc_units.m/madhoc; 

                /*Energy too*/
                cenma.E_new += hier->obj[0]->obj[mbhid]->Eint * cmc_units.E; 
                /*TODO: Note: this is not going to conserve the energy correctly, largely because
                 * we're doing the encounter in a vacuum (i.e. not the cluster potential), and making
                 * that transition already technicaly breaks energy conservation (since the binding
                 * energy of the outer particle is different in fewbody vs CMC).  But then again,
                 * accreting any fraction of the star other than 100% also breaks energy conservation
                 * since we're not tracking the gas...
                 *
                 * For now, it's the best we can do*/

                /* Destroy the original binary */
                destroy_obj(index);

                return 2; 
            }
        } else { /* If three objects then it's a triple!*/

            /*first find the binary*/
            if(hier->obj[0]->obj[0]->n == 1)
                binid = 1;
            else
                sinid = 1;

            if(hier->obj[0]->obj[sinid]->id[0] == 0) /*if the single is the MBH, then the binary is unchanged*/
                return 0;

            if(hier->obj[0]->obj[binid]->obj[1]->id[0] == 0)
                mbhid = 1;

            /* Create both stars and re-insert into cluster */
            knew1 = create_star(index, 0);
            knew2 = create_star(index, 0);
            star[knew1].id = hier->obj[0]->obj[sinid]->id[0];
            star[knew2].id = hier->obj[0]->obj[binid]->obj[1-mbhid]->id[0]; /*1-mbhid gives the id of the star bound to the MBH*/

            /* FIRST STAR (tertiary of triple) */
            /* Extract the position/velocity wrt the IMBH (assumed cluster center) from fewbody*/
            for(i=0; i<3; i++){
                r_imbh_frame[i] = hier->obj[0]->obj[sinid]->x[i] - hier->obj[0]->obj[binid]->obj[mbhid]->x[i];
                v_imbh_frame[i] = hier->obj[0]->obj[sinid]->v[i] - hier->obj[0]->obj[binid]->obj[mbhid]->v[i];
            }

            rmag = fb_mod(r_imbh_frame);
            for(i=0; i<3; i++) rhat[i] /= rmag;

            /*Set radial position*/
            star_r[get_global_idx(knew1)] = rmag*cmc_units.l; 

            /*Set velocities as well*/
            v_rmag = fb_dot(rhat,v_imbh_frame);
            for(i=0; i<3; i++) v_t[i] = v_imbh_frame[i] - rhat[i]*v_rmag;

            star[knew1].vr = v_rmag*cmc_units.v;
            star[knew1].vt = fb_mod(v_t)*cmc_units.v;

            /* SECOND STAR (inner binary star) */
            /* Extract the position/velocity wrt the IMBH (assumed cluster center) from fewbody*/
            for(i=0; i<3; i++){
                r_imbh_frame[i] = hier->obj[0]->obj[binid]->obj[1-mbhid]->x[i] - hier->obj[0]->obj[binid]->obj[mbhid]->x[i];
                v_imbh_frame[i] = hier->obj[0]->obj[binid]->obj[1-mbhid]->v[i] - hier->obj[0]->obj[binid]->obj[mbhid]->v[i];
            }

            rmag = fb_mod(r_imbh_frame);
            for(i=0; i<3; i++) rhat[i] /= rmag;

            /*Set radial position*/
            star_r[get_global_idx(knew2)] = rmag*cmc_units.l; 

            /*Set velocities as well*/
            v_rmag = fb_dot(rhat,v_imbh_frame);
            for(i=0; i<3; i++) v_t[i] = v_imbh_frame[i] - rhat[i]*v_rmag;

            star[knew2].vr = v_rmag*cmc_units.v;
            star[knew2].vt = fb_mod(v_t)*cmc_units.v;

            /*Set mass; this gets overwritten below*/
            star_m[get_global_idx(knew1)] = hier->obj[0]->obj[sinid]->m * cmc_units.m/madhoc;
            star_m[get_global_idx(knew2)] = hier->obj[0]->obj[binid]->obj[1-mbhid]->m * cmc_units.m/madhoc;

            /*set potential*/
            star_phi[get_global_idx(knew1)] = potential(star_r[get_global_idx(knew1)]);
            star_phi[get_global_idx(knew2)] = potential(star_r[get_global_idx(knew2)]);

            /* Calculate new energies by recomputing E = PE + KE using new velocity */
            set_star_EJ(knew1);
            set_star_EJ(knew2);

            /* set rnew, vrnew, vtnew */
            set_star_news(knew1);
            set_star_news(knew2);
            
            /* I don't even know if this is necessary */
            set_star_olds(knew1);
            set_star_olds(knew2);

            /* mark stars as interacted */
            star[knew1].interacted = 1;
            star[knew2].interacted = 1;

            /* Copy stellar evolution parameters from binary member */
            /* reusing binid here*/
            if(hier->obj[0]->obj[sinid]->id[0] == binary[star[index].binind].id1)
                binid = 0;
            else
                binid = 1;
            cp_SEvars_to_newstar(index, binid, knew1);
            cp_SEvars_to_newstar(index, 1-binid, knew2);

            /* Destroy the original binary */
            destroy_obj(index);

            return 1;
        }
    } else if (hier->nobj == 2){ /* Two top-level objects */ 
        if (hier->nstar == 3){ /*We have a binary and single; same code as above*/

            if(hier->obj[0]->n == 1) /*first find the binary*/
                binid = 1;
            else 
                sinid = 1;

            if (hier->obj[sinid]->id[0] == 0) /* If the single is the MBH, the binary was unchanged*/
                return 0;

            /* Otherwise we have an exchange*/
            if(hier->obj[binid]->obj[1]->id[0] == 0)
                mbhid = 1;

            /* Create both stars and re-insert into cluster */
            knew1 = create_star(index, 0);
            knew2 = create_star(index, 0);
            star[knew1].id = hier->obj[sinid]->id[0];
            star[knew2].id = hier->obj[binid]->obj[1-mbhid]->id[0]; /*1-mbhid gives the id of the star bound to the MBH*/

            /* FIRST STAR (tertiary of triple) */
            /* Extract the position/velocity wrt the IMBH (assumed cluster center) from fewbody*/
            for(i=0; i<3; i++){
                r_imbh_frame[i] = hier->obj[sinid]->x[i] - hier->obj[binid]->obj[mbhid]->x[i];
                v_imbh_frame[i] = hier->obj[sinid]->v[i] - hier->obj[binid]->obj[mbhid]->v[i];
            }

            rmag = fb_mod(r_imbh_frame);
            for(i=0; i<3; i++) rhat[i] /= rmag;

            /*Set radial position*/
            star_r[get_global_idx(knew1)] = rmag*cmc_units.l; 

            /*Set velocities as well*/
            v_rmag = fb_dot(rhat,v_imbh_frame);
            for(i=0; i<3; i++) v_t[i] = v_imbh_frame[i] - rhat[i]*v_rmag;

            star[knew1].vr = v_rmag*cmc_units.v;
            star[knew1].vt = fb_mod(v_t)*cmc_units.v;

            /* SECOND STAR (inner binary star) */
            /* Extract the position/velocity wrt the IMBH (assumed cluster center) from fewbody*/
            for(i=0; i<3; i++){
                r_imbh_frame[i] = hier->obj[binid]->obj[1-mbhid]->x[i] - hier->obj[binid]->obj[mbhid]->x[i];
                v_imbh_frame[i] = hier->obj[binid]->obj[1-mbhid]->v[i] - hier->obj[binid]->obj[mbhid]->v[i];
            }

            rmag = fb_mod(r_imbh_frame);
            for(i=0; i<3; i++) rhat[i] /= rmag;

            /*Set radial position*/
            star_r[get_global_idx(knew2)] = rmag*cmc_units.l; 

            /*Set velocities as well*/
            v_rmag = fb_dot(rhat,v_imbh_frame);
            for(i=0; i<3; i++) v_t[i] = v_imbh_frame[i] - rhat[i]*v_rmag;

            star[knew2].vr = v_rmag*cmc_units.v;
            star[knew2].vt = fb_mod(v_t)*cmc_units.v;

            /*Set mass*/
            star_m[get_global_idx(knew1)] = hier->obj[sinid]->m * cmc_units.m/madhoc;
            star_m[get_global_idx(knew2)] = hier->obj[binid]->obj[1-mbhid]->m * cmc_units.m/madhoc;

            /*set potential*/
            star_phi[get_global_idx(knew1)] = potential(star_r[get_global_idx(knew1)]);
            star_phi[get_global_idx(knew2)] = potential(star_r[get_global_idx(knew2)]);

            /* Calculate new energies by recomputing E = PE + KE using new velocity */
            set_star_EJ(knew1);
            set_star_EJ(knew2);

            /* set rnew, vrnew, vtnew */
            set_star_news(knew1);
            set_star_news(knew2);
            
            /* I don't even know if this is necessary */
            set_star_olds(knew1);
            set_star_olds(knew2);

            /* mark stars as interacted */
            star[knew1].interacted = 1;
            star[knew2].interacted = 1;

            /* Copy stellar evolution parameters from binary member */
            /* reusing binid here*/
            if(hier->obj[0]->obj[sinid]->id[0] == binary[star[index].binind].id1)
                binid = 0;
            else
                binid = 1;
            cp_SEvars_to_newstar(index, binid, knew1);
            cp_SEvars_to_newstar(index, 1-binid, knew2);

            /* Destroy the original binary */
            destroy_obj(index);

            return 1;
        } else { /*only two objects; must have been a merger (either binary or TDE)*/

            if ((hier->obj[0]->id[0] == 0) && (hier->obj[0]->ncoll == 1))
                sinid = 1;
            else if ((hier->obj[1]->id[0] == 0) && (hier->obj[1]->ncoll == 1))
                mbhid = 1;
            
            if ((mbhid == 1) || (binid == 1)){ /*The MBH is one of the top-level objects; the binary must have merged*/

                    /* Create a new star for the binary merger*/
                    knew = create_star(index, 0);

                    /*Extract the position/velocity wrt the IMBH (assumed cluster center) from fewbody*/
                    for(i=0; i<3; i++){
                        r_imbh_frame[i] = hier->obj[sinid]->x[i] - hier->obj[mbhid]->x[i];
                        v_imbh_frame[i] = hier->obj[sinid]->v[i] - hier->obj[mbhid]->v[i];
                    }

                    rmag = fb_mod(r_imbh_frame);
                    for(i=0; i<3; i++) rhat[i] /= rmag;

                    /*Set radial position*/
                    star_r[get_global_idx(knew)] = rmag*cmc_units.l; 

                    /*Set velocities as well*/
                    v_rmag = fb_dot(rhat,v_imbh_frame);
                    for(i=0; i<3; i++) v_t[i] = v_imbh_frame[i] - rhat[i]*v_rmag;

                    star[knew].vr = v_rmag*cmc_units.v;
                    star[knew].vt = fb_mod(v_t)*cmc_units.v;

                    /*Set mass; this gets overwritten below*/
                    star_m[get_global_idx(knew)] = hier->obj[sinid]->m * cmc_units.m/madhoc;

                    /*Binary merger, so set internal energy*/
                    star[knew].Eint = hier->obj[sinid]->Eint*cmc_units.E;

                    /*set potential*/
                    star_phi[get_global_idx(knew)] = potential(star_r[get_global_idx(knew)]);

                    /* Calculate new energies by recomputing E = PE + KE using new velocity */
                    set_star_EJ(knew);

                    /* set rnew, vrnew, vtnew */
                    set_star_news(knew);
                    
                    /* I don't even know if this is necessary */
                    set_star_olds(knew);

                    /* mark stars as interacted */
                    star[knew].interacted = 1;
                    
                    /* Finally, actually merge the two stars using COSMIC */
                    cp_SEvars_to_star(index, 0, &tempstar1);
                    cp_m_to_star(index, 0, &tempstar1);
                    cp_SEvars_to_star(index, 1, &tempstar2);
                    cp_m_to_star(index, 1, &tempstar2);
                    merge_two_stars(&tempstar1, &tempstar2, &(star[knew]), vs, curr_st);
                    star[knew].vr += vs[3] * 1.0e5 / (units.l/units.t);
                    vt_add_kick(&(star[knew].vt),vs[1],vs[2], curr_st);

                    /* Destroy the original binary */
                    destroy_obj(index);

                    return 4; 
                } else { /*Must have been a single TDE*/

                    /* find the remainind unTDE'd star and the MBH*/
                    if(hier->obj[0]->ncoll > 1)
                        sinid = 1;
                    else 
                        mbhid = 1;

                    knew = create_star(index, 0);
                    star[knew].id = hier->obj[sinid]->id[0];

                    /*Extract the position/velocity wrt the IMBH (assumed cluster center) from fewbody*/
                    for(i=0; i<3; i++){
                        r_imbh_frame[i] = hier->obj[sinid]->x[i] - hier->obj[mbhid]->x[i];
                        v_imbh_frame[i] = hier->obj[sinid]->v[i] - hier->obj[mbhid]->v[i];
                    }

                    rmag = fb_mod(r_imbh_frame);
                    for(i=0; i<3; i++) rhat[i] /= rmag;

                    /*Set radial position*/
                    star_r[get_global_idx(knew)] = rmag*cmc_units.l; 

                    /*Set velocities as well*/
                    v_rmag = fb_dot(rhat,v_imbh_frame);
                    for(i=0; i<3; i++) v_t[i] = v_imbh_frame[i] - rhat[i]*v_rmag;

                    star[knew].vr = v_rmag*cmc_units.v;
                    star[knew].vt = fb_mod(v_t)*cmc_units.v;

                    /*Set mass*/
                    star_m[get_global_idx(knew)] = hier->obj[sinid]->m * cmc_units.m/madhoc;

                    /*set potential*/
                    star_phi[get_global_idx(knew)] = potential(star_r[get_global_idx(knew)]);

                    /* Calculate new energies by recomputing E = PE + KE using new velocity */
                    set_star_EJ(knew);

                    /* set rnew, vrnew, vtnew */
                    set_star_news(knew);
                    
                    /* I don't even know if this is necessary */
                    set_star_olds(knew);

                    /* mark stars as interacted */
                    star[knew].interacted = 1;

                    /* Copy stellar evolution parameters from binary member */
                    if(hier->obj[sinid]->id[0] == binary[star[index].binind].id1)
                        binid = 0;
                    else
                        binid = 1;
                    cp_SEvars_to_newstar(index, binid, knew);

                    /*Add mass of other star to the MBH*/
                    cenma.m_new = hier->obj[mbhid]->m * cmc_units.m/madhoc; 

                    /*Energy too, note here we're incrementing the internal energy (which was zero initially)*/
                    cenma.E_new += hier->obj[mbhid]->Eint * cmc_units.E; 
                    /*Same caveat on energy conservation as above*/

                    /* Destroy the original binary */
                    destroy_obj(index);

                    return 2; 
                }
            } 
    } else if (hier->nobj == 3){ /* Three unbound objects; can only be binary disruption*/

        /*Here we can just cycle through the objects*/

        int sinids[2];
        int j=0;

        for(i=0 ; i < 3 ; i++){ /*first find the MBH*/
            if (hier->obj[i]->id[0] == 0)
                mbhid = i;
            else 
                sinids[j++];
        }

        for(j=0 ; j < 2 ; i++){

            /* Create star and re-insert into cluster */
            knew = create_star(index, 0);
            star[knew].id = hier->obj[sinids[j]]->id[0];

            /* Extract the position/velocity wrt the IMBH (assumed cluster center) from fewbody*/
            for(i=0; i<3; i++){
                r_imbh_frame[i] = hier->obj[sinids[j]]->x[i] - hier->obj[mbhid]->x[i];
                v_imbh_frame[i] = hier->obj[sinids[j]]->v[i] - hier->obj[mbhid]->v[i];
            }

            rmag = fb_mod(r_imbh_frame);
            for(i=0; i<3; i++) rhat[i] /= rmag;

            /*Set radial position*/
            star_r[get_global_idx(knew)] = rmag*cmc_units.l; 

            /*Set velocities as well*/
            v_rmag = fb_dot(rhat,v_imbh_frame);
            for(i=0; i<3; i++) v_t[i] = v_imbh_frame[i] - rhat[i]*v_rmag;

            star[knew].vr = v_rmag*cmc_units.v;
            star[knew].vt = fb_mod(v_t)*cmc_units.v;

            /*Set mass*/
            star_m[get_global_idx(knew)] = hier->obj[sinids[j]]->m * cmc_units.m/madhoc;

            /*set potential*/
            star_phi[get_global_idx(knew)] = potential(star_r[get_global_idx(knew)]);

            /* Calculate new energies by recomputing E = PE + KE using new velocity */
            set_star_EJ(knew);

            /* set rnew, vrnew, vtnew */
            set_star_news(knew);
            
            /* I don't even know if this is necessary */
            set_star_olds(knew);

            /* mark stars as interacted */
            star[knew].interacted = 1;

            /* Copy stellar evolution parameters from binary member */
            /* reusing binid here*/
            if(hier->obj[sinids[j]]->id[0] == binary[star[index].binind].id1)
                binid = 0;
            else
                binid = 1;
            cp_SEvars_to_newstar(index, binid, knew);
        }

        /* Destroy the original binary */
        destroy_obj(index);

        return 1;
    }
}
