/*******************************************************************/
/*  cyl:    Generate photon field for given radiation type,        */ 
/*          optical properties of meduim, and OM response.         */
/*                                                                 */
/*  written: R. Bay  2/98                                          */
/*                                                                 */
/*  based on similar codes by:     Bay, Jacobsen, Liubarsky,       */
/*                                 Lowder, Moorehead, Woschnagg,   */
/*                                 and others                      */
/*******************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include "cfortran.h"
#include "hbook.h"     /* define hbook routines */

#define SPEED      0.29979    /* speed of light (m/ns) */
#define N_ICE      1.32       /* refractive index of ice */
#define PI         3.1415927
#define EPS        1.0e-4     /* small distance to avoid errors */
#define MU_LEN     520        /* muon attenutation length */
#define CER_COS    0.743      /* cosine of cerenkov cone */
#define PMT_RAD    0.03        /* photomultiplier radius */
#define DYN_EFF    0.9        /* efficiency of first dynode */

#define MAX_STEP   1000000

#define Z_START    -50.0

#define DROP    101

#define TRUE       1
#define FALSE      0
#define ERR        (-1)

#define NUM_COEF   11         /* number of coef. for acceptance poly */

#define NUM_WVL    31         /* number of bins in wavelength and */
#define NUM_ANG    37         /* angle dependent curves  */

#define ABSORBED   1
#define CROSSED    2          /* photon interaction types */
#define LAYER      3        
#define SCATTERED  4

#define NBINS_X     100
#define XMIN        0.
#define XMAX       4000.

#define HBKSIZE    500000
int pawc_[HBKSIZE];
int icycle, rec_len=1024, istat;   /* PAW stuff */

/*
#define ETAGSIZ    6
typedef char HTAGTAB[9];
HTAGTAB etags[ETAGSIZ] = {"time", "cyl", "x", "y", "z", "weight"}; 
*/

time_t *dummy_time_pointer;
#define random            ((long)time(dummy_time_pointer))
/* rnd1:  random number on [0,1] inclusive */
#define rnd1()          (float) rand() / (float) RAND_MAX
/* rnd2:  random number on (0,1) exclusive */
#define rnd2()          (rand() + 1.0)/((double)(RAND_MAX)+2.0)

typedef int boolean;

typedef struct {float x,y,z;} vect;

typedef struct {     /* photon information */
  vect pos,dir;
  float length;
  float t_emission;
  float lambda;
  float lam_scat;    /* geometric scattering length */
  float scat_dist;   /* distance left until next scatter */
  float abs_dist;    /* distance left until absorption */  
  int region;
  int lay_reg;
} photon;

typedef struct {     /* shell information */
  float rad;
  float band_area; 
} cyl;

typedef struct {     /* layer information */
  float rad_in;
  float rad_out; 
  float ratio_abs;
  float ratio_scat;
} layer;

/* wavelength bins for medium optical properties, OM response */
static float wvlens[NUM_WVL] = {
  300., 310., 320., 330., 340., 350., 360., 370., 380., 390.,
  400., 410., 420., 430., 440., 450., 460., 470., 480., 490.,
  500., 510., 520., 530., 540., 550., 560., 570., 580., 590.,
  600.};


/* absorption length in medium versus wavelength */
/*static float abslens[NUM_WVL] = {
  100., 100., 100., 100., 100., 100., 100., 100., 100., 100.,
  80., 80., 80., 80., 80., 76., 76., 76., 76., 76., 50., 50.,
  50., 50., 50., 30., 30., 30., 30., 30., 10.}; */

/* absorption length in medium versus wavelength */
static float abslens[NUM_WVL] = {
  100., 100., 100., 100., 100., 100., 100., 100., 100., 100.,
  80., 80., 80., 80., 80., 76., 76., 76., 76., 76., 50., 23.,
  23., 23., 23., 23., 23., 23., 23., 23., 23.};

/* scattering length in medium versus wavelength */
static float scatlens[NUM_WVL] = {
  5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 
  5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 
  5.};

/* OM quantum efficiency versus wavelength */
static float quantum_eff[NUM_WVL] = {
  0.01, 0.22, 0.265, 0.29, 0.295, 0.3, 0.3, 0.3, 0.295, 0.29,
  0.285, 0.275, 0.27, 0.26, 0.25, 0.24, 0.225, 0.215, 0.20, 
  0.18, 0.16, 0.14, 0.12, 0.095, 0.08, 0.065, 0.05, 0.035,
  0.03, 0.02, 0.015};

/* Pressure vessel transmission versus wavelength */
#define VESSEL 1
#ifdef VESSEL         /* benthos pressure vessel */
  static float trans_eff[NUM_WVL] = {
    9.92e-4, 0.1, 0.2, 0.31, 0.48, 0.63, 0.66, 0.71, 0.73, 0.74,
    0.74, 0.74, 0.75, 0.75, 0.75, 0.75, 0.75, 0.76, 0.76, 0.77,
    0.78, 0.79, 0.81, 0.82, 0.83, 0.83, 0.83, 0.82, 0.81, 0.81,
    0.81};
#else                 /* billings pressure vessel */
  static float trans_eff[NUM_WVL] = {
    9.92e-4, 2.97e-3, 8.88e-3, 2.62e-2, 7.44e-2, 0.193, 0.409,
    0.653, 0.815, 0.888, 0.916, 0.925, 0.928, 0.93, 0.93, 0.93,
    0.93, 0.93, 0.93, 0.93, 0.93, 0.93, 0.93, 0.93, 0.93, 0.93,
    0.93, 0.93, 0.93, 0.93, 0.93};
#endif

/* OM fresnel correction versus impact angle */
static float fres_coef[NUM_ANG] = {
  1.067, 1.067, 1.066, 1.065, 1.063, 1.062, 1.061, 1.06, 1.059,
  1.059, 1.059, 1.059, 1.06, 1.06, 1.061, 1.062, 1.063, 1.065,
  1.067, 1.069, 1.072, 1.075, 1.079, 1.084, 1.090, 1.098, 1.107, 
  1.119, 1.134, 1.152, 1.176, 1.208, 1.249, 1.301, 1.371, 1.462,
  0.0};

/* function prototypes */
void init_hbook(char *, float, float, int);
void term_hbook();
int get_int(char *);
float get_float(char *);
boolean set_up_geometry(float, int, cyl *, float);
boolean set_up_layer(layer *, float, float, float, float);
boolean create_photon(photon *, layer *, float, int, int, float, float, 
		      float, float);
boolean polint(float[], float[], int, float, float *, float *);
boolean linint(float[], float[], int, float, float *);
boolean propagate_photon(cyl *, photon *, layer *, int, float, float, 
			 float, float *, int, float);
boolean scatter_photon(photon *, layer *, float);
boolean dist_to_layer(layer *, photon *, float *, int *);
boolean cross_layer(layer*, photon *, float , int );
boolean distance_to_boundary(cyl *, photon *, float *, int *);
boolean distance_to_sphere(photon *, float *, float);
boolean distance_to_cylinder(photon *, float *, float);
boolean cross_boundary(photon *, cyl *, float, float, int, float, float *, int,
		       float);
float cross_section(double);
int find_next_int(cyl *, photon *, layer *, float *, int *, float *, int *);
float dot_product(vect, vect);
float dot_prod2d(vect, vect);
float rangauss();
float *vector(int, int);
void free_vector(float *, int, int);
void die(char *);


/******************************************************************/
/******************************************************************/
main()
{
  register int loop, loop2; /* main photon counter */
  int book;          /* book paw histos? (0=no/1=yes) */
  int n_cyl;         /* number of cylinders */
  int n_photon;      /* number of photons to generate */
  int emission;      /* isotropic(0), costheta(1), or cerenkov(2) */
  int spectrum;      /* fixed wavelength(0) or cerenkov(1) */
  float num_hits;
  float log_length;    /* length of dust_logger */
  float band;
  float scale;       /* scaling factor for number of photons */
  float lambda;      /* wavelength of fixed wavelength source */
  float av_cos;      /* average cosine of mie scatter */
  float lay_begin;  
  float lay_end;    
  float r_abs;    
  float r_scat;    
  float cyl_length;  /* largest cylinder radius */
  float cyl_rho_max; /* largest cylinder radius */
  float sigma;       /* sigma of gaussian time schmear */
  photon p;          /* photon information */
  layer lay[3];
  cyl *shell;        /* shell information */


  /****************************************************************/
  /* executables...                                               */

  srand(random);     /* initialize random # generator */

  /* user inputs */
  book = get_int("Book PAW histograms? (1=yes/0=no): ");
  n_photon = get_int("Enter number of photons: ");
  scale = get_float("Enter photon number scaling factor: ");
  log_length = get_float("Enter dust logger length: ");
  lay_begin = get_float("Enter layer start postition: ");
  lay_end = get_float("Enter layer end position: ");
  r_abs = get_float("Enter layer absorption ratio: ");
  r_scat = get_float("Enter layer scattering ratio: ");
  cyl_length = get_float("Enter cylinder length: ");
  cyl_rho_max = get_float("Enter maximum cylinder radius: ");
  n_cyl = get_int("Enter number of cylinders: ");
  band = get_float("Enter width of bands: ");
  sigma = get_float("Enter sigma: ");
  emission = 
    get_int("Enter emission type [0=isotropic, 1=cos(theta), 2=cerenkov]: ");
  spectrum = 
    get_int("Enter spectral type [0=fixed wavelength, 1=cerenkov]: ");
  av_cos = get_float("Enter average cosine of Mie scattering: ");
  lambda = get_float("Enter wavelength for fixed wavelength source: ");

  if (book) init_hbook("cyl.hbk", band, cyl_length, n_cyl);

  /* allocate memory based on user inputs */
  if(!(shell = (cyl*)malloc((n_cyl+1)*(sizeof (cyl)))))
    die("Allocation problem.");

  /* set up dust layer  */
  if(!(set_up_layer(lay, lay_begin, lay_end, r_abs, r_scat)))
    die("Geometry definition failed.");
  
  /* set up shell geometry  */
  if(!(set_up_geometry(cyl_rho_max, n_cyl, shell, band)))
    die("Geometry definition failed.");
  
  /* create and propagate n photons */
  
  for (loop=0; loop<n_photon; loop++) {
    printf("Created\n");
    if(!(create_photon(&p, lay, log_length, emission, spectrum, lambda, lay_begin, 
		       lay_end, cyl_length))) 
      die("Photon creation failed.");
    if(!(propagate_photon(shell, &p, lay, n_cyl, av_cos, sigma, scale, &num_hits, book,
			  cyl_length))) 
      die("Photon propagation failed.");
  }
  
  /* clean up and end */    
  if (book) term_hbook();

  free(shell);

}
/****************************************************************/
/****************************************************************/


/****************************************************************/
/*  set_up_geometry:   set up shell radii based on the user     */
/*                     defined max radius and number of shells  */
/*                     <nbins_rho>.  Find area of cpsi bins for */
/*                     use in hit weighting. Store in <shell[]> */
boolean set_up_geometry(float cyl_rho_max, int n_cyl, cyl *shell,
			float band)
{
  int i, j;
  float d_rad, d_cos;

  d_rad = cyl_rho_max / (float)n_cyl;

  for (i=0; i<n_cyl; i++){
    shell[i].rad = (float)(i+1)*d_rad;
    shell[i].band_area = band*2.0*PI*shell[i].rad;
  }

  /*  shell[n_cyl].rad = 2.0 * cyl_rho_max;   active volume outer bound */

  shell[n_cyl].rad = 2000.0 * cyl_rho_max;  /* active volume outer bound */

  return TRUE;
}
/****************************************************************/


/****************************************************************/
boolean set_up_layer(layer *lay, float lay_begin, float lay_end,
		     float r_abs, float r_scat)

{
  lay[0].rad_in = -9.0e5;
  lay[0].rad_out = lay_begin;
  lay[1].rad_in = lay_begin;
  lay[1].rad_out = lay_end;
  lay[2].rad_in = lay_end;
  lay[2].rad_out = 9.0e5;

  lay[0].ratio_abs = r_abs;
  lay[0].ratio_scat = r_scat;

  return TRUE;
}
/****************************************************************/


/****************************************************************/
/*  create_photon:  Initialize photon position and direction.   */
/*                  User specifies isotropic, costheta, or      */
/*                  cerenkov radiation along a track of user    */
/*                  specified length.  Allows for fixed wave-   */
/*                  length or cerenkov (~lambda^(-2)) spectrum. */
boolean create_photon(photon *p, layer *lay, float log_length, int emission, int spectrum, 
		      float lambda, float lay_begin, float lay_end, 
		      float cyl_length)
{
  float costheta, sintheta, phi;
  float abs, e_abs, scat, e_scat;
  float scale_factor;

  p->length = 0.0;
  p->region = 0;

  p->pos.x = 0.0;
  p->pos.y = 0.0;


  /* do
     p->pos.z = -2.0 * cyl_length + -MU_LEN * log(rnd2());
     while(p->pos.z > 2.0 * cyl_length); */
  
  p->pos.z = log_length;

  /* p->t_emission = (2.0*cyl_length + p->pos.z)/SPEED; */

  switch(emission) {
  case 0:  /* isotropic source */
    costheta = 2.0*rnd1()-1.0;
    sintheta = sqrt(1.0-costheta*costheta);
    break;
  case 1:  /* costheta source */
    sintheta = sqrt(rnd1());
    costheta = sqrt(1.0-sintheta*sintheta);
    break;
  case 2:  /* cerenkov photon */
    costheta = CER_COS;    /* down-going atmos. muons */
    sintheta = sqrt(1.0-costheta*costheta);
    break;
  case 3:
    costheta = 0.0;
    sintheta = 1.0;
    break;
  default:
    return FALSE;
    break;
  }

  phi = 2.0*PI*rnd1();

  p->dir.x = sintheta*cos(phi);
  p->dir.y = sintheta*sin(phi);
  p->dir.z = costheta;

  /* determine wavelength, absorption length, scattering length */
  switch(spectrum) {
  case 0:  /* fixed wavelength */
    linint(wvlens, abslens, NUM_WVL, lambda, &abs);

    scat=0.5;      /* bubbly!!!! */
    /* linint(wvlens, scatlens, NUM_WVL, lambda, &scat); */

    /* polint(wvlens, abslens, NUM_WVL, lambda, &abs, &e_abs); */
    /* polint(wvlens, scatlens, NUM_WVL, lambda, &scat, &e_scat); */
    break;
  case 1:  /* cerenkov spectrum */
    scale_factor = 1.0/wvlens[0] - 1.0/wvlens[NUM_WVL-1];
    lambda = 1.0/(1.0/wvlens[0] - rnd2()*scale_factor);
    linint(wvlens, abslens, NUM_WVL, lambda, &abs);
    linint(wvlens, scatlens, NUM_WVL, lambda, &scat);
    /* polint(wvlens, abslens, NUM_WVL, lambda, &abs, &e_abs); */
    /* polint(wvlens, scatlens, NUM_WVL, lambda, &scat, &e_scat); */
    break;
  default:
    return FALSE;
    break;
  }


  p->lambda = lambda;
  p->lam_scat = scat;
  
  if (lay_begin < p->pos.z)
    if (lay_end < p->pos.z)
      p->lay_reg = 2;
    else
      p->lay_reg = 1;
  else
    p->lay_reg =0;

  switch(p->lay_reg) {
  case 0: case 2:
    p->abs_dist = -abs * log(rnd2());
    /*p->abs_dist = -abs * 85;*/
    p->scat_dist = -scat * log(rnd2());
    break;
  case 1:
    p->abs_dist = -(abs/lay[0].ratio_abs) * log(rnd2());
    /*p->abs_dist = -(abs/lay[0].ratio_abs) * 85;*/
    p->scat_dist = -(scat/lay[0].ratio_scat) * log(rnd2());
    break;
  default:
    break;
  }

  return TRUE;
}
/****************************************************************/


/****************************************************************/
/*  propagate photon:  Photon is tracked from one interaction   */
/*                     to the next, be it a scatter, absorption,*/
/*                     or shell crossing(detection).            */

boolean propagate_photon(cyl *shell, photon *p, layer *lay, int n_cyl, 
			 float av_cos, float sigma, float scale, 
			 float *num_hits, int book, float cyl_length)
{
  float bound_dist, layer_dist;
  int i, int_type, new_region, new_lay_reg;

  for (i=0; i<MAX_STEP; i++) {
    switch (int_type = find_next_int(shell, p, lay, 
				     &bound_dist, &new_region,
				     &layer_dist, &new_lay_reg)) {
    case ABSORBED:
      return;
      break;
    case CROSSED:  
      if (new_region > n_cyl)  /* photon escaped active volume */
	return;                /* so kill it */
      else 
        if(!(cross_boundary(p, shell, bound_dist, sigma, new_region,
			    scale, num_hits, book, cyl_length)))
	   die("Boundary crossing failed.");
      break;   
    case LAYER:  
      if(!(cross_layer(lay, p, layer_dist, new_lay_reg)))
	die("Layer crossing failed.");
      break;
    case SCATTERED:
      if(!(scatter_photon(p, lay, av_cos))) die("Scattering error.");
      break;
    }
  }
  die("Photon propagation did not terminate!");
}
/****************************************************************/
	  
   
/****************************************************************/
/*  find_next_int:  for the current photon position and         */
/*                  direction, find next interaction type and   */
/*                  where it occurs.                            */
int find_next_int(cyl *shell, photon *p, layer *lay, 
		  float *bound_dist, int *new_region, 
		  float *layer_dist, int *new_lay_reg)
{
  int int_type;

  int_type = ABSORBED;
  if (!distance_to_boundary(shell, p, bound_dist, new_region))
    die ("Boundary distance measurement failed.");
  if (!dist_to_layer(lay, p, layer_dist, new_lay_reg))
    die ("Layer distance measurement failed.");

  if ((*bound_dist < p->abs_dist) && (*bound_dist < p->scat_dist)
      && (*bound_dist < *layer_dist))
    int_type = CROSSED;
  if ((*layer_dist < p->abs_dist) && (*layer_dist < p->scat_dist)
      && (*layer_dist < *bound_dist))
    int_type = LAYER;
  if ((p->scat_dist < p->abs_dist) && 
      (p->scat_dist < *bound_dist) &&
      (p->scat_dist < *layer_dist))
    int_type = SCATTERED;
  return int_type;
}
/****************************************************************/


/****************************************************************/
/*  cross_boundary:  track photon to boundary and update table. */

boolean cross_boundary(photon *p, cyl *shell, float bound_dist, 
		       float sigma, int new_region, float scale, 
		       float *num_hits, int book, float cyl_length)
{
  int i_weight;
  int shell_num;              /* shell which was hit */
  int itheta;                 /* bin for fresnel correction */
  int position, id;
  float dist, hit_time;
  float distance;
  /*  float evn[ETAGSIZ];         paw ntuple  */ 
  float qu_eff, tr_eff;       /* quantum, transmission efficiency */
  float impact_dot;       /* dot product pmt dir. and photon dir. */
  float impact_theta;     /* theta of dot product */
  float head_on_area = PI * PMT_RAD * PMT_RAD;  /* pmt cross section */
  float eff_area;         /* effective pmt area */
  float norm;            
  float surf_proj;        /* surface projection correction */
  float weight;           /* total weight of hit */
  vect pmt, norm_pos;     /* photon unit position, pmt dir. vectors */

  dist = bound_dist + EPS;
  p->pos.x += dist * (p->dir.x);
  p->pos.y += dist * (p->dir.y);  /* track photon to boundary */
  p->pos.z += dist * (p->dir.z);

  p->abs_dist -= dist;
  p->scat_dist -= dist;   /* update photon history */
  p->length += dist;

  if ((p->pos.z > -cyl_length/2.0) && (p->pos.z < cyl_length/2.0))
    {
      /* calculate hit time */
      /*      hit_time = p->t_emission + (p->length) * N_ICE/SPEED;
	      hit_time += sigma*rangauss(); */

      /* determine which shell was crossed */
      if ((p->region) > new_region)
	shell_num = new_region;
      else
	shell_num = p->region;
      assert(abs(p->region-new_region)==1);
      
      /* linearly extrapolate efficiencies from array for given lambda */
      /*      linint(wvlens, quantum_eff, NUM_WVL, p->lambda, &qu_eff);
	      linint(wvlens, trans_eff, NUM_WVL, p->lambda, &tr_eff);*/
      

      /*impact_dot =  p->dir.z; */


      if (p->dir.z>0) {

	eff_area = p->dir.z*head_on_area;
	if (shell[shell_num].band_area == 0.0) printf("%d Uh-oh!\n", shell_num);

	/* get a normalized position vector and find surface correction */
	norm = sqrt(dot_prod2d(p->pos, p->pos));
	norm_pos.x = p->pos.x / norm;
	norm_pos.y = p->pos.y / norm;
	surf_proj = fabs(dot_prod2d(p->dir, norm_pos));
	surf_proj = (surf_proj > 0.05) ? surf_proj : 0.05;

	weight = scale*eff_area/(shell[shell_num].band_area*surf_proj);
	/* weight *= exp(-((23*85)-p->abs_dist)/23);                 */

	printf("Detected\n");

	*num_hits+=weight;
	/*printf("%0.16f %0.16f\n", weight, *num_hits);*/
	
      }
      
      /*      impact_theta = acos(impact_dot);
	      itheta = (int)(impact_theta*0.2+0.5) + 1;  bin for fresnel */
      
      /* find total weight of hit */

      /* cross_section = big for +1 */

      /*      eff_area = cross_section(impact_dot)*head_on_area*
	      fres_coef[itheta]*qu_eff*tr_eff*DYN_EFF; */


      /* update ntuples */
      if (book) {
	/*    evn[0] = hit_time;
	      evn[1] = shell_num;
	      evn[2] = p->pos.x;
	      evn[3] = p->pos.y;
	      evn[4] = p->pos.z;
	      evn[5] = weight;
	      HFN(1,evn); */

       	position = (int)((p->pos.z+cyl_length/2.0)*2.0*PI*shell[shell_num].rad
			  /shell[shell_num].band_area);

	distance = sqrt(p->pos.x*p->pos.x+p->pos.y*p->pos.y+(p->pos.z-Z_START)*(p->pos.z-Z_START));

	id = 1000*(shell_num+1) + position;
	HF1(id, hit_time, weight); 
	id = 10000 + shell_num + 1;
	HF1(id, p->pos.z, distance*weight);
	id = 20000 + shell_num + 1;
	HF1(id, distance, distance*weight);
      }
    } 
  
  p->region = new_region;  /* break on through */
  return TRUE;
}
/****************************************************************/


/****************************************************************/
/*  cross_section:  determine OM unit cross section given       */
/*                  incidence dot product.                      */
float cross_section(double dotp)
{
  int k;
  float s;
  float a[NUM_COEF] = {0.44787, 0.45995, -0.0644593, 0.777876, 
		       0.325881, -1.07997, -0.17585, -0.178416, 
		       -0.469926, 0.524479, 0.444801};

  s = 0.0; 
  for(k=0; k<NUM_COEF; k++) { 
    s += a[k] * pow(dotp,(double)k); 
  }
  return s;

}
/****************************************************************/


/****************************************************************/
/* scatter_photon:  Advance photon to scatter point and change  */
/*                  direction.  Written by L. Bergstrom.        */
boolean scatter_photon(photon *p, layer *lay, float av_cos)
{
  float coscorr, sincorr;
  float phi, costheta, sintheta;
  float norm, r2;
  vect rperp;

  /* advance photon to scatter point */
  p->pos.x += (p->scat_dist) * (p->dir.x);
  p->pos.y += (p->scat_dist) * (p->dir.y);
  p->pos.z += (p->scat_dist) * (p->dir.z);

  printf("%f %f %f\n", p->pos.x, p->pos.y, p->pos.z); 

  /* update abs_dist and length */
  p->abs_dist -= (p->scat_dist);
  p->length += (p->scat_dist);

  /* find next scatter distance */
  if (p->lay_reg == 1)
    p->scat_dist =  -(p->lam_scat/lay[0].ratio_scat) * log(rnd2());
  else
    p->scat_dist =  -(p->lam_scat) * log(rnd2());
  
  /* find new direction */
  r2 = rnd2();
  coscorr=(1.0-av_cos+av_cos*r2)*r2*(1.0+av_cos)*(1.0+av_cos)/
    ((1.0-av_cos+2.0*av_cos*r2)*(1.0-av_cos+2.0*av_cos*r2));
  coscorr = 2.0*coscorr - 1.0;
  sincorr = sqrt(1.0 - coscorr*coscorr);
  phi = 2.0*PI*rnd1();
  costheta = 2.0*rnd1() - 1.0;
  sintheta = sqrt(1.0 - costheta*costheta);
  rperp.x=-sin(phi) * sintheta * (p->dir.z) + costheta * (p->dir.y);
  rperp.y=cos(phi) * sintheta * (p->dir.z) - costheta * (p->dir.x);
  rperp.z=sin(phi) * sintheta * (p->dir.x) - cos(phi) * 
    sintheta * (p->dir.y);
  norm = dot_product(rperp, rperp);
  norm = sqrt(norm);
  rperp.x /= norm;
  rperp.y /= norm;
  rperp.z /= norm;
  p->dir.x = coscorr * (p->dir.x) + sincorr * rperp.x;
  p->dir.y = coscorr * (p->dir.y) + sincorr * rperp.y;
  p->dir.z = coscorr * (p->dir.z) + sincorr * rperp.z;

  return TRUE;
}
/****************************************************************/


/****************************************************************/
boolean dist_to_layer(layer *lay, photon *p, float *dist, 
			  int *new_lay_reg)
{
  float zdist;

  if (p->dir.z == 0.0) 
    *dist = 9.0e5;
  else {
    if (p->dir.z < 0.0 ) {
      zdist = p->pos.z - lay[p->lay_reg].rad_in;
      *new_lay_reg = p->lay_reg - 1;
    }
    else {
      zdist = lay[p->lay_reg].rad_out - p->pos.z; 
      *new_lay_reg = p->lay_reg + 1;
    }
    *dist = zdist/fabs(p->dir.z);
  }

  /*  if (*dist < 0.0) printf("*** negative distance to boundary ***\n");*/
}
/****************************************************************/
  

/****************************************************************/
boolean cross_layer(layer *lay, photon *p, float lay_dist, 
		    int new_lay_reg)
{
  float dist;

  dist = lay_dist + EPS;
  p->pos.x += dist * (p->dir.x);
  p->pos.y += dist * (p->dir.y);  /* track photon to boundary */
  p->pos.z += dist * (p->dir.z);

  p->abs_dist -= dist;
  p->scat_dist -= dist;   /* update photon history */
  p->length += dist;

  if(new_lay_reg==1 && (p->lay_reg==0 || p->lay_reg==2)) {
    p->abs_dist /= lay[0].ratio_abs;
    p->scat_dist /= lay[0].ratio_scat;
  }
  if((new_lay_reg==0 || new_lay_reg==2) && p->lay_reg==1) {
    p->abs_dist *= lay[0].ratio_abs;
    p->scat_dist *= lay[0].ratio_scat;
  }

  p->lay_reg = new_lay_reg;
}


/****************************************************************/
/*  distance to boundary:  find the distance to nearest shell   */
/*                         for a given photon region, position, */
/*                         and direction.                       */
boolean distance_to_boundary(cyl *shell, photon *p,
			     float *dist, int *new_region)
{
  if ( p->region != 0 ) 
    /* check inner cylinder */
    if (distance_to_cylinder(p, dist, shell[(p->region)-1].rad)) {
      *new_region = (p->region) - 1;
      return (*dist > 0.0);
    }
  /* check outer outer */
  distance_to_cylinder(p, dist, shell[p->region].rad);
  *new_region = (p->region) + 1;
  return (*dist > 0.0); 
}
/****************************************************************/


/****************************************************************/
/*  distance_to_sphere:  find the distance to a sphere centered */
/*                       on the origin for given position and   */
/*                       direction.                             */
boolean distance_to_sphere(photon *p, float *dist, float rad)
{
  int flag;
  float dot, pos2, root;
  float dist1, dist2;

  flag = FALSE;
  dot = dot_product(p->pos, p->dir);
  pos2 = dot_product(p->pos, p->pos);
  root = dot*dot - pos2 + rad*rad;

  if (root > 0.0) 
    {
      dist1 = -dot + sqrt(root);
      dist2 = -dot - sqrt(root);
      if (dist1 > 0.0)  /* at least one positive solution exists */ 
	{
	  flag = TRUE;
	  /* check for second positive solution, choose minimum */
	  *dist = (dist2 > 0.0) ? ((dist1 < dist2) ? dist1 : dist2) :
	    dist1;
	}
      else if (dist2 > 0.0) /* only one solution */
	{
	  flag = TRUE;
	  *dist = dist2;
	}
    }
  return flag;
}
/****************************************************************/


/****************************************************************/
/*  distance_to_cylinder:  find the distance to a cylinder      */
/*                         centered on the origin for given     */
/*                         position and direction.              */
boolean distance_to_cylinder(photon *p, float *dist, float rad)
{
  int flag;
  float root, a, b, c, t, t1, t2;

  flag = FALSE;
  a = dot_prod2d(p->dir,p->dir);
  b = 2.0*dot_prod2d(p->dir, p->pos);
  c = - rad*rad + dot_prod2d(p->pos, p->pos);
  root = b*b - 4.0*a*c;
  
  if(root<0.0) 
    /*     the photon misses the cylinder */
    *dist = 9.0e5;
  else {
    t1=(-b + sqrt(root))/(2.0*a);
    t2=(-b - sqrt(root))/(2.0*a);
    if (t1 > 0.0)         /* at least one positive solution exists */ 
      {
	flag = TRUE;
        /* check for second positive solution, choose minimum */
	t = (t2 > 0.0) ? ((t1 < t2) ? t1 : t2) : t1;
      }
    else if (t2 > 0.0)         /* only one solution */
      {
	flag = TRUE;
	t = t2;
      }
    
    if (flag)
      *dist=t*sqrt(dot_product(p->dir, p->dir));
    else
      *dist=9.0e5;
    
    return flag;
  }
}
/****************************************************************/


/****************************************************************/
/* polint: from Numerical Recipes.  Given two arrays xa[], ya[],*/
/*          with y=f(x).  For a given x return a value of y and */
/*          estimated error dy.                                 */
boolean polint(float xa[], float ya[], int n, float x, float *y, 
	       float *dy)
{
  int i, m, ns=0;
  float den, dif, dift, ho, hp, w;
  float *c, *d, *vector();
  void free_vector();

  dif = fabs(x-xa[0]);
  c=vector(0,n-1);
  d=vector(0,n-1);
  for (i=0; i<n; i++) {
    if ( (dift = fabs(x-xa[i])) < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  *y = ya[ns--];
  for (m=0; m<(n-1); m++) {
    for (i=0; i<n-m; i++) {
      ho = xa[i]-x;
      hp = xa[i+m]-x;
      w = c[i+1] - d[i];
      if ( (den=ho-hp) == 0.0) return FALSE;
      den = w/den;
      d[i] = hp*den;
      c[i] = ho*den;
    }
    *y += (*dy = (2*ns < (n-1-m) ? c[ns+1] : d[ns--]));
  }
  free_vector(d, 0, n-1);
  free_vector(c, 0, n-1);
  return TRUE;
}
/****************************************************************/


/****************************************************************/
/*  vector:  from Numerical Recipes.  Function needed by        */
/*           polint().  Allocates a float vector with range     */
/*           [nl..nh].                                          */
float *vector(int nl, int nh)
{
  float *v;
  
  v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
  if (!v) die("Vector allocation failure.");
  return v-nl;
}
/****************************************************************/


/****************************************************************/
/*  free_vector:  from Numerical Recipes.  Function needed by   */
/*                polint().  Frees a float vector allocated by  */
/*                vector().                                     */
void free_vector(float *v, int nl, int nh)
{
  free((char*) (v+nl));
}
/****************************************************************/


/****************************************************************/
/*  linint:  linear interpolation given y[n] = f( x[n] ) and a  */
/*           given x.                                           */
boolean linint(float xa[], float ya[], int n, float x, float *y)
{
  float bin_size, index, fac;
  double bin;

  if ((x < xa[0]) || (x > xa[n-1])) die("Invalid wavelength.");
  bin_size = (xa[n-1] - xa[0])/(n-1);
  index = (x - xa[0])/bin_size;
  fac = modf(index, &bin);
  *y = ya[(int)bin] + fac * (ya[(int)bin+1] - ya[(int)bin]);
  return TRUE;
}
/****************************************************************/


/****************************************************************/
/*  dot_product:  dot two vectors                               */
float dot_product(vect v1, vect v2)
{
  return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}
/****************************************************************/

/****************************************************************/
/*  dot_product:  dot two vectors                               */
float dot_prod2d(vect v1, vect v2)
{
  return (v1.x * v2.x + v1.y * v2.y);
}
/****************************************************************/


/****************************************************************/
/*  rangauss:  generate unit gaussian distributed spread.       */
float rangauss()
{
  float v1, v2, r, fac;
  static float gset;
  static int iset;
  
  if (!iset) {
    for (r=37.0; r>=1.0;) {
      v1 = 2.0 * rnd1() - 1.0;
      v2 = 2.0 * rnd1() - 1.0;
      r = (v1*v1) + (v2*v2);
    }
    fac = sqrt(-2. * log(r)/r);
    gset = v1 * fac;
    iset = 1;
    return v2*fac;
  }
  else {
    iset=0;
    return gset;}
}
/****************************************************************/


/****************************************************************/
/*  get_int:  get integer user input                            */
int get_int(char *istring)
{
  int input;
  char line[256];

  printf(istring);
  gets(line);
  sscanf(line,"%d",&input);
  printf("%d\n", input);
  return input;
}
/****************************************************************/


/****************************************************************/
/*  get_float:  get float user input                            */
float get_float(char *fstring)
{
  float input;
  char line[256];

  printf(fstring);
  gets(line);
  sscanf(line,"%f",&input);
  printf("%f\n", input);
  return input;
}
/****************************************************************/


/****************************************************************/
/*  init_hbook:   initializes hbook.                            */
void init_hbook(char *hnam, float band, float cyl_length, int n_cyl)
{
  int i, j, n_bins;
  int id;
  float x_min, x_max;
  char* histnam;

  HLIMIT(HBKSIZE);
/*  HBOOKN(1,"Cyl Output",ETAGSIZ,"",5000,etags); */
  HROPEN(1,"Cyl output",hnam,"N",rec_len,istat);
  if(istat != 0) printf("Warning: HROPEN gives status %d.\n", istat);

  n_bins = (int)(cyl_length / band); 

  histnam = malloc(120 * sizeof(char));
  for(i=1;i<=n_cyl;i++) {
    for(j=1;j<=n_bins;j++) {
      id = 1000*i + j;
      sprintf(histnam,"Hit times for cylinder %d level %d", i, j);
      HBOOK1(id,histnam,NBINS_X,XMIN,XMAX,0.);
    }
  }
  
  x_min = -cyl_length/2.0;
  x_max = cyl_length/2.0;
  
  for(i=1;i<=n_cyl;i++) {
    id = 10000 + i;
    sprintf(histnam,"D * Fluence vs Z positions for cylinder %d", i);
    HBOOK1(id, histnam, n_bins, x_min, x_max, 0.);
  }
  for(i=1;i<=n_cyl;i++) {
    id = 20000 + i;
    sprintf(histnam,"D * Fluence vs D for cylinder %d", i);
    HBOOK1(id, histnam, n_bins, 0, cyl_length, 0.);
  }

}
/****************************************************************/


/****************************************************************/
/*  term_hbook:  terminates hbook.                              */
void term_hbook()
{
  /*  HRPUT(0,hnam,"N");*/
  HROUT(0,icycle," ");
  HREND("Cyl output");

}
/****************************************************************/


/****************************************************************/
/* die: error detected. bite it.                                */
void die(char* diestr)
{
  fprintf(stderr, diestr);
  fprintf(stderr, "\n");
  exit(ERR);
}
/****************************************************************/






