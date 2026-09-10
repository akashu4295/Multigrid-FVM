/*
  2D SIMPLE + (optional) MG skeleton with OpenACC acceleration.
  - Uses persistent device data (enter/exit data) for all SoA arrays.
  - Adds OpenACC pragmas to the main kernels.
  - Updates host only for Tecplot output.

*/

#include <stdio.h>
#include <stdlib.h>
#include <openacc.h>
#include <math.h>
#include <time.h>

/*********************************
nfx = no. of control volumes in x direction
nfy = no. of control volumes in y direction
ngrid = no. of multi-grid levels
ig_obs = grid level where obstacles are prescribed
ig_inlt = grid level where inlets are prescribed
ig_outlt = grid level where outlets are prescribed
nvcycle = no. of v cycle
xl = domain length in x-direction
yl = domain length in y-direction
**********************************/

int nfx, nfy, ngrid, nvcycle, nswpm, nswpp, isolve, nobs, ninlt, noutlt, nvcycle_converged;
int ig_obs, ig_inlt, ig_outlt;
int ifobs[1000], ilobs[1000], jfobs[1000], jlobs[1000], jfinlt[10], jlinlt[10],jfoutlt[10],jloutlt[10];
int ibc_xm, ibc_xp, ibc_ym, ibc_yp;
int ibeg_obs;
double error = 0.0, convergence = 0.00001;
double xl, yl, amuf, prsc, relxm, relxp, omega,dpdxm, dpdym;
double u_xm, v_xm, p_xm, sc_xm, u_xp, v_xp, p_xp, sc_xp;
double u_ym, v_ym, p_ym, sc_ym, u_yp, v_yp, p_yp, sc_yp;
double u_guess, v_guess, p_guess, sc_guess;
double uinlt[10], vinlt [10], scinlt[10], uoutlt[10], voutlt[10], scoutlt[10];
double convergence_data[3000];

/* Struct holding per-grid-level arrays */
struct SoA {
    int npx, npy;
    int *isobs, *ieobs, *jsobs, *jeobs, *tag;
    int *jsinlt, *jeinlt, *jsoutlt, *jeoutlt;
    double dx, dy,amu;
    double *areax, *areay;

    double *ae, *aw, *an, *as, *ap;
    double *u, *v, *p, *pp, *sc, *su, *rs;
    double *apu, *apv, *spu, *spv;
    double *resu, *resv, *ressc;
    double *c;
    double *cx, *cy;
};

/* Prototypes */

struct SoA* AllocateMemorySoA(void);
void readInput(void);
void FreeMemorySoA(struct SoA *mySoA);
void init(struct SoA *mySoA);
void rest_geom(struct SoA *mySoA);

void EnterDataSoA(struct SoA *mySoA);
void ExitDataSoA(struct SoA *mySoA);
void UpdateHostForOutput(const struct SoA *s);

void fluxes(struct SoA *igSoA);
void solve_u(struct SoA *igSoA);
void solve_v(struct SoA *igSoA);
void solve_pprime(struct SoA *igSoA, int icycle, int ig);
void solve_sc(struct SoA *igSoA);

void coefu(struct SoA *igSoA);
void umom(struct SoA *igSoA);
void coefv(struct SoA *igSoA);
void vmom(struct SoA *igSoA);
void coefp(struct SoA *igSoA);
double masserror(struct SoA *igSoA);
void update(struct SoA *igSoA);
void pbound(struct SoA *igSoA);
void coefsc(struct SoA *igSoA);
void scalar(struct SoA *igSoA);

void solver(struct SoA *igSoA, double *phi, int nx, int ny);
void solvergs(struct SoA *igSoA, double *phi, int nx, int ny);

void restv(struct SoA *igSoA, struct SoA *igSoA1);
void restru(struct SoA *igSoA, struct SoA *igSoA1);
void restrv(struct SoA *igSoA, struct SoA *igSoA1);
void restrsc(struct SoA *igSoA, struct SoA *igSoA1);

void prolu(struct SoA *igSoA, struct SoA *igSoA1);
void prolv(struct SoA *igSoA, struct SoA *igSoA1);
void prolp(struct SoA *igSoA, struct SoA *igSoA1);
void prolsc(struct SoA *igSoA, struct SoA *igSoA1);

void write_tecplot_2d_cc(const char *fname, const struct SoA *s, int iter,
                         double xl, double yl);
void write_vtk_2d_cc(const char *fname, const struct SoA *s, int iter);
void write_vtk_2d_cc_square(const char *fname, const struct SoA *s, int iter);
void write_convergence_data_file(const char *fname, double* convergence_data);
void read_obstacle_data(const char *filename);


static void *xmalloc(size_t nbytes)
{
    void *ptr = malloc(nbytes);
    if (!ptr) {
        fprintf(stderr, "ERROR: malloc failed for %zu bytes\n", nbytes);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

/* -------------------- OpenACC data lifetime helpers -------------------- */

void EnterDataSoA(struct SoA *mySoA)
{
    #pragma acc enter data copyin(mySoA[0:ngrid])
    //#pragma acc enter data create(convergence_data[0:3000])
    for (int ig = 0; ig < ngrid; ++ig) {
        int nx = mySoA[ig].npx, ny = mySoA[ig].npy;
        int nxy = nx * ny;

        #pragma acc enter data copyin( \
            mySoA[ig].u[0:nxy], mySoA[ig].v[0:nxy],  mySoA[ig].p[0:nxy], \
            mySoA[ig].pp[0:nxy], mySoA[ig].sc[0:nxy],mySoA[ig].su[0:nxy], mySoA[ig].rs[0:nxy], \
            mySoA[ig].resu[0:nxy], mySoA[ig].resv[0:nxy], mySoA[ig].ressc[0:nxy],\
            mySoA[ig].c[0:nxy], \
            mySoA[ig].aw[0:nxy], mySoA[ig].ae[0:nxy], mySoA[ig].as[0:nxy], mySoA[ig].an[0:nxy], \
            mySoA[ig].ap[0:nxy], \
            mySoA[ig].apu[0:nxy], mySoA[ig].apv[0:nxy], \
	   mySoA[ig].spu[0:nxy], mySoA[ig].spv[0:nxy], \
            mySoA[ig].cx[0:nxy], mySoA[ig].cy[0:nxy], \
	   mySoA[ig].areax[0:nxy], mySoA[ig].areay[0:nxy], \
            mySoA[ig].isobs[0:nobs], mySoA[ig].ieobs[0:nobs],\
            mySoA[ig].jsobs[0:nobs], mySoA[ig].jeobs[0:nobs],\
            mySoA[ig].jsinlt[0:ninlt], mySoA[ig].jeinlt[0:ninlt],\
            mySoA[ig].jsoutlt[0:noutlt], mySoA[ig].jeoutlt[0:noutlt], mySoA[ig].tag[0:nxy])
    }
}

void ExitDataSoA(struct SoA *mySoA)
{
    for (int ig = 0; ig < ngrid; ++ig) {
        int nx = mySoA[ig].npx, ny = mySoA[ig].npy;
        int nxy = nx * ny;

        #pragma acc exit data delete( \
            mySoA[ig].u[0:nxy], mySoA[ig].v[0:nxy],  mySoA[ig].p[0:nxy], \
            mySoA[ig].pp[0:nxy], mySoA[ig].sc[0:nxy],mySoA[ig].su[0:nxy], mySoA[ig].rs[0:nxy], \
            mySoA[ig].resu[0:nxy], mySoA[ig].resv[0:nxy], mySoA[ig].ressc[0:nxy],\
            mySoA[ig].c[0:nxy], \
            mySoA[ig].aw[0:nxy], mySoA[ig].ae[0:nxy], mySoA[ig].as[0:nxy], mySoA[ig].an[0:nxy], \
            mySoA[ig].ap[0:nxy], \
            mySoA[ig].apu[0:nxy], mySoA[ig].apv[0:nxy], \
	   mySoA[ig].spu[0:nxy], mySoA[ig].spv[0:nxy], \
            mySoA[ig].cx[0:nxy], mySoA[ig].cy[0:nxy], \
	   mySoA[ig].areax[0:nxy], mySoA[ig].areay[0:nxy], \
            mySoA[ig].isobs[0:nobs], mySoA[ig].ieobs[0:nobs],\
            mySoA[ig].jsobs[0:nobs], mySoA[ig].jeobs[0:nobs],\
            mySoA[ig].jsinlt[0:ninlt], mySoA[ig].jeinlt[0:ninlt],\
            mySoA[ig].jsoutlt[0:noutlt], mySoA[ig].jeoutlt[0:noutlt], mySoA[ig].tag[0:nxy])
    }
    //#pragma acc exit data delete(convergence_data[0:3000])
    #pragma acc exit data delete(mySoA[0:ngrid])
}

void UpdateHostForOutput(const struct SoA *s)
{
    int nxy = s->npx * s->npy;
    #pragma acc update self(s->u[0:nxy], s->v[0:nxy], s->p[0:nxy], s->pp[0:nxy], s->sc[0:nxy])
}

/* -------------------- main -------------------- */

int main(void)
{
    readInput();

    struct SoA *mySoA = AllocateMemorySoA();
 
    init(mySoA);
    rest_geom(mySoA);
    
    /* Create device copies once (persistent) */
    
    
    EnterDataSoA(mySoA);
    clock_t start = clock();
    int icyc = 0;
    for (icyc = 0; icyc < nvcycle; ++icyc) {

        for (int ig = 0; ig < ngrid; ++ig) {

            struct SoA *igSoA = &mySoA[ig];
 //         for (int iter = 0; iter < 3; ++iter) {
            fluxes(igSoA);
            solve_u(igSoA);
            solve_v(igSoA);
            solve_pprime(igSoA, icyc, ig);
            pbound(igSoA);
            solve_sc(igSoA);
//           }
 
            /* restriction */

            if (ig != ngrid - 1) {

                struct SoA *igSoA1 = &mySoA[ig + 1];
                restv(igSoA, igSoA1);
                fluxes(igSoA1);
                fluxes(igSoA);
                restru(igSoA, igSoA1);
                restrv(igSoA, igSoA1);
                restrsc(igSoA, igSoA1);
            }
        }

        /* Prolongation and Relaxation */
     
        for (int ig = ngrid - 1; ig > 0; ig--) {
			
            struct SoA *igSoA1 = &mySoA[ig];
            struct SoA *igSoA  = &mySoA[ig - 1];

            prolu(igSoA, igSoA1);
            prolv(igSoA, igSoA1);
            prolp(igSoA, igSoA1);
            prolsc(igSoA, igSoA1);

            if (ig != 1) {
 //           for (int iter = 0; iter < 3; ++iter) {
                fluxes(igSoA);
                solve_u(igSoA);
                solve_v(igSoA);
                solve_pprime(igSoA, icyc, ig);
                pbound(igSoA);
                solve_sc(igSoA);
            }
//            }                                                  
        }
        nvcycle_converged = icyc; 
        if (convergence_data[icyc]<convergence) {break;}       
    }                                                     
    
    //char fname[128];
    struct SoA *igSoA = &mySoA[0];
    //snprintf(fname, sizeof(fname), "tec_acc.dat");
    UpdateHostForOutput(igSoA);
    
    clock_t end = clock();
    double time_taken = (double)( end- start)/CLOCKS_PER_SEC; 
    printf("simple solution took %f seconds to execute \n", time_taken);
    
    //write_tecplot_2d_cc(fname, igSoA, 0, xl, yl);         
    
    char fname1[128];
    snprintf(fname1, sizeof(fname1), "vtk_acc.vtk");
    
    char fname2[128];
    snprintf(fname2, sizeof(fname2), "vtk_acc_square.vtk");

    if (convergence_data[icyc]<convergence){
	write_vtk_2d_cc(fname1, igSoA, 0); 
    write_vtk_2d_cc_square(fname2, igSoA, 0);
	write_convergence_data_file("convergence.csv", convergence_data);
    }
    
    ExitDataSoA(mySoA);
 
    
    FreeMemorySoA(mySoA);
    
    return 0;
}

/* -------------------- IO + memory -------------------- */

void readInput(void)
{
    FILE *fptr = fopen("InputData", "r");
    if (!fptr) {
        perror("ERROR opening inputData");
        exit(EXIT_FAILURE);
    }

    if (fscanf(fptr, "%d %d %d %d %d %d %d %d %d %d %d %lf %lf",
               &nfx, &nfy, &ngrid, &ig_obs, &ibeg_obs, &ig_inlt, &ig_outlt, &nvcycle, &nswpm, &nswpp, &isolve, &omega, &convergence) != 13) {
        fprintf(stderr, "ERROR: failed reading first line of inputData\n");
        exit(EXIT_FAILURE);
    }
       printf("nx_fine = %d, ny_fine = %d, ngrid = %d, nvcycle = %d, "
           "nswpm = %d, nswpp = %d, isolve = %d,omega = %lf, convergence = %e\n",
           nfx, nfy, ngrid, nvcycle, nswpm, nswpp, isolve, omega, convergence);


    if (fscanf(fptr, "%lf %lf %lf %lf %lf %lf",
               &xl, &yl, &amuf, &prsc, &relxm, &relxp) != 6) {
        fprintf(stderr, "ERROR: failed reading second line of inputData\n");
        exit(EXIT_FAILURE);
    }
    
    printf("xlength = %lf, ylength = %lf, viscosity = %lf, scalar Prandtl number = %lf, "
           "relax_momentum = %lf, relax_press = %lf\n",
           xl, yl, amuf, prsc, relxm, relxp);
           
    if (fscanf(fptr, "%d %d %d %d",
               &ibc_xm, &ibc_xp, &ibc_ym, &ibc_yp) != 4) {
        fprintf(stderr, "ERROR: failed reading third line of inputData\n");
        exit(EXIT_FAILURE);
    }        
    printf("ibc_xm = %d, ibc_xp = %d, ibc_ym = %d, ibc_yp = %d\n",
           ibc_xm , ibc_xp,  ibc_ym,  ibc_yp);
		   
    if (fscanf(fptr, "%lf %lf %lf %lf",
               &u_xm, &v_xm, &p_xm, &sc_xm) != 4) {
        fprintf(stderr, "ERROR: failed reading fourth line of inputData\n");
        exit(EXIT_FAILURE);
    }        
    printf("u_xm = %lf, v_xm = %lf, p_xm = %lf, sc_xm = %lf\n",
           u_xm , v_xm,  p_xm,  sc_xm);

    if (fscanf(fptr, "%lf %lf %lf %lf",
               &u_xp, &v_xp, &p_xp, &sc_xp) != 4) {
        fprintf(stderr, "ERROR: failed reading fifth line of inputData\n");
        exit(EXIT_FAILURE);
    }        
    printf("u_xp = %lf, v_xp = %lf, p_xp = %lf, sc_xp = %lf\n",
           u_xp , v_xp,  p_xp,  sc_xp);
		   
    if (fscanf(fptr, "%lf %lf %lf %lf",
               &u_ym, &v_ym, &p_ym, &sc_ym) != 4) {
        fprintf(stderr, "ERROR: failed reading sixth line of inputData\n");
        exit(EXIT_FAILURE);
    }        
    printf("u_ym = %lf, v_ym = %lf, p_ym = %lf, sc_ym = %lf\n",
           u_ym , v_ym,  p_ym,  sc_ym);

    if (fscanf(fptr, "%lf %lf %lf %lf",
               &u_yp, &v_yp, &p_yp, &sc_yp) != 4) {
        fprintf(stderr, "ERROR: failed reading seventh line of inputData\n");
        exit(EXIT_FAILURE);
    }        
    printf("u_yp = %lf, v_yp = %lf, p_yp = %lf, sc_yp = %lf\n",
           u_yp , v_yp,  p_yp,  sc_yp);
		   
    if (fscanf(fptr, "%lf %lf %lf %lf",
               &u_guess, &v_guess, &p_guess, &sc_guess) != 4) {
        fprintf(stderr, "ERROR: failed reading eighth line of inputData\n");
        exit(EXIT_FAILURE);
    }        
    printf("u_guess = %lf, v_guess = %lf, p_guess = %lf, sc_guess = %lf\n",
           u_guess , v_guess,  p_guess,  sc_guess);
		   
  
    read_obstacle_data("pattern.txt");
    
/*    if (fscanf(fptr, "%d", &nobs) != 1) {
        fprintf(stderr, "ERROR: failed reading obstacle line of inputData\n");
        fclose(fptr);
        exit(EXIT_FAILURE);
    }*/
  
     if (fscanf(fptr, "%d", &ninlt) != 1) {
        fprintf(stderr, "ERROR: failed reading inlet and outlet line of inputData\n");
        exit(EXIT_FAILURE);
    }
         printf("ninlt %d\n", ninlt);      
		   
     if (ninlt > 0) {
      for (int k = 0; k < ninlt; ++k) {    
        fscanf(fptr, " %d %d %lf %lf %lf",&jfinlt[k], &jlinlt[k],&uinlt[k], &vinlt[k], &scinlt[k]);
    }
  }
    if (ninlt > 0) {
         for (int k = 0; k < ninlt; ++k) {    
        printf(" jfinlt = %d, jlinlt = %d, uinlt = %lf, vinlt = %lf,tinlt = %lf\n",jfinlt[k],
                 jlinlt[k],uinlt[k],vinlt[k], scinlt[k]);
    }
  }
      if (fscanf(fptr, "%d", &noutlt) != 1) {
        fprintf(stderr, "ERROR: failed reading inlet and outlet line of inputData\n");
        exit(EXIT_FAILURE);
    }
         printf("noutlt %d \n", noutlt);      
     if (noutlt > 0) {
      for (int k = 0; k < noutlt; ++k) {    
        fscanf(fptr, " %d %d %lf %lf %lf",&jfoutlt[k], &jloutlt[k],&uoutlt[k], &voutlt[k], &scoutlt[k]);
    }
  }
        
  if (noutlt > 0) {
         for (int k = 0; k < noutlt; ++k) {    
        printf(" jfoutlt = %d, jloutlt = %d, uoutlt = %lf, voutlt = %lf,scoutlt = %lf\n",   
		         jfoutlt[k],jloutlt[k],uoutlt[k],voutlt[k], scoutlt[k]);
    }
  }           

    fclose(fptr);
}
struct SoA* AllocateMemorySoA(void)
{
    struct SoA *mySoA = (struct SoA*)xmalloc((size_t)ngrid * sizeof(struct SoA));

    for (int ig = 0; ig < ngrid; ++ig) {

        int fac = 1 << ig;

        int nx = nfx / fac + 2;
        int ny = nfy / fac + 2;

        int nxy = nx * ny;

        mySoA[ig].npx = nx;
        mySoA[ig].npy = ny;

        mySoA[ig].dx = xl / (nx - 2);
        mySoA[ig].dy = yl / (ny - 2);
        
        mySoA[ig].amu = amuf;

        mySoA[ig].u    = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].v    = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].p    = (double*)xmalloc((size_t)nxy * sizeof(double));

        mySoA[ig].pp   = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].sc   = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].su   = (double*)xmalloc((size_t)nxy * sizeof(double));

        mySoA[ig].rs   = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].resu = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].resv = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].ressc = (double*)xmalloc((size_t)nxy * sizeof(double));

        mySoA[ig].c    = (double*)xmalloc((size_t)nxy * sizeof(double));

        mySoA[ig].aw   = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].ae   = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].as   = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].an   = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].ap   = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].apu  = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].apv  = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].spu  = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].spv  = (double*)xmalloc((size_t)nxy * sizeof(double));

        mySoA[ig].cx   = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].cy   = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].tag  = (int*)xmalloc((size_t)nxy * sizeof(int));
		
        mySoA[ig].areax    = (double*)xmalloc((size_t)nxy * sizeof(double));
        mySoA[ig].areay    = (double*)xmalloc((size_t)nxy * sizeof(double));
        
        mySoA[ig].isobs= (int*)xmalloc((size_t)nobs * sizeof(int));
        mySoA[ig].ieobs= (int*)xmalloc((size_t)nobs * sizeof(int));
        mySoA[ig].jsobs= (int*)xmalloc((size_t)nobs * sizeof(int));
        mySoA[ig].jeobs= (int*)xmalloc((size_t)nobs * sizeof(int));

        mySoA[ig].jsinlt=(int*)xmalloc((size_t)ninlt * sizeof(int));
        mySoA[ig].jeinlt=(int*)xmalloc((size_t)ninlt * sizeof(int));

        mySoA[ig].jsoutlt=(int*)xmalloc((size_t)noutlt * sizeof(int));
        mySoA[ig].jeoutlt=(int*)xmalloc((size_t)noutlt * sizeof(int));

    }
     if (nobs >  0) {
 
      for (int k= 0; k<nobs;++k ){
        for (int ig = 0; ig < ig_obs; ++ig) {
	     int fac = 1 << (ig_obs - 1 - ig);	
        	     mySoA[ig].isobs[k] = (ifobs[k] + ibeg_obs - 1)* fac + 1 ;
              mySoA[ig].ieobs[k] = (ilobs[k] + ibeg_obs) * fac ;
              mySoA[ig].jsobs[k] = (jfobs[k]-1)* fac + 1 ;
              mySoA[ig].jeobs[k] = jlobs[k]* fac ;
              //printf(" isobs = %d, ieobs = %d\n",mySoA[ig].isobs[k], mySoA[ig].ieobs[k]);
	}
      }
    }
     
    if (ninlt >  0) {
 
      for (int k= 0; k<ninlt;++k ){
        for (int ig = 0; ig < ig_inlt; ++ig) {
	    int fac = 1 << (ig_inlt - 1 - ig);	
        mySoA[ig].jsinlt[k] = (jfinlt[k]-1)* fac + 1 ;
        mySoA[ig].jeinlt[k] = jlinlt[k]* fac ;
	    }
      }
    }

    if (noutlt >  0) {
 
      for (int k= 0; k<noutlt;++k ){
        for (int ig = 0; ig < ig_outlt; ++ig) {
	    int fac = 1 << (ig_outlt - 1 - ig);	
        mySoA[ig].jsoutlt[k] = (jfoutlt[k]-1)* fac + 1 ;
        mySoA[ig].jeoutlt[k] = jloutlt[k]* fac ;
		
	    }
      }
    }


    return mySoA;
}

void FreeMemorySoA(struct SoA *mySoA)
{
    if (!mySoA) return;

    for (int ig = 0; ig < ngrid; ++ig) {
        free(mySoA[ig].u);
        free(mySoA[ig].v);
        free(mySoA[ig].p);
        free(mySoA[ig].sc);

        free(mySoA[ig].pp);
        free(mySoA[ig].su);

        free(mySoA[ig].rs);
        free(mySoA[ig].resu);
        free(mySoA[ig].resv);
        free(mySoA[ig].ressc);
        free(mySoA[ig].c);
        free(mySoA[ig].aw);
        free(mySoA[ig].ae);
        free(mySoA[ig].as);
        free(mySoA[ig].an);
        free(mySoA[ig].ap);
        free(mySoA[ig].apu);
        free(mySoA[ig].apv);
        free(mySoA[ig].spu);
        free(mySoA[ig].spv);

        free(mySoA[ig].cx);
        free(mySoA[ig].cy);
        free(mySoA[ig].areax);
        free(mySoA[ig].areay);

    }

    free(mySoA);
}

void init(struct SoA *mySoA)
{

    for (int ig = 0; ig < ngrid; ++ig) {

            for (int j = 0; j < mySoA[ig].npy; ++j) {
                for (int i = 0; i < mySoA[ig].npx; ++i) {

                    int ij = j*mySoA[ig].npx + i;
                    mySoA[ig].u[ij]    = u_guess;
                    mySoA[ig].v[ij]    = v_guess;
                    mySoA[ig].sc[ij]   = sc_guess;
                    mySoA[ig].p[ij]    = p_guess;
                    mySoA[ig].pp[ij]   = 0.0;
                    mySoA[ig].su[ij]   = 0.0;
                    mySoA[ig].c[ij]    = 0.0;
                    mySoA[ig].ap[ij]   = 0.0;
                    mySoA[ig].resu[ij] = 0.0;
                    mySoA[ig].resv[ij] = 0.0;
                    mySoA[ig].ressc[ij]= 0.0;
                    mySoA[ig].rs[ij]   = 0.0;
                    mySoA[ig].tag[ij]  = 1;	 
                    mySoA[ig].apu[ij]  = 1.0e30;
                    mySoA[ig].apv[ij]  = 1.0e30;
	           mySoA[ig].spu[ij]  = 0.0;
                    mySoA[ig].spv[ij]  = 0.0;
		  mySoA[ig].areax[ij] = mySoA[ig].dy;
                    mySoA[ig].areay[ij] = mySoA[ig].dx;
                    mySoA[ig].cx[ij]   = 0.0;
                    mySoA[ig].cy[ij]   = 0.0;
                }
            }
    }

    /* boundary conditions on all boundaries */
    for (int ig = 0; ig < ngrid; ++ig) {
	
	    for (int j = 0; j < mySoA[ig].npy; ++j) {
                int ij = (j)*mySoA[ig].npx + 0;
                mySoA[ig].u[ij]  = u_xm;
	       mySoA[ig].v[ij]  = v_xm; 
	       mySoA[ig].p[ij]  = p_xm; 
	       mySoA[ig].sc[ij] = sc_xm;
        }
	    for (int j = 0; j < mySoA[ig].npy; ++j) {
                int ij = (j)*mySoA[ig].npx + mySoA[ig].npx-1;;
                mySoA[ig].u[ij-1]= u_xp;
	       mySoA[ig].v[ij]  = v_xp; 
	       mySoA[ig].p[ij]  = p_xp; 
	       mySoA[ig].sc[ij] = sc_xp;
        }
	    for (int i = 0; i < mySoA[ig].npx; ++i) {
                int ij = i;
                mySoA[ig].u[ij]  = u_ym;
	       mySoA[ig].v[ij]  = v_ym; 
                mySoA[ig].p[ij]  = p_ym; 
                mySoA[ig].sc[ij] = sc_ym;
        }
           for (int i = 0; i < mySoA[ig].npx; ++i) {
                int ij = (mySoA[ig].npy - 1)*mySoA[ig].npx + i;
                mySoA[ig].u[ij]  = u_yp;
	       mySoA[ig].v[ij]  = v_yp; 
	       mySoA[ig].p[ij]  = p_yp; 
	       mySoA[ig].sc[ij] = sc_yp;
        }
		
    }
    
	/*  zero velocities in obstacle region */
   
   if (nobs > 0) {
      for (int ig = 0; ig < ig_obs; ++ig) {
        for (int iobs = 0; iobs < nobs; ++iobs) {
	 for (int j = mySoA[ig].jsobs[iobs]; j <= mySoA[ig].jeobs[iobs]; ++j) {
           for (int i = mySoA[ig].isobs[iobs]; i <= mySoA[ig].ieobs[iobs]; ++i) {
          
             int ij = j*mySoA[ig].npx + i;
				
		     mySoA[ig].u[ij] = 0.0; 
		     mySoA[ig].v[ij] = 0.0;  
		     mySoA[ig].u[ij-1] = 0.0; 
		     mySoA[ig].v[ij-mySoA[ig].npx] = 0.0;  
	        //      mySoA[ig].areax[ij] = 0.0; 
		//     mySoA[ig].areay[ij] = 0.0;  
		//     mySoA[ig].areax[ij-1] = 0.0; 
		//     mySoA[ig].areay[ij-mySoA[ig].npx] = 0.0;  
		     mySoA[ig].spu[ij] = 1.0e30; 
		     mySoA[ig].spv[ij] = 1.0e30;  
		     mySoA[ig].spu[ij-1] = 1.0e30; 
		     mySoA[ig].spv[ij-mySoA[ig].npx] = 1.0e30;
		     mySoA[ig].tag[ij] = 0;				

           }
          }
        }
      }
 }
 
	/*  inlet velocities in inlet region */
   
   if (ninlt > 0) {
      for (int ig = 0; ig < ig_inlt; ++ig) {
        for (int iinlt = 0; iinlt < ninlt; ++iinlt) {
	  for (int j = mySoA[ig].jsinlt[iinlt]; j <= mySoA[ig].jeinlt[iinlt]; ++j) {
          
             int ij = j*mySoA[ig].npx;			
             mySoA[ig].u[ij] = uinlt[iinlt]; 
	    mySoA[ig].v[ij] = vinlt[iinlt]; 
             mySoA[ig].sc[ij]= scinlt[iinlt];		

           }
        }
      }
    }
	
	/*  outlet velocities in outlet region */

    if (noutlt > 0) {
      for (int ig = 0; ig < ig_outlt; ++ig) {
        for (int ioutlt = 0; ioutlt < noutlt; ++ioutlt) {
	 for (int j = mySoA[ig].jsoutlt[ioutlt]; j <= mySoA[ig].jeoutlt[ioutlt]; ++j) {
          
             int ij = j*mySoA[ig].npx + mySoA[ig].npx - 2;			
             mySoA[ig].u[ij] = uoutlt[ioutlt];
	    mySoA[ig].v[ij] = voutlt[ioutlt];
	    mySoA[ig].sc[ij+1]= scoutlt[ioutlt];  			

           }
        }
      }
    }	

}

void fluxes(struct SoA *igSoA){

    int nx = igSoA->npx, ny = igSoA->npy;
    int nxy = nx * ny;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->u[0:nxy], igSoA->v[0:nxy], \
        igSoA->cx[0:nxy], igSoA->cy[0:nxy], \
        igSoA->areax[0:nxy], igSoA->areay[0:nxy])
        for (int j = 0; j < ny-1; ++j) {
            for (int i = 0; i < nx-1; ++i) {
                int ij = j * nx + i;
                igSoA->cx[ij] = igSoA->u[ij] * igSoA->areax[ij];
                igSoA->cy[ij] = igSoA->v[ij] * igSoA->areay[ij];
            }
        }
}

void solve_u(struct SoA *igSoA) {

    coefu(igSoA);
    umom(igSoA);
    
    if (isolve == 1) solver  (igSoA, igSoA->u, igSoA->npx - 2, igSoA->npy - 1);
    if (isolve == 2) solvergs(igSoA, igSoA->u, igSoA->npx - 2, igSoA->npy - 1);
}

void solve_v(struct SoA *igSoA) {

    coefv(igSoA);
    vmom(igSoA);
    
    if (isolve == 1) solver  (igSoA, igSoA->v, igSoA->npx - 1, igSoA->npy - 2);
    if (isolve == 2) solvergs(igSoA, igSoA->v, igSoA->npx - 1, igSoA->npy - 2);

}

void solve_pprime(struct SoA *igSoA, int icycle, int ig) {

    double err_level = masserror(igSoA);
    coefp(igSoA);

    if (ig == 0) {
    	printf("icycle = %d error = %e\n", icycle, err_level);
    	convergence_data[icycle] = err_level;
    }
    if (isolve == 1) solver  (igSoA, igSoA->pp, igSoA->npx - 1, igSoA->npy - 1);
    if (isolve == 2) solvergs(igSoA, igSoA->pp, igSoA->npx - 1, igSoA->npy - 1);
    
    update(igSoA);
}
void solve_sc(struct SoA *igSoA) {
   
    coefsc(igSoA);
    scalar(igSoA);
    
    if (isolve == 1) solver  (igSoA, igSoA->sc, igSoA->npx - 1, igSoA->npy - 1);
    if (isolve == 2) solvergs(igSoA, igSoA->sc, igSoA->npx - 1, igSoA->npy - 1);
}

/* -------------------- momentum assembly -------------------- */

void umom(struct SoA* igSoA)
{
    double relxm1 = 1.0 - relxm;
    int nx = igSoA->npx, ny = igSoA->npy;
    int nxy = nx * ny;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->p[0:nxy], igSoA->su[0:nxy], igSoA->ap[0:nxy], \
        igSoA->apu[0:nxy], igSoA->u[0:nxy], igSoA->resu[0:nxy],\
        igSoA->areax[0:nxy], igSoA->areay[0:nxy]) copyin(dpdxm)
		
        for (int j = 1; j < ny-1; ++j) {
            for (int i = 1; i < nx-2; ++i) {
                int ij = j*nx + i;

                igSoA->su[ij]  = igSoA->areax[ij] * (igSoA->p[ij] - igSoA->p[ij+1]);
                igSoA->ap[ij]  = igSoA->ap[ij] / relxm;
                igSoA->su[ij]  = igSoA->su[ij]
                                + relxm1 * igSoA->ap[ij] * igSoA->u[ij]
                                + igSoA->resu[ij] + dpdxm*igSoA->areax[ij]*igSoA->areay[ij];
								
	       igSoA->ap [ij] = igSoA->ap[ij] + igSoA->spu[ij];
	       igSoA->apu[ij] = igSoA->ap[ij];
	       
            }
        }
}

void vmom(struct SoA *igSoA)
{
    double relxm1 = 1.0 - relxm;
    int nx = igSoA->npx, ny = igSoA->npy;
    int nxy= nx * ny;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->p[0:nxy], igSoA->su[0:nxy], igSoA->ap[0:nxy], \
        igSoA->apv[0:nxy], igSoA->v[0:nxy], igSoA->resv[0:nxy],\
        igSoA->areax[0:nxy], igSoA->areay[0:nxy]) copyin(dpdym)

        for (int j = 1; j < ny-2; ++j) {
            for (int i = 1; i < nx-1; ++i) {
                int ij = j*nx + i;

                igSoA->su[ij]  = igSoA->areay[ij] * (igSoA->p[ij] - igSoA->p[ij + nx]);
                igSoA->ap[ij]  = igSoA->ap[ij] / relxm;
                igSoA->su[ij]  = igSoA->su[ij]
                                + relxm1 * igSoA->ap[ij] * igSoA->v[ij]
                                + igSoA->resv[ij] + dpdym*igSoA->areax[ij]*igSoA->areay[ij];	
                igSoA->ap[ij] = igSoA->ap[ij] + igSoA->spv[ij];
	       igSoA->apv[ij] = igSoA->ap[ij];
            }
        }
}

void scalar(struct SoA *igSoA)
{
    int nx = igSoA->npx, ny = igSoA->npy;
    int nxy = nx * ny;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->p[0:nxy], igSoA->su[0:nxy], igSoA->ap[0:nxy], \
        igSoA->apu[0:nxy], igSoA->sc[0:nxy], igSoA->ressc[0:nxy])
		
        for (int j = 1; j < ny-1; ++j) {
            for (int i = 1; i < nx-1; ++i) {
                int ij = j*nx + i;
                igSoA->su[ij] = igSoA->ressc[ij];
								
            }
        }
}


void coefu(struct SoA *igSoA)
{
    int nx = igSoA->npx, ny = igSoA->npy;
    int nxy = nx * ny;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->cx[0:nxy], igSoA->cy[0:nxy], \
        igSoA->aw[0:nxy], igSoA->ae[0:nxy], igSoA->as[0:nxy], igSoA->an[0:nxy], \
        igSoA->ap[0:nxy], igSoA->areax[0:nxy], igSoA->areay[0:nxy])
		
        for (int j = 1; j < ny-1; ++j) {
            for (int i = 1; i < nx-2; ++i) {
                int ij  = j*nx + i;
                int ij1 = ij - nx;

                double cxm = 0.5 * (igSoA->cx[ij]  + igSoA->cx[ij-1]);
                double cxp = 0.5 * (igSoA->cx[ij]  + igSoA->cx[ij+1]);
                double cym = 0.5 * (igSoA->cy[ij1] + igSoA->cy[ij1+1]);
                double cyp = 0.5 * (igSoA->cy[ij]  + igSoA->cy[ij+1]);

                double dxm = igSoA->amu * (igSoA->areax[ij] / igSoA->dx);
                double dym = igSoA->amu * (igSoA->areay[ij] / igSoA->dy);
                
                igSoA->aw[ij] = fmax(fabs(cxm), dxm) + cxm;
                igSoA->ae[ij] = fmax(fabs(cxp), dxm) - cxp;
                igSoA->as[ij] = fmax(fabs(cym), dym) + cym;
                igSoA->an[ij] = fmax(fabs(cyp), dym) - cyp; 

                igSoA->ap[ij] = igSoA->aw[ij] + igSoA->ae[ij]
                              + igSoA->as[ij] + igSoA->an[ij];
            }
        }

    /* y-boundary modification */
    #pragma acc parallel loop collapse(1) present(igSoA[0:1], igSoA->an[0:nxy], igSoA->as[0:nxy], igSoA->ap[0:nxy])
        for (int i = 1; i < nx-2; ++i) {
            int ij  = (ny-2)*nx + i;
            int ij1 = (1)*nx    + i;

            double add = igSoA->amu * (igSoA->areay[ij] / igSoA->dy);
            
            igSoA->an[ij]  += add;  igSoA->ap[ij]  += add;
            igSoA->as[ij1] += add;  igSoA->ap[ij1] += add;
        }

}
void coefsc(struct SoA *igSoA)
{
    int nx = igSoA->npx, ny = igSoA->npy;
    int nxy = nx * ny;
    double amusc = amuf/prsc;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->cx[0:nxy], igSoA->cy[0:nxy], \
        igSoA->aw[0:nxy], igSoA->ae[0:nxy], igSoA->as[0:nxy], igSoA->an[0:nxy], \
        igSoA->ap[0:nxy], igSoA->areax[0:nxy], igSoA->areay[0:nxy])
		
        for (int j = 1; j < ny-1; ++j) {
            for (int i = 1; i < nx-1; ++i) {
                int ij  = j*nx + i;
                int ij1 = ij - nx;

                double cxm = igSoA->cx[ij-1];
                double cxp = igSoA->cx[ij];
                double cym = igSoA->cy[ij1];
                double cyp = igSoA->cy[ij];

                double dxm = amusc * (igSoA->areax[ij-1] / igSoA->dx);
                double dym = amusc * (igSoA->areay[ij1] / igSoA->dy);
                double dxp = amusc * (igSoA->areax[ij] / igSoA->dx);
                double dyp = amusc * (igSoA->areay[ij] / igSoA->dy);

                
                igSoA->aw[ij] = fmax(fabs(cxm), dxm) + cxm;
                igSoA->ae[ij] = fmax(fabs(cxp), dxp) - cxp;
                igSoA->as[ij] = fmax(fabs(cym), dym) + cym;
                igSoA->an[ij] = fmax(fabs(cyp), dyp) - cyp; 

                igSoA->ap[ij] = igSoA->aw[ij] + igSoA->ae[ij]
                              + igSoA->as[ij] + igSoA->an[ij];
            }
        }
        
    // x-boundary mods: i=1 and i=nx-2; zero flux boundary
    #pragma acc parallel loop collapse(1) present(igSoA[0:1], igSoA->ae[0:nxy], igSoA->aw[0:nxy], igSoA->ap[0:nxy])
 
        for (int j = 1; j < ny-1; ++j) {
            int ij  = j*nx + (nx-2);
            int ij1 = j*nx + 1;
            int flag = 0;
            for (int iinlt = 0; iinlt < ninlt; ++iinlt) {
	     if ((j >= igSoA->jsinlt[iinlt]) && (j <= igSoA->jeinlt[iinlt])) {
		flag = 1;    
              }
            }
            if (flag==0){
            	igSoA->ap[ij1] = igSoA->ap[ij1] - igSoA->aw[ij1];
            	igSoA->aw[ij1]  = 0.0;
            }
            flag = 0;
            for (int ioutlt = 0; ioutlt < noutlt; ++ioutlt) {
   	     if ((j >= igSoA->jsoutlt[ioutlt]) && (j <= igSoA->jeoutlt[ioutlt])) { 
		flag = 1;
	     }
	   }
	   if (flag==0){
	         igSoA->ap[ij] = igSoA->ap[ij] - igSoA->ae[ij];
            	igSoA->ae[ij]  = 0.0;
            }
          }  
        
    // y-boundary modification : zero flux boundary
    #pragma acc parallel loop collapse(1) present(igSoA[0:1], igSoA->an[0:nxy], igSoA->as[0:nxy], igSoA->ap[0:nxy])
        for (int i = 1; i < nx-1; ++i) {
            int ij  = (ny-2)*nx + i; 
            int ij1 = (1)*nx    + i;
 
            igSoA->ap[ij] = igSoA->ap[ij] - igSoA->an[ij];
            igSoA->ap[ij1] = igSoA->ap[ij1] - igSoA->as[ij1];
            igSoA->an[ij]  = 0.0;
            igSoA->as[ij1] = 0.0;           
            
        }

        if (nobs > 0) {
        for (int iobs = 0; iobs<nobs; ++iobs){
          for(int i = igSoA->isobs[iobs];i<=igSoA->ieobs[iobs]; ++i){
          for (int j = igSoA->jsobs[iobs];j<=igSoA->jeobs[iobs]; ++j){
            int ij  = j*nx + i;
           igSoA->ae[ij-1]  = 0.0;
           igSoA->aw[ij+1]  = 0.0;
           igSoA->an[ij-nx] = 0.0;
           igSoA->as[ij+nx] = 0.0;
          }
          } 
        }
        }  
}     

void coefv(struct SoA *igSoA)
{
    int nx = igSoA->npx, ny = igSoA->npy;
    int nxy = nx * ny;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->cx[0:nxy], igSoA->cy[0:nxy], \
        igSoA->aw[0:nxy], igSoA->ae[0:nxy], igSoA->as[0:nxy], igSoA->an[0:nxy], \
        igSoA->ap[0:nxy], igSoA->areax[0:nxy], igSoA->areay[0:nxy])
		
        for (int j = 1; j < ny-2; ++j) {
            for (int i = 1; i < nx-1; ++i) {
                int ij  = j*nx + i;
                int ij1 = ij + nx;
                int ij2 = ij - nx;

                double cxm = 0.5 * (igSoA->cx[ij-1] + igSoA->cx[ij1-1]);
                double cxp = 0.5 * (igSoA->cx[ij]   + igSoA->cx[ij1]);
                double cym = 0.5 * (igSoA->cy[ij]   + igSoA->cy[ij2]);
                double cyp = 0.5 * (igSoA->cy[ij]   + igSoA->cy[ij1]);

                double dxm = igSoA->amu * (igSoA->areax[ij] / igSoA->dx);
                double dym = igSoA->amu * (igSoA->areay[ij] / igSoA->dy);
                double dxp = igSoA->amu * (igSoA->areax[ij] / igSoA->dx);
                double dyp = igSoA->amu * (igSoA->areay[ij] / igSoA->dy);

                
                igSoA->aw[ij] = fmax(fabs(cxm), dxm) + cxm;
                igSoA->ae[ij] = fmax(fabs(cxp), dxp) - cxp;
                igSoA->as[ij] = fmax(fabs(cym), dym) + cym;
                igSoA->an[ij] = fmax(fabs(cyp), dyp) - cyp; 
                
                igSoA->ap[ij] = igSoA->aw[ij] + igSoA->ae[ij]
                              + igSoA->as[ij] + igSoA->an[ij];
            }
        }

    /* boundary mods: i=1 and i=nx-2 */
    #pragma acc parallel loop collapse(1) present(igSoA[0:1], igSoA->ae[0:nxy], igSoA->aw[0:nxy], igSoA->ap[0:nxy], igSoA->areay[0:nxy])
 
        for (int j = 1; j < ny-2; ++j) {
            int ij  = j*nx + (nx-2);
            int ij1 = j*nx + 1;

            double add = igSoA->amu * (igSoA->areay[ij] / igSoA->dy);
           
            igSoA->ae[ij]  += add;  igSoA->ap[ij]  += add;
            igSoA->aw[ij1] += add;  igSoA->ap[ij1] += add;
            
        }
}

void coefp(struct SoA *igSoA)
{
    int nx = igSoA->npx, ny = igSoA->npy;
    int nxy = nx * ny;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->apu[0:nxy], igSoA->apv[0:nxy], \
        igSoA->aw[0:nxy], igSoA->ae[0:nxy], igSoA->as[0:nxy], igSoA->an[0:nxy], igSoA->ap[0:nxy], \
        igSoA->areax[0:nxy], igSoA->areay[0:nxy], igSoA->tag[0:nxy])
  
        for (int j = 1; j < ny-1; ++j) {
            for (int i = 1; i < nx-1; ++i) {
                int ij  = j*nx + i;
                int ij1 = ij - nx;

                igSoA->aw[ij] = igSoA->areax[ij-1] * igSoA->areax[ij-1] / igSoA->apu[ij - 1];
                igSoA->ae[ij] = igSoA->areax[ij] * igSoA->areax[ij] / igSoA->apu[ij];

                igSoA->as[ij] = igSoA->areay[ij1] * igSoA->areay[ij1] / igSoA->apv[ij1];
                igSoA->an[ij] = igSoA->areay[ij] * igSoA->areay[ij] / igSoA->apv[ij];

                igSoA->ap[ij] = igSoA->aw[ij] + igSoA->ae[ij]
                              + igSoA->as[ij] + igSoA->an[ij] + (1 - igSoA->tag[ij])*1.0e30;
               
            }
        }		
         double diag = 1.0e30;
	if (ibc_xm == 2) {
		
       #pragma acc parallel loop collapse(1) present(igSoA[0:1], \
        igSoA->ap[0:nxy], igSoA->su[0:nxy])
  
        for (int j = 1; j < ny-1; ++j) {
                int ij  = j*nx + 1;
                igSoA->ap[ij] = igSoA->ap[ij] + diag;
	       igSoA->su[ij] = 0.0;				
            }
        }			
	if (ibc_xp == 2) {
		
	#pragma acc parallel loop collapse(1) present(igSoA[0:1], \
        igSoA->ap[0:nxy], igSoA->su[0:nxy])
  
        for (int j = 1; j < ny-1; ++j) {
                int ij  = j*nx + nx-2;
                igSoA->ap[ij] = igSoA->ap[ij] + diag;
	       igSoA->su[ij] = 0.0;
				
            }
        }	
	    if (ibc_ym == 2) {
		
	    #pragma acc parallel loop collapse(1) present(igSoA[0:1], \
        igSoA->ap[0:nxy], igSoA->su[0:nxy])
  
        for (int i = 1; i < nx-1; ++i) {
                int ij  = 1*nx + i;
                igSoA->ap[ij] = igSoA->ap[ij] + diag;
	       igSoA->su[ij] = 0.0;				
            }
        }	
		
	    if (ibc_yp == 2) {
		
	    #pragma acc parallel loop collapse(1) present(igSoA[0:1], \
        igSoA->ap[0:nxy], igSoA->su[0:nxy])
  
        for (int i = 1; i < nx-1; ++i) {
                int ij  = (ny-2)*nx + i;
                igSoA->ap[ij] = igSoA->ap[ij] + diag;
	       igSoA->su[ij] = 0.0;				
            }
        }	
		
}

/* -------------------- mass error + update -------------------- */

double masserror(struct SoA *igSoA)
{
    int nx = igSoA->npx, ny = igSoA->npy;
    int nxy = nx * ny;

    double err = 0.0;

    #pragma acc parallel loop collapse(2) reduction(+:err) present(igSoA[0:1], \
        igSoA->u[0:nxy], igSoA->v[0:nxy], \
        igSoA->su[0:nxy], igSoA->pp[0:nxy], igSoA->areax[0:nxy], igSoA->areay[0:nxy])

        for (int j = 1; j < ny-1; ++j) {
            for (int i = 1; i < nx-1; ++i) {
                int ij = j*nx + i;
                int ij1 = ij - nx;

                double div = igSoA->areax[ij] * igSoA->u[ij] - igSoA->areax[ij]* igSoA->u[ij-1] +
                             igSoA->areay[ij] * igSoA->v[ij] - igSoA->areay[ij] * igSoA->v[ij1];

                igSoA->su[ij] = -div;
                igSoA->pp[ij] =  0.0;

                err += fabs(div);
 
 //             printf("i =  %d, j = %d, div = %lf \n", i,j,div);  
            }
        }

    return err;
}

void update(struct SoA *igSoA)
{
    int nx = igSoA->npx, ny = igSoA->npy;
    int nxy = nx * ny;

    /* your code sets ppref = 0 anyway */
    double ppref = 0.0;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->u[0:nxy], igSoA->v[0:nxy], igSoA->p[0:nxy], igSoA->pp[0:nxy], \
        igSoA->apu[0:nxy], igSoA->apv[0:nxy], igSoA->areax[0:nxy],igSoA->areay[0:nxy])

        for (int j = 1; j < ny-1; ++j) {
            for (int i = 1; i < nx-1; ++i) {

                int ij  = j*nx + i;
                int ij1 = ij + 1;
                int ij2 = ij + nx;

                igSoA->u[ij] += igSoA->areax[ij] * (igSoA->pp[ij] - igSoA->pp[ij1]) / igSoA->apu[ij];
                igSoA->v[ij] += igSoA->areay[ij] * (igSoA->pp[ij] - igSoA->pp[ij2]) / igSoA->apv[ij];
                igSoA->p[ij] += relxp * (igSoA->pp[ij] - ppref);
            }
        }
        
        
       if (nobs > 0) {
       for (int iobs = 0; iobs < nobs; ++iobs) {
	 int j = igSoA->jsobs[iobs]; 
           for (int i = igSoA->isobs[iobs]; i <= igSoA->ieobs[iobs]; ++i) {
             int ij = j*igSoA->npx + i; 
             igSoA->p[ij] = igSoA->p[ij-nx];
           }
            j = igSoA->jeobs[iobs]; 
           for (int i = igSoA->isobs[iobs]; i <= igSoA->ieobs[iobs]; ++i) {
             int ij = j*igSoA->npx + i; 
             igSoA->p[ij] = igSoA->p[ij+nx];
           }
           
           int i = igSoA->isobs[iobs]; 
           for (int j = igSoA->jsobs[iobs]; j <= igSoA->jeobs[iobs]; ++j) {
             int ij = j*igSoA->npx + i; 
             igSoA->p[ij] = igSoA->p[ij-1];
           }
               i = igSoA->ieobs[iobs]; 
           for (int j = igSoA->jsobs[iobs]; j <= igSoA->jeobs[iobs]; ++j) {
             int ij = j*igSoA->npx + i; 
             igSoA->p[ij] = igSoA->p[ij+1];
           }
       }
       }
}
void pbound(struct SoA *igSoA)
{
    int nx = igSoA->npx, ny = igSoA->npy;
    int nxy = nx * ny;

   
 if (ibc_xm == 3) {
  #pragma acc parallel loop present(igSoA[0:1], \
        igSoA->u[0:nxy], igSoA->v[0:nxy])
        for (int j = 1; j < ny-1; ++j) { 
             int ij  = j*nx ;
	    int ij1 = ij + nx - 3;
	    igSoA->u[ij] = igSoA->u[ij1];
             igSoA->u[ij1+1] = igSoA->u[ij+1];
	    igSoA->v[ij] = igSoA->v[ij1];
	    igSoA->v[ij1+1] = igSoA->v[ij+1];      
            }
         }
         
  if (ibc_ym == 3) {
   #pragma acc parallel loop present(igSoA[0:1], \
        igSoA->u[0:nxy], igSoA->v[0:nxy] )
        for (int i = 1; i < nx-1; ++i) { 
                int ij  = i ;
                int ij1 = ij + ny-3;
	     igSoA->u[ij] = igSoA->u[ij1];
              igSoA->u[ij1+nx] = igSoA->u[ij+nx];
              igSoA->v[ij] = igSoA->v[ij1];
              igSoA->v[ij1+nx] = igSoA->v[ij+nx];
            }
         }
}
/* -------------------------jacobi-------------------------- */

void solver(struct SoA *igSoA, double *phi, int nx, int ny)
{
    int npx  = igSoA->npx;
    int npy  = igSoA->npy;
    int npxy = npx * npy;
    int nxy  = igSoA->npx * igSoA->npy;

    int nsw = (phi == igSoA->pp) ? nswpp : nswpm;

    /* Weighted Jacobi helps a lot for PP */
    double wj =1.0; (phi == igSoA->pp) ? 0.67 : 1.0;

    /* allocate workspace ONCE per call (better: store per grid and reuse) */
    double *phi_new = (double*)malloc((size_t)nxy * sizeof(double));
    if (!phi_new) { fprintf(stderr,"malloc phi_new failed\n"); exit(1); }

    #pragma acc data create(phi_new[0:nxy]) present(igSoA[0:1], \
        igSoA->aw[0:nxy], igSoA->ae[0:nxy], igSoA->as[0:nxy], igSoA->an[0:nxy], \
        igSoA->ap[0:nxy], igSoA->su[0:nxy], \
        phi[0:nxy])
    {
        /* initialize phi_new = phi once (optional; helps first iteration) */
        #pragma acc parallel loop present(phi[0:nxy], phi_new[0:nxy])
        for (int idx = 0; idx < nxy; ++idx) phi_new[idx] = phi[idx];

        double *cur  = phi;
        double *next = phi_new;

        for (int it = 0; it < nsw; ++it) {

            #pragma acc parallel loop collapse(2) present(cur[0:nxy], next[0:nxy], \
                igSoA->aw[0:nxy], igSoA->ae[0:nxy], igSoA->as[0:nxy], igSoA->an[0:nxy], \
                igSoA->ap[0:nxy], igSoA->su[0:nxy])
				
                for (int j = 1; j < ny; ++j) {
                    for (int i = 1; i < nx; ++i) {

                        int ij = j*npx + i;

//                        double ap = igSoA->ap[ij];
//                        if (ap == 0.0) { next[ij] = cur[ij]; continue; }

                        double sum_nb =
                              igSoA->aw[ij] * cur[ij - 1]
                            + igSoA->ae[ij] * cur[ij + 1]
                            + igSoA->as[ij] * cur[ij - npx]
                            + igSoA->an[ij] * cur[ij + npx];

                        double phi_gs = (sum_nb + igSoA->su[ij]) / (igSoA->ap[ij]+1.0e-30);
                        next[ij] = (1.0 - wj)*cur[ij] + wj*phi_gs;
                    }
                }

            /* swap pointers (host-side swap is fine; device arrays unchanged) */
            double *tmp = cur; cur = next; next = tmp;
        }

        /* if final result ended up in phi_new, copy back once */
        if (cur != phi) {
            #pragma acc parallel loop present(phi[0:nxy], phi_new[0:nxy])
            for (int idx = 0; idx < nxy; ++idx) phi[idx] = phi_new[idx];
        }
    }

    free(phi_new);
}
void solvergs(struct SoA *igSoA, double *phi, int nx, int ny)
{
    int npx  = igSoA->npx;
    int npy  = igSoA->npy;
    int npxy = npx * npy;
    int nxy  = igSoA->npx * igSoA->npy;

    int nsw = (phi == igSoA->pp) ? nswpp : nswpm;



        for (int it = 0; it < nsw; ++it) {

            #pragma acc parallel loop collapse(2) present(phi[0:nxy],\
             igSoA->aw[0:nxy], igSoA->ae[0:nxy], igSoA->as[0:nxy], igSoA->an[0:nxy], \
                igSoA->ap[0:nxy], igSoA->su[0:nxy])
				
                for (int j = 1; j < ny; ++j) {
                    for (int i = 1; i < nx; ++i) {
                      if ((i+j)%2 == 0) {
                        int ij = j*npx + i;

                        double ap = igSoA->ap[ij];

                        double sum_nb =
                              igSoA->aw[ij] * phi[ij - 1]
                            + igSoA->ae[ij] * phi[ij + 1]
                            + igSoA->as[ij] * phi[ij - npx]
                            + igSoA->an[ij] * phi[ij + npx];
                        phi[ij] = (sum_nb + igSoA->su[ij]) / ap;
                      }  
                    }
                }
                            #pragma acc parallel loop collapse(2) present(phi[0:nxy],\
             igSoA->aw[0:nxy], igSoA->ae[0:nxy], igSoA->as[0:nxy], igSoA->an[0:nxy], \
                igSoA->ap[0:nxy], igSoA->su[0:nxy])
				
                for (int j = 1; j < ny; ++j) {
                    for (int i = 1; i < nx; ++i) {
                      if ((i+j)%2 == 1) {
                        int ij = j*npx + i;

                        double ap = igSoA->ap[ij];

                        double sum_nb =
                              igSoA->aw[ij] * phi[ij - 1]
                            + igSoA->ae[ij] * phi[ij + 1]
                            + igSoA->as[ij] * phi[ij - npx]
                            + igSoA->an[ij] * phi[ij + npx];
                        phi[ij] = (sum_nb + igSoA->su[ij]) / ap;
                      }  
                    }
                }
        }
}

/* -------------------- restriction/prolongation -------------------- */
void restv(struct SoA *igSoA, struct SoA *igSoA1)
{
    int nxf = igSoA->npx, nyf = igSoA->npy;
    int nxc = igSoA1->npx, nyc = igSoA1->npy;

    int nxyf  = nxf * nyf;
    int nxyc  = nxc * nyc;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], \
        igSoA->u[0:nxyf], igSoA->v[0:nxyf], igSoA->p[0:nxyf], igSoA->sc[0:nxyf],\
        igSoA1->u[0:nxyc], igSoA1->v[0:nxyc], igSoA1->p[0:nxyc], igSoA1->sc[0:nxyc],\
        igSoA->areax[0:nxyf], igSoA->areay[0:nxyf],igSoA1->areax[0:nxyc], igSoA1->areay[0:nxyc])
 
        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {

                int ijc = jc*nxc + ic;

                int i = 2*ic;
                int j = 2*jc;

                int ij = j*nxf + i;
                int ijmx = ij - 1;
                int ijmy = ij - nxf;

                igSoA1->u[ijc] = (igSoA->u[ij] * igSoA->areax[ij] +
			       igSoA->u[ijmy] * igSoA->areax[ijmy])/(igSoA1->areax[ijc]+1.0e-30);
	       igSoA1->v[ijc] = (igSoA->v[ij] * igSoA->areay[ij] +
			       igSoA->v[ijmx] * igSoA->areay[ijmx])/(igSoA1->areay[ijc]+1.0e-30);
                igSoA1->p[ijc] = 0.25 * (igSoA->p[ij]    +  igSoA->p[ijmx]
                                       + igSoA->p[ijmy]  +  igSoA->p[ijmy - 1]);
                igSoA1->sc[ijc] = 0.25 * (igSoA->sc[ij]   +  igSoA->sc[ijmx]
                                       +  igSoA->sc[ijmy] +  igSoA->sc[ijmy - 1]);
            }
        }
}

void rest_geom(struct SoA *mySoA)
{
	
// if (nobs>0) {
    
 for (int ig = 1; ig < ngrid; ++ig) {
 
    int nxf = mySoA[ig-1].npx, nyf = mySoA[ig-1].npy;
    int nxc = mySoA[ig].npx, nyc = mySoA[ig].npy;

    int nxyf  = nxf * nyf;
    int nxyc  = nxc * nyc;
	
 
        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 0; ic < nxc-1; ++ic) {
                int ijc = jc*nxc + ic;
                int i = 2*ic;
                int j = 2*jc;
                int ij  = j*nxf + i;
                mySoA[ig].areax[ijc] = mySoA[ig-1].areax[ij] + mySoA[ig-1].areax[ij-nxf];
            }
        }

 
        for (int jc = 0; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {
                int ijc = jc*nxc + ic;
                int i = 2*ic;
                int j = 2*jc;
                int ij  = j*nxf + i;
                mySoA[ig].areay[ijc] = mySoA[ig-1].areay[ij] + mySoA[ig-1].areay[ij-1];
            }
        }
		
    }
 //}
	
 // if (ninlt>0); {
		
   for (int ig = 1; ig < ngrid; ++ig) {
    	
    int nxf = mySoA[ig-1].npx, nyf = mySoA[ig-1].npy;
    int nxc = mySoA[ig].npx, nyc = mySoA[ig].npy;

    int nxyf  = nxf * nyf;
    int nxyc  = nxc * nyc;
		
 
        for (int jc = 1; jc < nyc-1; ++jc) {
                int ijc = jc*nxc + 0;
                int j = 2*jc;
                int ij    = j*nxf + 0;
                mySoA[ig].u[ijc] = 0.5 * (mySoA[ig-1].u[ij] + mySoA[ig-1].u[ij-nxf]);
	       mySoA[ig].v[ijc] = 0.5 * (mySoA[ig-1].v[ij] + mySoA[ig-1].v[ij-1]);
        }		
    }
 // }
  
 // if (noutlt>0); {
	
  for (int ig = 1; ig < ngrid; ++ig) {
	
    int nxf = mySoA[ig-1].npx, nyf = mySoA[ig-1].npy;
    int nxc = mySoA[ig].npx, nyc = mySoA[ig].npy;

    int nxyf  = nxf * nyf;
    int nxyc  = nxc * nyc;

 
        for (int jc = 1; jc < nyc-1; ++jc) {
          
                int ijc = jc*nxc + (nxc-2);
                int i = 2*(nxc-2);
                int j = 2*jc;
                int ij  = j*nxf + i;
                mySoA[ig].u[ijc] = 0.5 * (mySoA[ig-1].u[ij] + mySoA[ig-1].u[ij-nxf]);
	       mySoA[ig].v[ijc] = 0.5 * (mySoA[ig-1].v[ij] + mySoA[ig-1].v[ij-1]);
        }
		
    }
//   }   
}
void restru(struct SoA *igSoA, struct SoA *igSoA1)
{
    coefu(igSoA);

    int nxf = igSoA->npx, nyf = igSoA->npy;
    int nxyf = nxf * nyf;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
         igSoA->rs[0:nxyf])
		
        for (int j = 0; j < nyf; ++j) {
            for (int i = 0; i < nxf; ++i) {
                int ij = j*nxf + i;
                igSoA->rs[ij] =  0.0;
            }
        }
    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->aw[0:nxyf], igSoA->ae[0:nxyf], igSoA->as[0:nxyf], igSoA->an[0:nxyf], igSoA->ap[0:nxyf], \
        igSoA->u[0:nxyf], igSoA->p[0:nxyf], igSoA->rs[0:nxyf], igSoA->resu[0:nxyf],\
        igSoA->areax[0:nxyf], igSoA->areay[0:nxyf])

	
        for (int j = 1; j < nyf-1; ++j) {
            for (int i = 1; i < nxf-2; ++i) {
                int ij = j*nxf + i;

                igSoA->rs[ij] =
                      igSoA->aw[ij] * igSoA->u[ij-1]
                    + igSoA->ae[ij] * igSoA->u[ij+1]
                    + igSoA->as[ij] * igSoA->u[ij-nxf]
                    + igSoA->an[ij] * igSoA->u[ij+nxf]
                    - igSoA->ap[ij] * igSoA->u[ij]
                    + igSoA->areax[ij] * (igSoA->p[ij] - igSoA->p[ij+1])
                    + igSoA->resu[ij] + dpdxm * igSoA->areax[ij] * igSoA->areay[ij];
                if (igSoA->areax[ij] == 0.0) igSoA->rs[ij] = 0.0;
            }
        }

    coefu(igSoA1);

    int nxc = igSoA1->npx, nyc = igSoA1->npy;
    int nxyc = nxc * nyc;

    #pragma acc parallel loop collapse(2) present(igSoA1[0:1], \
         igSoA1->rs[0:nxyc])
		
        for (int j = 0; j < nyc; ++j) {
            for (int i = 0; i < nxc; ++i) {
                int ij = j*nxc + i;
	        igSoA1->rs[ij] =  0.0;
	        }
        }			
				
    #pragma acc parallel loop collapse(2) present(igSoA1[0:1], \
        igSoA1->aw[0:nxyc], igSoA1->ae[0:nxyc], igSoA1->as[0:nxyc], igSoA1->an[0:nxyc], igSoA1->ap[0:nxyc], \
        igSoA1->u[0:nxyc], igSoA1->p[0:nxyc], igSoA1->rs[0:nxyc],\
        igSoA1->areax[0:nxyc], igSoA1->areay[0:nxyc])

        for (int j = 1; j < nyc-1; ++j) {
            for (int i = 1; i < nxc-2; ++i) {
                int ij = j*nxc + i;

                igSoA1->rs[ij] =
                      igSoA1->aw[ij] * igSoA1->u[ij-1]
                    + igSoA1->ae[ij] * igSoA1->u[ij+1]
                    + igSoA1->as[ij] * igSoA1->u[ij-nxc]
                    + igSoA1->an[ij] * igSoA1->u[ij+nxc]
                    - igSoA1->ap[ij] * igSoA1->u[ij]
                    + igSoA1->areax[ij] * (igSoA1->p[ij] - igSoA1->p[ij+1])
                    + dpdxm * igSoA1->areax[ij] * igSoA1->areay[ij];
            }
        }

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], \
        igSoA->rs[0:nxyf], igSoA1->rs[0:nxyc], igSoA1->resu[0:nxyc])
 
        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-2; ++ic) {
                int ijc = jc*nxc + ic;

                int i = 2*ic, j = 2*jc;
                int ij  = j*nxf + i;

                igSoA1->resu[ijc] =(igSoA->rs[ij]    + igSoA->rs[ij - nxf])						
                           + 0.5 * (igSoA->rs[ij - 1]      + igSoA->rs[ij - 1 - nxf] +
                                    igSoA->rs[ij + 1]      + igSoA->rs[ij + 1 - nxf])						
                           -        igSoA1->rs[ijc];
					
            }
        }
}
void restrsc(struct SoA *igSoA, struct SoA *igSoA1)
{
    coefsc(igSoA);

    int nxf = igSoA->npx, nyf = igSoA->npy;
    int nxyf = nxf * nyf;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
         igSoA->rs[0:nxyf])
		
        for (int j = 0; j < nyf; ++j) {
            for (int i = 0; i < nxf; ++i) {
                int ij = j*nxf + i;
                igSoA->rs[ij] =  0.0;
            }
        }
    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->aw[0:nxyf], igSoA->ae[0:nxyf], igSoA->as[0:nxyf], igSoA->an[0:nxyf], igSoA->ap[0:nxyf], \
        igSoA->sc[0:nxyf], igSoA->rs[0:nxyf], igSoA->ressc[0:nxyf])

	
        for (int j = 1; j < nyf-1; ++j) {
            for (int i = 1; i < nxf-1; ++i) {
                int ij = j*nxf + i;

                igSoA->rs[ij] =
                      igSoA->aw[ij] * igSoA->sc[ij-1]
                    + igSoA->ae[ij] * igSoA->sc[ij+1]
                    + igSoA->as[ij] * igSoA->sc[ij-nxf]
                    + igSoA->an[ij] * igSoA->sc[ij+nxf]
                    - igSoA->ap[ij] * igSoA->sc[ij]
                    + igSoA->ressc[ij];
            }
        }

    coefsc(igSoA1);

    int nxc = igSoA1->npx, nyc = igSoA1->npy;
    int nxyc = nxc * nyc;

    #pragma acc parallel loop collapse(2) present(igSoA1[0:1], \
         igSoA1->rs[0:nxyc])
		
        for (int j = 0; j < nyc; ++j) {
            for (int i = 0; i < nxc; ++i) {
                int ij = j*nxc + i;
	        igSoA1->rs[ij] =  0.0;
	        }
        }			
				
    #pragma acc parallel loop collapse(2) present(igSoA1[0:1], \
        igSoA1->aw[0:nxyc], igSoA1->ae[0:nxyc], igSoA1->as[0:nxyc], igSoA1->an[0:nxyc], igSoA1->ap[0:nxyc], \
        igSoA1->sc[0:nxyc], igSoA1->p[0:nxyc], igSoA1->rs[0:nxyc])

        for (int j = 1; j < nyc-1; ++j) {
            for (int i = 1; i < nxc-1; ++i) {
                int ij = j*nxc + i;

                igSoA1->rs[ij] =
                      igSoA1->aw[ij] * igSoA1->sc[ij-1]
                    + igSoA1->ae[ij] * igSoA1->sc[ij+1]
                    + igSoA1->as[ij] * igSoA1->sc[ij-nxc]
                    + igSoA1->an[ij] * igSoA1->sc[ij+nxc]
                    - igSoA1->ap[ij] * igSoA1->sc[ij];
            }
        }

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], \
        igSoA->rs[0:nxyf], igSoA1->rs[0:nxyc], igSoA1->ressc[0:nxyc])
 
        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {
                int ijc = jc*nxc + ic;

                int i = 2*ic, j = 2*jc;
                int ij  = j*nxf + i;

                igSoA1->ressc[ijc] =(igSoA->rs[ij]    + igSoA->rs[ij - nxf]						
                           +        igSoA->rs[ij - 1] + igSoA->rs[ij - 1 - nxf])						
                           -        igSoA1->rs[ijc];
					
            }
        }
}

void restrv(struct SoA *igSoA, struct SoA *igSoA1)
{
    coefv(igSoA);

    int nxf = igSoA->npx, nyf = igSoA->npy;
    int nxyf = nxf * nyf;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->rs[0:nxyf])
		
        for (int j = 0; j < nyf; ++j) {
            for (int i = 0; i < nxf; ++i) {
                int ij = j*nxf + i;
                igSoA->rs[ij] =  0.0;
            }
        }
	
    #pragma acc parallel loop collapse(2) present(igSoA[0:1], \
        igSoA->aw[0:nxyf], igSoA->ae[0:nxyf], igSoA->as[0:nxyf], igSoA->an[0:nxyf], igSoA->ap[0:nxyf], \
        igSoA->v[0:nxyf], igSoA->p[0:nxyf], igSoA->rs[0:nxyf], igSoA->resv[0:nxyf], \
        igSoA->areax[0:nxyf], igSoA->areay[0:nxyf])

        for (int j = 1; j < nyf-2; ++j) {
            for (int i = 1; i < nxf-1; ++i) {
                int ij = j*nxf + i;

                igSoA->rs[ij] =
                      igSoA->aw[ij] * igSoA->v[ij-1]
                    + igSoA->ae[ij] * igSoA->v[ij+1]
                    + igSoA->as[ij] * igSoA->v[ij-nxf]
                    + igSoA->an[ij] * igSoA->v[ij+nxf]
                    - igSoA->ap[ij] * igSoA->v[ij]
                    + igSoA->areay[ij] * (igSoA->p[ij] - igSoA->p[ij + nxf])
                    + igSoA->resv[ij] + dpdym * igSoA->areax[ij] * igSoA->areay[ij];
                if (igSoA->areay[ij] == 0.0) igSoA->rs[ij] = 0.0;
            }
        }

    coefv(igSoA1);

    int nxc = igSoA1->npx, nyc = igSoA1->npy;
    int nxyc = nxc * nyc;

    #pragma acc parallel loop collapse(2) present(igSoA1[0:1], \
        igSoA1->rs[0:nxyc])
		
        for (int j = 0; j < nyc; ++j) {
            for (int i = 0; i < nxc; ++i) {
                int ij = j*nxc + i;
		igSoA1->rs[ij] =  0.0;
	        }
        }

    #pragma acc parallel loop collapse(2) present(igSoA1[0:1], \
        igSoA1->aw[0:nxyc], igSoA1->ae[0:nxyc], igSoA1->as[0:nxyc], igSoA1->an[0:nxyc], igSoA1->ap[0:nxyc], \
        igSoA1->v[0:nxyc], igSoA1->p[0:nxyc], igSoA1->rs[0:nxyc], \
        igSoA1->areax[0:nxyc], igSoA1->areay[0:nxyc])
		
        for (int j = 1; j < nyc-2; ++j) {
            for (int i = 1; i < nxc-1; ++i) {
                int ij = j*nxc + i;

                igSoA1->rs[ij] =
                      igSoA1->aw[ij] * igSoA1->v[ij-1]
                    + igSoA1->ae[ij] * igSoA1->v[ij+1]
                    + igSoA1->as[ij] * igSoA1->v[ij-nxc]
                    + igSoA1->an[ij] * igSoA1->v[ij+nxc]
                    - igSoA1->ap[ij] * igSoA1->v[ij]
                    + igSoA1->areay[ij] * (igSoA1->p[ij] - igSoA1->p[ij + nxc]) 
                    + dpdym * igSoA1->areax[ij] * igSoA1->areay[ij];
            }
        }

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], \
        igSoA->rs[0:nxyf], igSoA1->rs[0:nxyc], igSoA1->resv[0:nxyc])

        for (int jc = 1; jc < nyc-2; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {
                int ijc = jc*nxc + ic;

                int i = 2*ic, j = 2*jc;
                int ij  = j*nxf + i;

                igSoA1->resv[ijc] =(igSoA->rs[ij]    + igSoA->rs[ij - 1])				
                           + 0.5 * (igSoA->rs[ij - nxf]  + igSoA->rs[ij - nxf - 1] +
                                    igSoA->rs[ij + nxf] + igSoA->rs[ij + nxf - 1])						
                                  - igSoA1->rs[ijc];
					
            }
        }
}


/* ******************  prolongations ********************/

void prolu(struct SoA *igSoA, struct SoA *igSoA1)
{

    int nxf = igSoA->npx, nyf = igSoA->npy;
    int nxc = igSoA1->npx, nyc = igSoA1->npy;

    int nxyf = nxf * nyf;
    int nxyc = nxc * nyc;

    #pragma acc parallel loop collapse(2) present(igSoA1[0:1], igSoA1->c[0:nxyc])
        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {
                int ijc = jc*nxc + ic;
                igSoA1->c[ijc] = 0.0;
            }
        }

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], \
        igSoA->u[0:nxyf], igSoA1->u[0:nxyc], igSoA1->c[0:nxyc],\
        igSoA->areax[0:nxyf], igSoA1->areax[0:nxyc])

        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-2; ++ic) {
				
                int ijc = jc*nxc + ic;
                int i = 2*ic;
                int j = 2*jc;				
                int ij  = j*nxf + i;

                double Ru = ( igSoA->u[ij] * igSoA->areax[ij] 
                          +   igSoA->u[ij - nxf] * igSoA->areax[ij-nxf])/
                            ( igSoA1->areax[ijc] + 1.0E-30);
                igSoA1->c[ijc] = igSoA1->u[ijc] - Ru;
            }
        }

    double f1 = 3./4., f2 = 1./4.;
    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], \
        igSoA->u[0:nxyf], igSoA1->c[0:nxyc])
		
        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-2; ++ic) {

                int ijc = jc*nxc + ic;

                int i = 2*ic;
                int j = 2*jc;
                int ij  = j*nxf + i;
                     
                double corr1             = f1 * igSoA1->c[ijc] + f2 * igSoA1->c[ijc + nxc];
                double corr2             = f1 * igSoA1->c[ijc] + f2 * igSoA1->c[ijc - nxc];

                if (igSoA->spu[ij] > 1.0e10) corr1 = 0.0;
                if (igSoA->spu[ij] > 1.0e10) corr2 = 0.0;
      
                #pragma acc atomic update
                igSoA->u[ij]           += corr1;
		#pragma acc atomic update
                igSoA->u[ij - nxf]     += corr2;
		#pragma acc atomic update

		#pragma acc atomic update
                igSoA->u[ij-1]         += 0.5 * corr1;
		#pragma acc atomic update
                igSoA->u[ij - nxf-1]   += 0.5 * corr2;
		#pragma acc atomic update
                igSoA->u[ij+1]         += 0.5 * corr1;
		#pragma acc atomic update
                igSoA->u[ij - nxf+1]   += 0.5 * corr2;

            }
        }
    return;

}

void prolv(struct SoA *igSoA, struct SoA *igSoA1)
{
    int nxf = igSoA->npx, nyf = igSoA->npy;
    int nxc = igSoA1->npx, nyc = igSoA1->npy;

    int nxyf = nxf * nyf;
    int nxyc = nxc * nyc;

    #pragma acc parallel loop collapse(2) present(igSoA1[0:1], igSoA1->c[0:nxyc])

        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {
                int ijc = jc*nxc + ic;
                igSoA1->c[ijc] = 0.0;
            }
        }
		
    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], \
        igSoA->v[0:nxyf], igSoA1->v[0:nxyc], igSoA1->c[0:nxyc], \
        igSoA->areay[0:nxyf], igSoA1->areay[0:nxyc])
  
        for (int jc = 1; jc < nyc-2; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {

                int ijc = jc*nxc + ic;

                int i = 2*ic, j = 2*jc;
                int ij  = j*nxf + i;

                double Rv = (igSoA->v[ij] * igSoA->areay[ij] 
                           + igSoA->v[ij - 1] * igSoA->areay[ij-1])/
                            (igSoA1->areay[ijc] + 1.0E-30);
                igSoA1->c[ijc] = igSoA1->v[ijc] - Rv;
            }
        }

    double f1 = 3./4., f2 = 1./4.;
    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], \
        igSoA->v[0:nxyf], igSoA1->c[0:nxyc])

        for (int jc = 1; jc < nyc-2; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {

                int ijc = jc*nxc + ic;

                int i = 2*ic;
                int j = 2*jc;
                int ij  = j*nxf + i;

                double corr1  = f1 * igSoA1->c[ijc] + f2 * igSoA1->c[ijc + 1];
                double corr2  = f1 * igSoA1->c[ijc] + f2 * igSoA1->c[ijc - 1];

                if (igSoA->spv[ij] > 1.0e10) corr1 = 0.0;
                if (igSoA->spv[ij] > 1.0e10) corr2 = 0.0;
              
        #pragma acc atomic update
                igSoA->v[ij]           += corr1;
		#pragma acc atomic update
                igSoA->v[ij - 1]       += corr2;	
		#pragma acc atomic update
                igSoA->v[ij-nxf]       += 0.5 * corr1;
		#pragma acc atomic update
                igSoA->v[ij-1-nxf]     += 0.5 * corr2;
		#pragma acc atomic update
                igSoA->v[ij+nxf]       += 0.5 * corr1;
		#pragma acc atomic update
                igSoA->v[ij-1+nxf]     += 0.5 * corr2;

            }
        }
    return;
}

void prolp(struct SoA *igSoA, struct SoA *igSoA1)
{
    int nxf = igSoA->npx,  nyf = igSoA->npy;
    int nxc = igSoA1->npx, nyc = igSoA1->npy;

    int nxyf = nxf * nyf;
    int nxyc = nxc * nyc;
	
    #pragma acc parallel loop collapse(2) present(igSoA1[0:1], igSoA1->c[0:nxyc])
 
        for (int jc = 0; jc < nyc; ++jc) {
            for (int ic = 0; ic < nxc; ++ic) {
                int ijc = jc*nxc + ic;
                igSoA1->c[ijc] = 0.0;
            }
        }

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], \
        igSoA->p[0:nxyf], igSoA1->p[0:nxyc], igSoA1->c[0:nxyc])
		
        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {
                int ijc = jc*nxc + ic;

                int i = 2*ic, j = 2*jc;
                int ij  = j*nxf + i;
                int ij1 = ij - nxf;

                double Rp = 0.25 * ( igSoA->p[ij]    + igSoA->p[ij - 1]
                                   + igSoA->p[ij1]   + igSoA->p[ij1 - 1]);

                igSoA1->c[ijc] = igSoA1->p[ijc] - Rp;
                if (igSoA->spu[ij] > 1.0e10) igSoA1->c[ijc] = 0.0;
            }
        }

    double f1 = 9./16., f2 = 3./16., f3 = 3./16., f4 = 1./16.;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], igSoA->p[0:nxyf], igSoA1->c[0:nxyc])
		
        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {

                int ijc  = jc*nxc + ic;

                int i = 2*ic;
                int j = 2*jc;

                int ij  = j*nxf + i;
                
		#pragma acc atomic update
                igSoA->p[ij]     += f1 * igSoA1->c[ijc]        + f2 * igSoA1->c[ijc + 1] +
                                    f3 * igSoA1->c[ijc + nxc]  + f4 * igSoA1->c[ijc + 1 + nxc];
                #pragma acc atomic update
                igSoA->p[ij-1]   += f1 * igSoA1->c[ijc]        + f2 * igSoA1->c[ijc - 1] +
                                    f3 * igSoA1->c[ijc + nxc]  + f4 * igSoA1->c[ijc - 1 + nxc];
                #pragma acc atomic update
                igSoA->p[ij-nxf] += f1 * igSoA1->c[ijc]        + f2 * igSoA1->c[ijc + 1] +
                                    f3 * igSoA1->c[ijc - nxc]  + f4 * igSoA1->c[ijc + 1 - nxc];
                #pragma acc atomic update
                igSoA->p[ij-nxf-1]+= f1 * igSoA1->c[ijc]        + f2 * igSoA1->c[ijc - 1] +
                                     f3 * igSoA1->c[ijc - nxc]  + f4 * igSoA1->c[ijc - 1 - nxc];

            }
        }
    return;

}
void prolsc(struct SoA *igSoA, struct SoA *igSoA1)
{
    int nxf = igSoA->npx,  nyf = igSoA->npy;
    int nxc = igSoA1->npx, nyc = igSoA1->npy;

    int nxyf = nxf * nyf;
    int nxyc = nxc * nyc;
	
    #pragma acc parallel loop collapse(2) present(igSoA1[0:1], igSoA1->c[0:nxyc])
 
        for (int jc = 0; jc < nyc; ++jc) {
            for (int ic = 0; ic < nxc; ++ic) {
                int ijc = jc*nxc + ic;
                igSoA1->c[ijc] = 0.0;
            }
        }

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], \
        igSoA->sc[0:nxyf], igSoA1->sc[0:nxyc], igSoA1->c[0:nxyc])
		
        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {
                int ijc = jc*nxc + ic;

                int i = 2*ic, j = 2*jc;
                int ij  = j*nxf + i;
                int ij1 = ij - nxf;

                double Rp = 0.25 * ( igSoA->sc[ij]    + igSoA->sc[ij - 1]
                                   + igSoA->sc[ij1]   + igSoA->sc[ij1 - 1]);

                igSoA1->c[ijc] = igSoA1->sc[ijc] - Rp;
            }
        }

    double f1 = 9./16., f2 = 3./16., f3 = 3./16., f4 = 1./16.;

    #pragma acc parallel loop collapse(2) present(igSoA[0:1], igSoA1[0:1], igSoA->sc[0:nxyf], igSoA1->c[0:nxyc])
		
        for (int jc = 1; jc < nyc-1; ++jc) {
            for (int ic = 1; ic < nxc-1; ++ic) {

                int ijc  = jc*nxc + ic;

                int i = 2*ic;
                int j = 2*jc;

                int ij  = j*nxf + i;
                
		#pragma acc atomic update
                igSoA->sc[ij]     += f1 * igSoA1->c[ijc]       + f2 * igSoA1->c[ijc + 1] +
                                    f3 * igSoA1->c[ijc + nxc]  + f4 * igSoA1->c[ijc + 1 + nxc];
                #pragma acc atomic update
                igSoA->sc[ij-1]   += f1 * igSoA1->c[ijc]       + f2 * igSoA1->c[ijc - 1] +
                                    f3 * igSoA1->c[ijc + nxc]  + f4 * igSoA1->c[ijc - 1 + nxc];
                #pragma acc atomic update
                igSoA->sc[ij-nxf] += f1 * igSoA1->c[ijc]       + f2 * igSoA1->c[ijc + 1] +
                                    f3 * igSoA1->c[ijc - nxc]  + f4 * igSoA1->c[ijc + 1 - nxc];
                #pragma acc atomic update
                igSoA->sc[ij-nxf-1]+= f1 * igSoA1->c[ijc]       + f2 * igSoA1->c[ijc - 1] +
                                     f3 * igSoA1->c[ijc - nxc]  + f4 * igSoA1->c[ijc - 1 - nxc];

            }
        }
    return;

}

/* -------------------- Tecplot writer (host) -------------------- */

void write_tecplot_2d_cc(const char *fname, const struct SoA *s, int iter,
                         double xl_, double yl_) {
    (void)xl_; (void)yl_; 

    int nx = s->npx;
    int ny = s->npy;

    int I = nx - 2;
    int J = ny - 2;

    int nxy = nx * ny;

    FILE *fp = fopen(fname, "w");
    if (!fp) {
        perror("write_tecplot_2d_cc: fopen");
        return;
    }

    fprintf(fp, "TITLE = \"2D SIMPLE/MG output\"\n");
    fprintf(fp, "VARIABLES = \"X\",\"Y\",\"U\",\"V\",\"P\",\"PP\"\n");
    fprintf(fp, "ZONE T=\"iter_%d\", I=%d, J=%d, DATAPACKING=BLOCK\n",
            iter, I, J);

    
	
        double dx = s->dx, dy = s->dy;

         
            for (int j = 1; j <= J; ++j){
                for (int i = 1; i <= I; ++i) {
                    double x = (i - 0.5) * dx;
                    fprintf(fp, "%.8e\n", x);
				}   
		    }

         
            for (int j = 1; j <= J; ++j){
                for (int i = 1; i <= I; ++i) {
                    double y = (j - 0.5) * dy;
                    fprintf(fp, "%.8e\n", y);
                }
            }
			
    
        for (int j = 1; j <= J; ++j){
            for (int i = 1; i <= I; ++i) {
                int ij  = j*nx + i;
                int ijm = ij - 1;
                double ucc = 0.5 * (s->u[ij] + s->u[ijm]);
                fprintf(fp, "%.8e\n", ucc);
            }
        }

     
 
        for (int j = 1; j <= J; ++j){
            for (int i = 1; i <= I; ++i) {
                int ij  = j*nx + i;
                int ijm = ij - nx;
                double vcc = 0.5 * (s->v[ij] + s->v[ijm]);
                fprintf(fp, "%.8e\n", vcc);
            }
	}

    
        for (int j = 1; j <= J; ++j){
            for (int i = 1; i <= I; ++i) {
                int ij =  j*nx + i;
                fprintf(fp, "%.8e\n", s->p[ij]);
            }
        }

     

        for (int j = 1; j <= J; ++j){
            for (int i = 1; i <= I; ++i) {
                int ij = j*nx + i;
                fprintf(fp, "%.8e\n", s->pp[ij]);
            }
		}

    fclose(fp);
}

void write_vtk_2d_cc(const char *fname, const struct SoA *s, int iter)
{
    // 1. Determine active simulation grid size
    int nx = s->npx;
    int ny = s->npy;
    int I = nx - 2; // Active interior cells along X
    int J = ny - 2; // Active interior cells along Y

    // Target grid dimensions
    const int TARGET_I = 3200;
    const int TARGET_J = 640;
    const int TARGET_N = TARGET_I * TARGET_J;

    // Determine spatial scaling ratio between target and simulation grid
    // E.g., if I=1600, scale_x = 2; if I=3200, scale_x = 1; if I=6400, scale_x = 0.5
    double stride_x = (double)I / (double)TARGET_I;
    double stride_y = (double)J / (double)TARGET_J;

    // Use target dx/dy based on standard domain dimensions
    double target_dx = s->dx * stride_x;
    double target_dy = s->dy * stride_y;

    double sum_vx = 0.0, sum_mag = 0.0;

    FILE *fp = fopen(fname, "w");
    if (!fp) {
        perror("fopen");
        return;
    }

    // 2. VTK Header setup for fixed 3200 x 640 grid
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "2D SIMPLE output iter %d\n", iter);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d 1\n", TARGET_I, TARGET_J);
    fprintf(fp, "POINTS %d float\n", TARGET_N);

    // Write Coordinates for 3200 x 640 grid
    for (int tj = 1; tj <= TARGET_J; ++tj) {
        for (int ti = 1; ti <= TARGET_I; ++ti) {
            double x = (ti - 0.5) * target_dx;
            double y = (tj - 0.5) * target_dy;
            fprintf(fp, "%f %f 0.0\n", x, y);
        }
    }

    fprintf(fp, "\nPOINT_DATA %d\n", TARGET_N);

    // 3. Velocity Vectors (Mapped to target grid)
    fprintf(fp, "VECTORS Velocity float\n");
    for (int tj = 1; tj <= TARGET_J; ++tj) {
        // Map target j -> active simulation cell j
        int j = 1 + (int)((tj - 0.5) * stride_y);
        if (j > J) j = J;

        for (int ti = 1; ti <= TARGET_I; ++ti) {
            // Map target i -> active simulation cell i
            int i = 1 + (int)((ti - 0.5) * stride_x);
            if (i > I) i = I;

            int ij    = j * nx + i;
            int ijm_x = ij - 1;
            int ijm_y = ij - nx;

            double ucc = 0.5 * (s->u[ij] + s->u[ijm_x]);
            double vcc = 0.5 * (s->v[ij] + s->v[ijm_y]);

            sum_vx += ucc;
            sum_mag += sqrt(ucc * ucc + vcc * vcc);

            fprintf(fp, "%f %f 0.0\n", ucc, vcc);
        }
    }

    // 4. Pressure
    fprintf(fp, "\nSCALARS P float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int tj = 1; tj <= TARGET_J; ++tj) {
        int j = 1 + (int)((tj - 0.5) * stride_y);
        if (j > J) j = J;

        for (int ti = 1; ti <= TARGET_I; ++ti) {
            int i = 1 + (int)((ti - 0.5) * stride_x);
            if (i > I) i = I;

            int ij = j * nx + i;
            fprintf(fp, "%f\n", s->p[ij]);
        }
    }

    // 5. Scalar
    fprintf(fp, "\nSCALARS SC float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int tj = 1; tj <= TARGET_J; ++tj) {
        int j = 1 + (int)((tj - 0.5) * stride_y);
        if (j > J) j = J;

        for (int ti = 1; ti <= TARGET_I; ++ti) {
            int i = 1 + (int)((ti - 0.5) * stride_x);
            if (i > I) i = I;

            int ij = j * nx + i;
            fprintf(fp, "%f\n", s->sc[ij]);
        }
    }

    // 6. Obstacle Scalar
    fprintf(fp, "\nSCALARS OBS float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int tj = 1; tj <= TARGET_J; ++tj) {
        int j = 1 + (int)((tj - 0.5) * stride_y);
        if (j > J) j = J;

        for (int ti = 1; ti <= TARGET_I; ++ti) {
            int i = 1 + (int)((ti - 0.5) * stride_x);
            if (i > I) i = I;

            int flag = 0;
            if (nobs > 0) {
                for (int iobs = 0; iobs < nobs; ++iobs) {
                    if (i >= s->isobs[iobs] && i <= s->ieobs[iobs] &&
                        j >= s->jsobs[iobs] && j <= s->jeobs[iobs]) {
                        flag = 1;
                        break;
                    }
                }
            }
            fprintf(fp, "%f\n", flag ? 1.0 : 0.0);
        }
    }

    fclose(fp);

    // Diagnostics calculated over target point count
    int N_sim = I * J;
    printf("Average velocity magnitude at iter %d: %e\n", iter, sum_mag / TARGET_N);
    printf("Average x-velocity at iter %d: %e\n", iter, sum_vx / TARGET_N);
    if (sum_vx != 0.0) {
        printf("Tortuosity at iter %d: %e\n", iter, sum_mag / sum_vx);
    }
}

void write_vtk_2d_cc_square(const char *fname, const struct SoA *s, int iter)
{
    // 1. Determine active simulation grid size
    int nx = s->npx;
    int ny = s->npy;
    int I = nx - 2; // Active interior cells along X
    int J = ny - 2; // Active interior cells along Y

    // Define square target grid dimensions (using TARGET_J for a square aspect ratio)
    const int TARGET_SZ = 640; 
    const int TARGET_N = TARGET_SZ * TARGET_SZ;

    // Define X-range sub-domain from ny to 2*ny
    int x_start = ny;
    int x_end = 2 * ny;
    int x_span = x_end - x_start;

    // Determine spatial scaling ratio for the square sub-domain
    double stride_x = (double)x_span / (double)TARGET_SZ;
    double stride_y = (double)J / (double)TARGET_SZ;

    double target_dx = s->dx * stride_x;
    double target_dy = s->dy * stride_y;

    double sum_vx = 0.0, sum_mag = 0.0;

    FILE *fp = fopen(fname, "w");
    if (!fp) {
        perror("fopen");
        return;
    }

    // 2. VTK Header setup for square grid
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "2D SQUARE output iter %d\n", iter);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d 1\n", TARGET_SZ, TARGET_SZ);
    fprintf(fp, "POINTS %d float\n", TARGET_N);

    // Write Coordinates for square grid
    for (int tj = 1; tj <= TARGET_SZ; ++tj) {
        for (int ti = 1; ti <= TARGET_SZ; ++ti) {
            double x = (ti - 0.5) * target_dx;
            double y = (tj - 0.5) * target_dy;
            fprintf(fp, "%f %f 0.0\n", x, y);
        }
    }

    fprintf(fp, "\nPOINT_DATA %d\n", TARGET_N);

    // 3. Velocity Vectors (Mapped to square target grid with X offset from ny to 2*ny)
    fprintf(fp, "VECTORS Velocity float\n");
    for (int tj = 1; tj <= TARGET_SZ; ++tj) {
        int j = 1 + (int)((tj - 0.5) * stride_y);
        if (j > J) j = J;

        for (int ti = 1; ti <= TARGET_SZ; ++ti) {
            int i = x_start + (int)((ti - 0.5) * stride_x);
            if (i > I) i = I;

            int ij    = j * nx + i;
            int ijm_x = ij - 1;
            int ijm_y = ij - nx;

            double ucc = 0.5 * (s->u[ij] + s->u[ijm_x]);
            double vcc = 0.5 * (s->v[ij] + s->v[ijm_y]);

            sum_vx += ucc;
            sum_mag += sqrt(ucc * ucc + vcc * vcc);

            fprintf(fp, "%f %f 0.0\n", ucc, vcc);
        }
    }

    // 4. Pressure
    fprintf(fp, "\nSCALARS P float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int tj = 1; tj <= TARGET_SZ; ++tj) {
        int j = 1 + (int)((tj - 0.5) * stride_y);
        if (j > J) j = J;

        for (int ti = 1; ti <= TARGET_SZ; ++ti) {
            int i = x_start + (int)((ti - 0.5) * stride_x);
            if (i > I) i = I;

            int ij = j * nx + i;
            fprintf(fp, "%f\n", s->p[ij]);
        }
    }

    // 5. Scalar
    fprintf(fp, "\nSCALARS SC float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int tj = 1; tj <= TARGET_SZ; ++tj) {
        int j = 1 + (int)((tj - 0.5) * stride_y);
        if (j > J) j = J;

        for (int ti = 1; ti <= TARGET_SZ; ++ti) {
            int i = x_start + (int)((ti - 0.5) * stride_x);
            if (i > I) i = I;

            int ij = j * nx + i;
            fprintf(fp, "%f\n", s->sc[ij]);
        }
    }

    // 6. Obstacle Scalar
    fprintf(fp, "\nSCALARS OBS float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int tj = 1; tj <= TARGET_SZ; ++tj) {
        int j = 1 + (int)((tj - 0.5) * stride_y);
        if (j > J) j = J;

        for (int ti = 1; ti <= TARGET_SZ; ++ti) {
            int i = x_start + (int)((ti - 0.5) * stride_x);
            if (i > I) i = I;

            int flag = 0;
            if (nobs > 0) {
                for (int iobs = 0; iobs < nobs; ++iobs) {
                    if (i >= s->isobs[iobs] && i <= s->ieobs[iobs] &&
                        j >= s->jsobs[iobs] && j <= s->jeobs[iobs]) {
                        flag = 1;
                        break;
                    }
                }
            }
            fprintf(fp, "%f\n", flag ? 1.0 : 0.0);
        }
    }

    fclose(fp);

    printf("Square domain - Average velocity magnitude at iter %d: %e\n", iter, sum_mag / TARGET_N);
    printf("Square domain - Average x-velocity at iter %d: %e\n", iter, sum_vx / TARGET_N);
}


/*
void write_vtk_2d_cc(const char *fname, const struct SoA *s, int iter)
{
    int nx = s->npx;
    int ny = s->npy;

    int I = nx - 2;
    int J = ny - 2;

    double sum_vx = 0.0, sum_mag = 0.0;
    
    FILE *fp = fopen(fname, "w");
    if (!fp) {
        perror("fopen");
        return;
    }

    int N = I * J;

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "2D SIMPLE output iter %d\n", iter);
    fprintf(fp, "ASCII\n");

    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d 1\n", I, J);
    fprintf(fp, "POINTS %d float\n", N);

    // Coordinates
    for (int j = 1; j <= J; ++j) {
        for (int i = 1; i <= I; ++i) {
            double x = (i - 0.5) * s->dx;
            double y = (j - 0.5) * s->dy;
            fprintf(fp, "%f %f 0.0\n", x, y);
        }
    }

    fprintf(fp, "\nPOINT_DATA %d\n", N);

    // Velocity vector
    fprintf(fp, "VECTORS Velocity float\n");
    for (int j = 1; j <= J; ++j) {
        for (int i = 1; i <= I; ++i) {
            int ij  = j*nx + i;
            int ijm_x = ij - 1;
            int ijm_y = ij - nx;

            double ucc = 0.5 * (s->u[ij] + s->u[ijm_x]);
            double vcc = 0.5 * (s->v[ij] + s->v[ijm_y]);
            
            sum_vx += ucc;
            sum_mag += sqrt(ucc*ucc + vcc*vcc);

            fprintf(fp, "%f %f 0.0\n", ucc, vcc);
        }
    }

    // Pressure
    fprintf(fp, "\nSCALARS P float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int j = 1; j <= J; ++j) {
        for (int i = 1; i <= I; ++i) {
            int ij = j*nx + i;
            fprintf(fp, "%f\n", s->p[ij]);
        }
    }
    
    // scalar
    fprintf(fp, "\nSCALARS SC float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int j = 1; j <= J; ++j) {
        for (int i = 1; i <= I; ++i) {
            int ij = j*nx + i;
            fprintf(fp, "%f\n", s->sc[ij]);
        }
    }

    // obstacle scalar
    fprintf(fp, "\nSCALARS OBS float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    if (nobs > 0) {
        for (int j = 1; j <= J; ++j) {
            for (int i = 1; i <= I; ++i) {
                int flag = 0;
                for (int iobs = 0; iobs<nobs; ++iobs){
                    if (i >= s->isobs[iobs] && i <= s->ieobs[iobs] &&
                        j >= s->jsobs[iobs] && j <= s->jeobs[iobs]) {
                        flag = 1;
                        break;
                    } 
                }
                if (flag == 1)
                    fprintf(fp, "%f\n", 1.0);
                else
                    fprintf(fp, "%f\n", 0.0);
            }
        }
    }

    fclose(fp);

    printf("Average velocity magnitude at iter %d: %e\n", iter, sum_mag / N);
    printf("Average x-velocity at iter %d: %e\n", iter, sum_vx / N);
    printf("Tortuosity at iter %d: %e\n", iter, (sum_mag) / (sum_vx));
}
*/

/////////////////// Read obstacle from file
void read_obstacle_data(const char *filename) {
    FILE *fptr = fopen(filename, "r");
    if (fptr == NULL) {
        fprintf(stderr, "ERROR: Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    if (fscanf(fptr, "%d", &nobs) != 1) {
        fprintf(stderr, "ERROR: failed reading obstacle line of inputData\n");
        fclose(fptr);
        exit(EXIT_FAILURE);
    }
    
    printf("nobs %d \n", nobs);	
    
    if (nobs > 0) {
        for (int k = 0; k < nobs; ++k) {    
            // Added a quick check here to prevent infinite loops on bad file formats
            if (fscanf(fptr, " %d %d %d %d", &ifobs[k], &ilobs[k], &jfobs[k], &jlobs[k]) != 4) {
                fprintf(stderr, "ERROR: failed reading coordinates for obstacle %d\n", k);
                fclose(fptr);
                exit(EXIT_FAILURE);
            }
        }
    }      
     /*
    if (nobs > 0) {
        for (int k = 0; k < nobs; ++k) {    
            printf(" ifobs = %d, ilobs = %d, jfobs = %d, jlobs = %d\n", ifobs[k], ilobs[k], jfobs[k], jlobs[k]);
        }
    }
	*/
    fclose(fptr);
}


////////////////////////////////////////////////////////////////////
void write_convergence_data_file(const char *fname, double* convergence_data) {
    FILE *fp = fopen(fname, "w");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: Could not open convergence_data file %s for writing.\n", fname);
        return;
    }
    fprintf(fp, "Iteration,convergence_data\n");

    for (int i = 0; i < nvcycle_converged; i++) {
        // i + 1 represents the 1-based iteration number (icycle)
        fprintf(fp, "%d,%e\n", i + 1, convergence_data[i]);
    }
    fclose(fp);
    printf("Successfully wrote %d convergence_data records to %s\n", nvcycle, fname);
}
