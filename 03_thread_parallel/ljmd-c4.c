/* 
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * OpenMP parallel baseline c version.
 * optimization 1: apply serial improvements except newtons 3rd law
 * optimization 2: use newton's 3rd law. using 'omp critical' directive.
 * optimization 3: use newton's 3rd law. using 'omp atomic' directive.
 * optimization 4: work around having to use 'omp atomic' 
 *                 by working on per thread force arrays.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information 
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    int nthreads;
};
typedef struct _mdsys mdsys_t;

/* helper function: read a line and then return
   the first string with whitespace stripped off */
static int get_me_a_line(FILE *fp, char *buf)
{
    char tmp[BLEN], *ptr;

    /* read a line and cut of comments and blanks */
    if (fgets(tmp,BLEN,fp)) {
        int i;

        ptr=strchr(tmp,'#');
        if (ptr) *ptr= '\0';
        i=strlen(tmp); --i;
        while(isspace(tmp[i])) {
            tmp[i]='\0';
            --i;
        }
        ptr=tmp;
        while(isspace(*ptr)) {++ptr;}
        i=strlen(ptr);
        strcpy(buf,tmp);
        return 0;
    } else {
        perror("problem reading input");
        return -1;
    }
    return 0;
}
 
/* helper function: zero out an array */
static void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        *d++=0.0;
    }
}

/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2, const double box)
{
    while (x >  boxby2) x -= box;
    while (x < -boxby2) x += box;
    return x;
}

/* compute kinetic energy */
static void ekin(mdsys_t *sys)
{   
    int i;

    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += sys->vx[i]*sys->vx[i] 
            + sys->vy[i]*sys->vy[i] 
            + sys->vz[i]*sys->vz[i];
    }
    sys->ekin *= 0.5*mvsq2e*sys->mass;
    sys->temp  = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}

/* compute forces */
static void force(mdsys_t *sys) 
{
    double epot;
    epot = 0.0;

#if defined(_OPENMP)
#pragma omp parallel reduction(+:epot)
#endif
    {
        double c12,c6,boxby2,rcsq;
        double *fx, *fy, *fz;
        int i, tid, fromidx, toidx;

        /* precompute some constants */
        c12 = 4.0*sys->epsilon*pow(sys->sigma,12.0);
        c6  = 4.0*sys->epsilon*pow(sys->sigma, 6.0);
        rcsq= sys->rcut * sys->rcut;
        boxby2 = 0.5*sys->box;
        epot = 0.0;

        /* let each thread operate on a different
           part of the enlarged force array */
#if defined(_OPENMP)
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        fx=sys->fx + (tid*sys->natoms);
        azzero(fx,sys->natoms);
        fy=sys->fy + (tid*sys->natoms);
        azzero(fy,sys->natoms);
        fz=sys->fz + (tid*sys->natoms);
        azzero(fz,sys->natoms);

        for(i=0; i < (sys->natoms -1); i += sys->nthreads) {
            int ii,j;
            double rx1,ry1,rz1;

            ii = i + tid;
            if (ii >= (sys->natoms -1)) break;
            rx1=sys->rx[ii];
            ry1=sys->ry[ii];
            rz1=sys->rz[ii];
        
            for(j=ii+1; j < sys->natoms; ++j) {
                double rx,ry,rz,rsq;

                /* get distance between particle i and j */
                rx=pbc(rx1 - sys->rx[j], boxby2, sys->box);
                ry=pbc(ry1 - sys->ry[j], boxby2, sys->box);
                rz=pbc(rz1 - sys->rz[j], boxby2, sys->box);
                rsq = rx*rx + ry*ry + rz*rz;

                /* compute force and energy if within cutoff */
                if (rsq < rcsq) {
                    double r6,rinv,ffac;

                    rinv=1.0/rsq;
                    r6=rinv*rinv*rinv;
                    
                    ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                    epot += r6*(c12*r6 - c6);

                    fx[ii] += rx*ffac;
                    fy[ii] += ry*ffac;
                    fz[ii] += rz*ffac;
                    fx[j] -= rx*ffac;
                    fy[j] -= ry*ffac;
                    fz[j] -= rz*ffac;
                }
            }
        }

        /* before reducing the forces, we have to make sure 
           that all threads are done adding to them. */
#if defined (_OPENMP)
#pragma omp barrier
#endif
        /* set equal chunks of index ranges */
        i = 1 + (sys->natoms / sys->nthreads);
        fromidx = tid * i;
        toidx = fromidx + i;
        if (toidx > sys->natoms) toidx = sys->natoms;

        /* reduce forces from threads with tid != 0 into
           the storage of the first thread. since we have
           threads already spawned, we do this in parallel. */
        for (i=1; i < sys->nthreads; ++i) {
            int offs, j;

            offs = i*sys->natoms;

            for (j=fromidx; j < toidx; ++j) {
                sys->fx[j] += sys->fx[offs+j];
                sys->fy[j] += sys->fy[offs+j];
                sys->fz[j] += sys->fz[offs+j];
            }
        }
    }
    sys->epot = epot;
}

/* velocity verlet */
static void velverlet(mdsys_t *sys)
{
    int i;
    double dtmf;
    dtmf = 0.5*sys->dt / mvsq2e / sys->mass;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += dtmf * sys->fx[i];
        sys->vy[i] += dtmf * sys->fy[i];
        sys->vz[i] += dtmf * sys->fz[i];
        sys->rx[i] += sys->dt*sys->vx[i];
        sys->ry[i] += sys->dt*sys->vy[i];
        sys->rz[i] += sys->dt*sys->vz[i];
    }

    /* compute forces and potential energy */
    force(sys);

    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += dtmf * sys->fx[i];
        sys->vy[i] += dtmf * sys->fy[i];
        sys->vz[i] += dtmf * sys->fz[i];
    }
}

/* append data to output. */
static void output(mdsys_t *sys, FILE *erg, FILE *traj)
{
    int i;
    
    printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin+sys->epot);
    for (i=0; i<sys->natoms; ++i) {
        fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i], sys->rz[i]);
    }
}


/* main */
int main(int argc, char **argv) 
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;

#if defined(_OPENMP)
#pragma omp parallel
    {
        if(0 == omp_get_thread_num()) {
            sys.nthreads=omp_get_num_threads();
            printf("Running OpenMP version using %d threads\n", sys.nthreads);
        }
    }
#else
    sys.nthreads=1;
#endif

    /* read input file */
    if(get_me_a_line(stdin,line)) return 1;
    sys.natoms=atoi(line);
    if(get_me_a_line(stdin,line)) return 1;
    sys.mass=atof(line);
    if(get_me_a_line(stdin,line)) return 1;
    sys.epsilon=atof(line);
    if(get_me_a_line(stdin,line)) return 1;
    sys.sigma=atof(line);
    if(get_me_a_line(stdin,line)) return 1;
    sys.rcut=atof(line);
    if(get_me_a_line(stdin,line)) return 1;
    sys.box=atof(line);
    if(get_me_a_line(stdin,restfile)) return 1;
    if(get_me_a_line(stdin,trajfile)) return 1;
    if(get_me_a_line(stdin,ergfile)) return 1;
    if(get_me_a_line(stdin,line)) return 1;
    sys.nsteps=atoi(line);
    if(get_me_a_line(stdin,line)) return 1;
    sys.dt=atof(line);
    if(get_me_a_line(stdin,line)) return 1;
    nprint=atoi(line);

    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.nthreads*sys.natoms*sizeof(double));
    sys.fy=(double *)malloc(sys.nthreads*sys.natoms*sizeof(double));
    sys.fz=(double *)malloc(sys.nthreads*sys.natoms*sizeof(double));

    /* read restart */
    fp=fopen(restfile,"r");
    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
        }
        fclose(fp);
        azzero(sys.fx, sys.nthreads*sys.natoms);
        azzero(sys.fy, sys.nthreads*sys.natoms);
        azzero(sys.fz, sys.nthreads*sys.natoms);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys);
    ekin(&sys);
    
    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if ((sys.nfi % nprint) == 0) 
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet(&sys);
        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done.\n");
    fclose(erg);
    fclose(traj);
    
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);

    return 0;
}
