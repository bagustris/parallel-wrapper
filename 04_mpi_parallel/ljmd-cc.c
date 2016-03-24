/* 
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * MPI parallel baseline c version.
 * optimization 1: include serial code optimizations _except_ newton's 3rd law.
 * optimization 2: include serial code optimizations including newton's 3rd law.
 * optimization 3: use cell lists to reduce number of pairs to look at.
 * optimization 4: merge arrays to reduce communication overhead.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* ratio between cutoff radius and length of a cell */
const double cellrat=2.0;
/* number of MD steps between cell list updates */
const int cellfreq=2;

/* structure for cell-list data */
struct _cell {
    int natoms;                 /* number of atoms in this cell */
    int owner;                  /* task/thread id that owns this cell */
    int *idxlist;               /* list of atom indices */
};
typedef struct _cell cell_t;
    
/* structure to hold the complete information 
 * about the MD system */
struct _mdsys {
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp, _pad1;
    double *pos, *vel, *frc, *buf;
    cell_t *clist;
    int *plist, _pad2;
    int natoms, nfi, nsteps;
    int ngrid, ncell, npair, nidx, _pad3;
    double delta;

    int mpisize, mpirank;
    MPI_Comm mpicomm;
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
__attribute__((always_inline))
static void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

/* helper function: apply minimum image convention */
__attribute__((always_inline,pure))
static double pbc(double x, const double boxby2, const double box)
{
    while (x >  boxby2) x -= box;
    while (x < -boxby2) x += box;
    return x;
}

/* build and update cell list */
static void updcells(mdsys_t *sys)
{
    int i, ngrid, ncell, npair, midx, natoms;
    double delta, boxby2, boxoffs;
    boxby2 = 0.5 * sys->box;
    natoms = sys->natoms;
        
    if (sys->clist == NULL) {
        int nidx;
        
        ngrid  = floor(cellrat * sys->box / sys->rcut);
        ncell  = ngrid*ngrid*ngrid;
        delta  = sys->box / ngrid;
        boxoffs= boxby2 - 0.5*delta;
        
        sys->delta = delta;
        sys->ngrid = ngrid;
        sys->ncell = ncell;

        /* allocate cell list storage */
        sys->clist = (cell_t *) malloc(ncell*sizeof(cell_t));
        sys->plist = (int *) malloc(2*ncell*ncell*sizeof(int));

        /* allocate index lists within cell. cell density < 2x avg. density */
        nidx = 2*natoms / ncell + 2;
        nidx = ((nidx/2) + 1) * 2;
        sys->nidx = nidx;
        for (i=0; i<ncell; ++i) {
            sys->clist[i].idxlist = (int *) malloc(nidx*sizeof(int));
        }

        /* build cell pair list, assuming newtons 3rd law. */
        npair = 0;
        for (i=0; i < ncell-1; ++i) {
            int j,k;
            double x1,y1,z1;
            
            k  = i/ngrid/ngrid;
            x1 = k*delta - boxoffs;
            y1 = ((i-(k*ngrid*ngrid))/ngrid)*delta - boxoffs;
            z1 = (i % ngrid)*delta - boxoffs;

            for (j=i+1; j<ncell; ++j) {
                double x2,y2,z2,rx,ry,rz;
                
                k  = j/ngrid/ngrid;
                x2 = k*delta - boxoffs;
                y2 = ((j-(k*ngrid*ngrid))/ngrid)*delta - boxoffs;
                z2 = (j % ngrid)*delta - boxoffs;

                rx=pbc(x1 - x2, boxby2, sys->box);
                ry=pbc(y1 - y2, boxby2, sys->box);
                rz=pbc(z1 - z2, boxby2, sys->box);

                /* check for cells on a line that are too far apart */
                if (fabs(rx) > sys->rcut+delta) continue;
                if (fabs(ry) > sys->rcut+delta) continue;
                if (fabs(rz) > sys->rcut+delta) continue;

                /* check for cells in a plane that are too far apart */
                if (sqrt(rx*rx+ry*ry) > (sys->rcut+sqrt(2.0)*delta)) continue;
                if (sqrt(rx*rx+rz*rz) > (sys->rcut+sqrt(2.0)*delta)) continue;
                if (sqrt(ry*ry+rz*rz) > (sys->rcut+sqrt(2.0)*delta)) continue;

                /* other cells that are too far apart */
                if (sqrt(rx*rx + ry*ry + rz*rz) > (sqrt(3.0)*delta+sys->rcut)) continue;
                
                /* cells are close enough. add to list */
                sys->plist[2*npair  ] = i;
                sys->plist[2*npair+1] = j;
                ++npair;
            }
        }
        sys->npair = npair;
        if (sys->mpirank == 0)
            printf("Cell list has %dx%dx%d=%d cells with %d/%d pairs and "
                   "%d atoms/celllist.\n", ngrid, ngrid, ngrid, sys->ncell, 
                   sys->npair, ncell*(ncell-1)/2, nidx);
    }

    /* reset cell list and sort atoms into cells */
    ncell=sys->ncell;
    delta=sys->delta;
    ngrid=sys->ngrid;
    
    for (i=0; i < sys->ncell; ++i) {
        sys->clist[i].natoms=0;
    }

    boxoffs= boxby2 - 0.5*delta;
    midx=0;
    for (i=0; i < natoms; ++i) {
        int idx,j,k,m,n;
        
        k=floor((pbc(sys->pos[i], boxby2, sys->box)+boxby2)/delta);
        m=floor((pbc(sys->pos[natoms + i], boxby2, sys->box)+boxby2)/delta);
        n=floor((pbc(sys->pos[2*natoms + i], boxby2, sys->box)+boxby2)/delta);
        j = ngrid*ngrid*k+ngrid*m+n;

        idx = sys->clist[j].natoms;
        sys->clist[j].idxlist[idx]=i;
        ++idx;
        sys->clist[j].natoms = idx;
        if (idx > midx) midx=idx;
    }
    if (midx > sys->nidx) {
        if (sys->mpirank == 0)
          printf("overflow in cell list: %d/%d atoms/cells.\n", midx, sys->nidx);
        MPI_Abort(sys->mpicomm, 1);
        exit(1);
    }
    return;
}


/* release cell list storage */
static void free_cell_list(mdsys_t *sys)
{
    int i;
    
    if (sys->clist == NULL) 
        return;
    
    for (i=0; i < sys->ncell; ++i) {
        free(sys->clist[i].idxlist);
    }
    
    free(sys->clist);
    sys->clist = NULL;
    sys->ncell = 0;
}


/* compute kinetic energy */
static void ekin(mdsys_t *sys)
{   
    int i;

    sys->ekin=0.0;
    for (i=0; i< 3*sys->natoms; ++i) {
        sys->ekin += sys->vel[i]*sys->vel[i];
    }
    sys->ekin *= 0.5*mvsq2e*sys->mass;
    sys->temp  = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}

/* compute forces */
static void force(mdsys_t *sys) 
{
    double c12,c6,boxby2,rcsq,epot;
    double *fx,*fy,*fz;
    const double *rx,*ry,*rz;

    int n, natoms;
    natoms = sys->natoms;

    /* zero energy and forces. broadcast positions. */
    epot=0.0;
    azzero(sys->buf,3*natoms);
    MPI_Bcast(sys->pos, 3*natoms, MPI_DOUBLE, 0, sys->mpicomm);

    /* update cell list, if needed. */
    if ((sys->clist == NULL) || ((sys->nfi % cellfreq) == 0)) 
        updcells(sys);

    /* create pointers to the array blocks. since those are 
       const offsets, they can be optimized better than when
       (re-)computing the offsets within the loops. we know
       they are always the same, but the compiler cannot 
       figure it out that easily. */
    fx = sys->buf;
    fy = sys->buf+natoms;
    fz = sys->buf+2*natoms;
    rx = sys->pos;
    ry = sys->pos+natoms;
    rz = sys->pos+2*natoms;

    /* precompute some constants */
    c12 = 4.0*sys->epsilon*pow(sys->sigma,12.0);
    c6  = 4.0*sys->epsilon*pow(sys->sigma, 6.0);
    rcsq= sys->rcut * sys->rcut;
    boxby2 = 0.5*sys->box;

    /* self interaction of atoms in cell */
    for(n=0; n < sys->ncell; n += sys->mpisize) {
        int i,j;
        const cell_t *c1;
        
        i = n + sys->mpirank;
        if (i >= sys->ncell) break;
        c1=sys->clist + i;

        for (j=0; j < c1->natoms-1; ++j) {
            int ii, k;
            double rx1, ry1, rz1;

            ii=c1->idxlist[j];
            rx1=rx[ii];
            ry1=ry[ii];
            rz1=rz[ii];
        
            for(k=j+1; k < c1->natoms; ++k) {
                int jj;
                double rx2,ry2,rz2,rsq;
                
                jj=c1->idxlist[k];
                
                /* get distance between particle i and j */
                rx2=pbc(rx1 - rx[jj], boxby2, sys->box);
                ry2=pbc(ry1 - ry[jj], boxby2, sys->box);
                rz2=pbc(rz1 - rz[jj], boxby2, sys->box);
                rsq = rx2*rx2 + ry2*ry2 + rz2*rz2;
                
                /* compute force and energy if within cutoff */
                if (rsq < rcsq) {
                    double r6,rinv,ffac;

                    rinv=1.0/rsq;
                    r6=rinv*rinv*rinv;
                    
                    ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                    epot += r6*(c12*r6 - c6);

                    fx[ii] += rx2*ffac;
                    fy[ii] += ry2*ffac;
                    fz[ii] += rz2*ffac;
                    fx[jj] -= rx2*ffac;
                    fy[jj] -= ry2*ffac;
                    fz[jj] -= rz2*ffac;
                }
            }
        }
    }    

    /* interaction of atoms in different cells */
    for(n=0; n < sys->npair; n += sys->mpisize) {
        int i,j;
        const cell_t *c1, *c2;
        
        i = n + sys->mpirank;
        if (i >= sys->npair) break;
        c1=sys->clist + sys->plist[2*i];
        c2=sys->clist + sys->plist[2*i+1];
        
        for (j=0; j < c1->natoms; ++j) {
            int ii, k;
            double rx1, ry1, rz1;

            ii=c1->idxlist[j];
            rx1=rx[ii];
            ry1=ry[ii];
            rz1=rz[ii];
        
            for(k=0; k < c2->natoms; ++k) {
                int jj;
                double rx2,ry2,rz2,rsq;
                
                jj=c2->idxlist[k];
                
                /* get distance between particle i and j */
                rx2=pbc(rx1 - rx[jj], boxby2, sys->box);
                ry2=pbc(ry1 - ry[jj], boxby2, sys->box);
                rz2=pbc(rz1 - rz[jj], boxby2, sys->box);
                rsq = rx2*rx2 + ry2*ry2 + rz2*rz2;
                
                /* compute force and energy if within cutoff */
                if (rsq < rcsq) {
                    double r6,rinv,ffac;

                    rinv=1.0/rsq;
                    r6=rinv*rinv*rinv;
                    
                    ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                    epot += r6*(c12*r6 - c6);

                    fx[ii] += rx2*ffac;
                    fy[ii] += ry2*ffac;
                    fz[ii] += rz2*ffac;
                    fx[jj] -= rx2*ffac;
                    fy[jj] -= ry2*ffac;
                    fz[jj] -= rz2*ffac;
                }
            }
        }
    }

    /* sum the partial forces into the sys->fx/y/z arrays */
    MPI_Reduce(sys->buf, sys->frc, 3*natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
}

/* velocity verlet */
static void velverlet(mdsys_t *sys)
{
    int i;
    double dtmf;
    dtmf = 0.5*sys->dt / mvsq2e / sys->mass;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<3*sys->natoms; ++i) {
        sys->vel[i] += dtmf * sys->frc[i];
        sys->pos[i] += sys->dt*sys->vel[i];
    }

    /* compute forces and potential energy */
    force(sys);

    /* second part: propagate velocities by another half step */
    for (i=0; i<3*sys->natoms; ++i) {
        sys->vel[i] += dtmf * sys->frc[i];
    }
}

/* append data to output. */
static void output(mdsys_t *sys, FILE *erg, FILE *traj)
{
    int i,natoms;
    natoms=sys->natoms;
    
    printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin+sys->epot);
    for (i=0; i<natoms; ++i) {
        fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->pos[i], sys->pos[natoms+i], sys->pos[2*natoms+i]);
    }
}


/* main */
int main(int argc, char **argv) 
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;
    double starttime, endtime;

    MPI_Init(&argc, &argv);
    sys.mpicomm = MPI_COMM_WORLD;
    MPI_Comm_size(sys.mpicomm, &sys.mpisize);
    MPI_Comm_rank(sys.mpicomm, &sys.mpirank);
        
    /* only MPI master reads the input file */
    if (sys.mpirank == 0) {
        printf("running in parallel mode with %d MPI tasks\n", sys.mpisize);
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
    }
    
    /* broadcast input to all nodes */
    MPI_Bcast(&sys, sizeof(mdsys_t), MPI_BYTE, 0, sys.mpicomm);
    /* re-initialize the mpi rank that just got overwritten */
    MPI_Comm_rank(sys.mpicomm, &sys.mpirank);
    
    /* allocate memory (needs to be done by all) */
    sys.pos=(double *)malloc(3*sys.natoms*sizeof(double));
    sys.vel=(double *)malloc(3*sys.natoms*sizeof(double));
    sys.frc=(double *)malloc(3*sys.natoms*sizeof(double));
    /* extra storage for temporary force data */
    sys.buf=(double *)malloc(3*sys.natoms*sizeof(double));

    /* read restart, only master. */
    if (sys.mpirank == 0) {

        fp=fopen(restfile,"r");
        if(fp) {
            int natoms;
            natoms=sys.natoms;

            for (i=0; i<natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",sys.pos+i, 
                       sys.pos+natoms+i, sys.pos+2*natoms+i);
            }
            for (i=0; i<natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",sys.vel+i, 
                       sys.vel+natoms+i, sys.vel+2*natoms+i);
            }
            fclose(fp);
            azzero(sys.frc, 3*natoms);
        } else {
            perror("cannot read restart file");
            MPI_Abort(sys.mpicomm, 3);
            return 3;
        }
    } else {
        azzero(sys.vel, 3*sys.natoms);
        azzero(sys.frc, 3*sys.natoms);
    }

    /* create initial cell list */
    sys.clist = NULL;
    sys.plist = NULL;

    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys);
    ekin(&sys);
    
    /* output only for master */
    if (sys.mpirank == 0) {
        erg=fopen(ergfile,"w");
        traj=fopen(trajfile,"w");

        printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
        printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
        output(&sys, erg, traj);
    }
    

    /**************************************************/
    /* main MD loop */
    starttime = MPI_Wtime();
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if (sys.mpirank == 0)
            if ((sys.nfi % nprint) == 0) 
                output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet(&sys);
        ekin(&sys);
    }
    endtime = MPI_Wtime();
    /**************************************************/

    /* clean up: close files, free memory */
    if (sys.mpirank == 0) {
        printf("Simulation Done. Loop time: %.2f seconds\n", endtime-starttime);
        fclose(erg);
        fclose(traj);
    }
    
    free(sys.pos);
    free(sys.vel);
    free(sys.frc);
    free(sys.buf);
    free_cell_list(&sys);

    MPI_Barrier(sys.mpicomm);
    MPI_Finalize();
    return 0;
}
