#include "../src/spreadinterp.h"
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void usage()
{
  printf("usage: spreadtestnd_real [dims [M [N [tol [sort [flags [debug [kerpad [kerevalmeth]]]]]]]]]\n\twhere dims=1,2 or 3\n\tM=# nonuniform pts\n\tN=# uniform pts\n\ttol=requested accuracy\n\tsort=0 (don't sort NU pts), 1 (do), or 2 (maybe sort; default)\n\tflags : expert timing flags (see cnufftspread.h)\n\tdebug=0 (less text out), 1 (more), 2 (lots)\n\tkerpad=0 (no pad to mult of 4), 1 (do)\n\tkerevalmeth=0 (direct), 1 (Horner ppval)\n\nexample: ./spreadtestnd 1 1e6 1e6 1e-6 2 0 1\n");
}

int main(int argc, char* argv[])
/* Test executable for the 1D, 2D, or 3D C++ spreader, both directions.
 * It checks speed, and basic correctness via the grid sum of the result.
 * See usage() for usage.
 *
 * Example: spreadtestnd 3 8e6 8e6 1e-6 2 0 1
 *
 * Compilation (also check ../makefile):
 *    g++ spreadtestnd.cpp ../src/spreadinterp.o ../src/utils.o -o spreadtestnd -fPIC -Ofast -funroll-loops -fopenmp
 *
 * Magland; expanded by Barnett 1/14/17. Better cmd line args 3/13/17
 * indep setting N 3/27/17. parallel rand() & sort flag 3/28/17
 * timing_flags 6/14/17. debug control 2/8/18. sort=2 opt 3/5/18, pad 4/24/18
 */
{
  int d = 3;            // Cmd line args & their defaults:  default #dims
  double tol = 1e-6;    // default (eg 1e-6 has nspread=7)
  BIGINT M = 1e6;       // default # NU pts
  BIGINT roughNg = 1e6; // default # U pts
  int sort = 2;         // spread_sort
  int flags = 0;        // default
  int debug = 0;        // default
  int kerpad = 0;       // default
  int kerevalmeth = 1;  // default: Horner
  if (argc<=1) { usage(); return 0; }
  sscanf(argv[1],"%d",&d);
  if (d<1 || d>3) {
    printf("d must be 1, 2 or 3!\n"); usage(); return 1;
  }
  if (argc>2) {
    double w; sscanf(argv[2],"%lf",&w); M = (BIGINT)w;  // so can read 1e6 right!
    if (M<1) {
      printf("M (# NU pts) must be positive!\n"); usage(); return 1;
    }
  }
  if (argc>3) {
    double w; sscanf(argv[3],"%lf",&w); roughNg = (BIGINT)w;
    if (roughNg<1) {
      printf("N (# U pts) must be positive!\n"); usage(); return 1;
    }
  }
  if (argc>4) {
    sscanf(argv[4],"%lf",&tol);
    if (tol<=0.0) {
      printf("tol must be positive!\n"); usage(); return 1;
    }
  }
  if (argc>5) {
    sscanf(argv[5],"%d",&sort);
    if ((sort!=0) && (sort!=1) && (sort!=2)) {
      printf("sort must be 0, 1 or 2!\n"); usage(); return 1;
    }
  }
  if (argc>6)
    sscanf(argv[6],"%d",&flags);
  if (argc>7) {
    sscanf(argv[7],"%d",&debug);
    if ((debug<0) || (debug>2)) {
      printf("debug must be 0, 1 or 2!\n"); usage(); return 1;
    }
  }
  if (argc>8) {
    sscanf(argv[8],"%d",&kerpad);
    if ((kerpad<0) || (kerpad>1)) {
      printf("kerpad must be 0 or 1!\n"); usage(); return 1;
    }
  }
  if (argc>9) {
    sscanf(argv[9],"%d",&kerevalmeth);
    if ((kerevalmeth<0) || (kerevalmeth>1)) {
      printf("kerevalmeth must be 0 or 1!\n"); usage(); return 1;
    }
  }
  if (argc>10) {
    usage(); return 1;
  }

  int dodir1 = true;                        // control if dir=1 tested at all
  BIGINT N = (BIGINT)round(pow(roughNg,1.0/d));     // Fourier grid size per dim
  BIGINT Ng = (BIGINT)pow(N,d);                     // actual total grid points
  BIGINT N2 = (d>=2) ? N : 1, N3 = (d==3) ? N : 1;  // the y and z grid sizes
  std::vector<FLT> kx(M),ky(1),kz(1),d_nonuniform(2*M);    // NU, Re & Im
  if (d>1) ky.resize(M);                           // only alloc needed coords
  if (d>2) kz.resize(M);
  std::vector<FLT> d_uniform(2*Ng);                        // Re and Im

  spread_opts opts;
  FLT upsampfac = 2.0; // big to test spreader error alone, or 2 or 1.25 for kerevalmeth=1
  setup_spreader(opts,(FLT)tol,upsampfac,kerevalmeth);
  opts.pirange = 0;  // crucial, since the below has NU pts on [0,Nd] in each dim
  opts.chkbnds = 0;  // only for debug, since below code has correct bounds
  opts.debug = debug;   // print more diagnostics?
  opts.sort = sort;
  opts.flags = flags;
  opts.kerpad = kerpad;
  //opts.max_subproblem_size = 1e4; // default is 1e5; minimal difference
  FLT maxerr, ansmod;
  
  // spread a single source, only for reference accuracy check...
  opts.spread_direction=1;
  d_nonuniform[0] = 1.0; d_nonuniform[1] = 0.0;   // unit strength
  kx[0] = ky[0] = kz[0] = N/2.0;                  // at center
  int ier = spreadinterp_real(N,N2,N3,d_uniform.data(),1,kx.data(),ky.data(),kz.data(),d_nonuniform.data(),opts);          // vector::data officially C++11 but works
  if (ier!=0) {
    printf("error when spreading M=1 pt for ref acc check (ier=%d)!\n",ier);
    return ier;
  }
  FLT kersumre = 0.0, kersumim = 0.0;  // sum kernel on uniform grid
  for (BIGINT i=0;i<Ng;++i) {
    kersumre += d_uniform[i]; 
  }

  // now do the large-scale test w/ random sources..
  printf("making random data...\n");
  FLT strre = 0.0, strim = 0.0;          // also sum the strengths
#pragma omp parallel
  {
    unsigned int se=MY_OMP_GET_THREAD_NUM();  // needed for parallel random #s
#pragma omp for schedule(dynamic,1000000) reduction(+:strre,strim)
    for (BIGINT i=0; i<M; ++i) {
      kx[i]=rand01r(&se)*N;
      //kx[i]=2.0*kx[i] - 50.0;      //// to test folding within +-1 period
      if (d>1) ky[i]=rand01r(&se)*N;      // only fill needed coords
      if (d>2) kz[i]=rand01r(&se)*N;
      d_nonuniform[i]=randm11r(&se);
      strre += d_nonuniform[i]; 
    }
  }
  CNTime timer;
  double t;
  if (dodir1) {   // test direction 1 (NU -> U spreading) ......................
    printf("spreadinterp %dD, %.3g U pts, dir=%d, tol=%.3g: nspread=%d\n",d,(double)Ng,opts.spread_direction,tol,opts.nspread);
    timer.start();
    ier = spreadinterp_real(N,N2,N3,d_uniform.data(),M,kx.data(),ky.data(),kz.data(),d_nonuniform.data(),opts);
    t=timer.elapsedsec();
    if (ier!=0) {
      printf("error (ier=%d)!\n",ier);
      return ier;
    } else
      printf("    %.3g NU pts in %.3g s \t%.3g pts/s \t%.3g spread pts/s\n",(double)M,t,M/t,pow(opts.nspread,d)*M/t);
  
    FLT sumre = 0.0, sumim = 0.0;   // check spreading accuracy, wrapping
#pragma omp parallel for reduction(+:sumre,sumim)
    for (BIGINT i=0;i<Ng;++i) {
      sumre += d_uniform[i]; 
    }
    FLT pre = kersumre*strre - kersumim*strim;   // pred ans, complex mult
    FLT pim = kersumim*strre + kersumre*strim;
    FLT maxerr = std::max(fabs(sumre-pre), fabs(sumim-pim));
    FLT ansmod = sqrt(sumre*sumre+sumim*sumim);
    printf("    rel err in total over grid:      %.3g\n",maxerr/ansmod);
    // note this is weaker than below dir=2 test, but is good indicator that
    // periodic wrapping is correct
  }

  // test direction 2 (U -> NU interpolation) ..............................
  printf("making more random NU pts...\n");
  for (BIGINT i=0;i<Ng;++i) {     // unit grid data
    d_uniform[i] = 1.0;
  }
#pragma omp parallel
  {
    unsigned int s=MY_OMP_GET_THREAD_NUM();  // needed for parallel random #s
#pragma omp for schedule(dynamic,1000000)
      for (BIGINT i=0; i<M; ++i) {       // random target pts
        //kx[i]=10+.9*rand01r(&s)*N;   // or if want to keep ns away from edges
	kx[i]=rand01r(&s)*N;
	if (d>1) ky[i]=rand01r(&s)*N;
	if (d>2) kz[i]=rand01r(&s)*N;
      }
  }

  opts.spread_direction=2;
  printf("spreadinterp %dD, %.3g U pts, dir=%d, tol=%.3g: nspread=%d\n",d,(double)Ng,opts.spread_direction,tol,opts.nspread);
  timer.restart();
  ier = spreadinterp_real(N,N2,N3,d_uniform.data(),M,kx.data(),ky.data(),kz.data(),d_nonuniform.data(),opts);
  t=timer.elapsedsec();
  if (ier!=0) {
    printf("error (ier=%d)!\n",ier);
    return 1;
  } else
    printf("    %.3g NU pts in %.3g s \t%.3g pts/s \t%.3g spread pts/s\n",(double)M,t,M/t,pow(opts.nspread,d)*M/t);

  // math test is worst-case error from pred value (kersum) on interp pts:
  maxerr = 0.0;
  for (BIGINT i=0;i<M;++i) {
    FLT err = fabs(d_nonuniform[i]-kersumre);
    if (err>maxerr) maxerr=err;
  }
  ansmod = sqrt(kersumre*kersumre+kersumim*kersumim);
  printf("    max rel err in values at NU pts: %.3g\n",maxerr/ansmod);
  // this is stronger test than for dir=1, since it tests sum of kernel for
  // each NU pt. However, it cannot detect reading
  // from wrong grid pts (they are all unity)

  return 0;
}
