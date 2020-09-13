/* lrslib.hpp (vertex enumeration using lexicographic reverse search) */
#define TITLE "lrslib "
#define VERSION "v.7.0 2019.7.10"
#define AUTHOR "*Copyright (C) 1995,2019, David Avis   avis@cs.mcgill.ca "

/* This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   a__int32_t with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
 */
/*Ver 6.1   major change is new lrsnash driver and library coded by Terje Lensberg */
/*Ver 6.0   major change is mplrs wrapper for multithreading coded by Skip Jordan  */
/*Ver 5.0   major change is plrs wrapper for multithreading coded by Gary Roumanis */
/*Ver 4.2*  library version                                      */
/******************************************************************************/
/*  See http://cgm.cs.mcgill.ca/~avis/C/lrs.html for usage instructions         */
/******************************************************************************/

#ifdef PLRS
#include <fstream>
#include <sstream>
#include <string>
#endif

#define lrs_main suf(lrs_main)
#define lrs_overflow suf(lrs_overflow)
#define cache_misses suf(cache_misses)
#define cache_tries suf(cache_tries)
#define checkcobasic suf(checkcobasic)
#define checkindex suf(checkindex)
#define checkpo__int32_t suf(checkpoint)
#define checkredund suf(checkredund)
#define copy_dict suf(copy_dict)
#define dan_selectpivot suf(dan_selectpivot)
#define dict_count suf(dict_count)
#define dict_limit suf(dict_limit)
#define die_gracefully suf(die_gracefully)
#define digits_overflow suf(digits_overflow)
#define getabasis suf(getabasis)
#define getnextoutput suf(getnextoutput)
#define infile suf(infile)
#define infileLen suf(infileLen)
#define infilename suf(infilename)
#define ismin suf(ismin)
#define lexmin suf(lexmin)
#define lprat suf(lprat)
#define lreadrat suf(lreadrat)
#define lrs_alloc_dat suf(lrs_alloc_dat)
#define lrs_alloc_dic suf(lrs_alloc_dic)
#define lrs_cache_to_file suf(lrs_cache_to_file)
#define lrs_cfp suf(lrs_cfp)
#define lrs_checkbound suf(lrs_checkbound)
#define lrs_checkpoint_seconds suf(lrs_checkpoint_seconds)
#define lrs_close suf(lrs_close)
#define lrs_close_outputblock suf(lrs_close_outputblock)
#define lrs_degenerate suf(lrs_degenerate)
#define lrs_dump_state suf(lrs_dump_state)
#define lrs_estimate suf(lrs_estimate)
#define lrs_exit suf(lrs_exit)
#define lrs_stdin_to_file suf(lrs_stdin_to_file)
#define lrs_file_to_cache suf(lrs_file_to_cache)
#define lrs_free_all_memory suf(lrs_free_all_memory)
#define lrs_free_dat suf(lrs_free_dat)
#define lrs_free_dic suf(lrs_free_dic)
#define lrs_free_dic2 suf(lrs_free_dic2)
#define lrs_getdic suf(lrs_getdic)
#define lrs_getfirstbasis suf(lrs_getfirstbasis)
#define lrs_getinput suf(lrs_getinput)
#define lrs_getnextbasis suf(lrs_getnextbasis)
#define lrs_getray suf(lrs_getray)
#define lrs_getsolution suf(lrs_getsolution)
#define lrs_getvertex suf(lrs_getvertex)
#define lrs_global_count suf(lrs_global_count)
#define lrs_global_list suf(lrs_global_list)
#define lrs_init suf(lrs_init)
#define lrs_leaf suf(lrs_leaf)
#define lrs_lpoutput suf(lrs_lpoutput)
#define lrs_open_outputblock suf(lrs_open_outputblock)
#define lrs_printcobasis suf(lrs_printcobasis)
#define lrs_print_header suf(lrs_print_header)
#define lrs_printoutput suf(lrs_printoutput)
#define lrs_printrow suf(lrs_printrow)
#define lrs_printtotals suf(lrs_printtotals)
#define lrs_ratio suf(lrs_ratio)
#define lrs_read_dat suf(lrs_read_dat)
#define lrs_read_dic suf(lrs_read_dic)
#define lrs_return_unexplored suf(lrs_return_unexplored)
#define lrs_set_digits suf(lrs_set_digits)
#define lrs_set_obj suf(lrs_set_obj)
#define lrs_set_obj_mp suf(lrs_set_obj_mp)
#define lrs_set_row suf(lrs_set_row)
#define lrs_set_row_mp suf(lrs_set_row_mp)
#define lrs_solvelp suf(lrs_solvelp)
#define lrs_solve_lp suf(lrs_solve_lp)
#define outfilename suf(outfilename)
#define overflow suf(overflow)
#define pivoting suf(pivoting)
#define phaseone suf(phaseone)
#define pimat suf(pimat)
#define pivot suf(pivot)
#define primalfeasible suf(primalfeasible)
#define printA suf(printA)
#define print_basis suf(print_basis)
#define readfacets suf(readfacets)
#define readlinearity suf(readlinearity)
#define redund_main suf(redund_main)
#define removecobasicindex suf(removecobasicindex)
#define reorder suf(reorder)
#define reorder1 suf(reorder1)
#define rescaledet suf(rescaledet)
#define rescalevolume suf(rescalevolume)
#define resize suf(resize)
#define restartpivots suf(restartpivots)
#define reverse suf(reverse)
#define save_basis suf(save_basis)
#define selectpivot suf(selectpivot)
#define timecheck suf(timecheck)
#define tmpfilename suf(tmpfilename)
#define update suf(update)
#define updatevolume suf(updatevolume)

#ifdef LRSLONG
#define ARITH "lrslong.h"    /* lrs __int32_t integer arithmetic package */
#else
#if defined(GMP) || defined(FLINT)
#define ARITH "lrsgmp.h"     /* lrs wrapper for gmp multiple precsion arithmetic    */
#else
#define ARITH "lrsmp.h"      /* lrs multiple precsion arithmetic    */
#define MP
#endif
#endif

#include ARITH

#ifndef SIGNALS
#include <signal.h>
#include <unistd.h>
#define errcheck(s,e) if ((long)(e)==-1L){  perror(s);exit(1);}
#endif

#define CALLOC(n,s) xcalloc(n,s,__LINE__,__FILE__)

/*********************/
/*global constants   */
/*********************/
#define MAX_LRS_GLOBALS 10000L  /* number of allocated dictionaries */
#define MAXIMIZE 1L         /* maximize the lp  */
#define MINIMIZE 0L         /* maximize the lp  */
#define GE 1L               /* constra__int32_t is >= */
#define EQ 0L               /* constra__int32_t is linearity */

/*************/
/* typedefs  */
/*************/

/******************************************************************************/
/*                   Indexing after initialization                            */
/*               Basis                                    Cobasis             */
/*   ---------------------------------------    ----------------------------- */
/*  |  i  |0|1| .... |lastdv|lastdv+1|...|m|   | j  | 0 | 1 | ... |d-1|  d  | */
/*  |-----|+|+|++++++|++++++|--------|---|-|   |----|---|---|-----|---|+++++| */
/*  |B[i] |0|1| .... |lastdv|lastdv+1|...|m|   |C[j]|m+1|m+2| ... |m+d|m+d+1| */
/*   -----|+|+|++++++|++++++|????????|???|?|    ----|???|???|-----|???|+++++| */
/*                                                                            */
/* Row[i] is row location for B[i]         Col[j] is column location for C[j] */
/*  -----------------------------              -----------------------------  */
/* |   i   |0|1| ..........|m-1|m|            | j    | 0 | 1 | ... |d-1| d  | */
/* |-------|+|-|-----------|---|-|            |------|---|---|---  |---|++++| */
/* |Row[i] |0|1|...........|m-1|m|            |Col[j]| 1 | 2 | ... | d |  0 | */
/* --------|+|*|***********|***|*|             ------|***|***|*****|***|++++| */
/*                                                                            */
/*  + = remains invariant   * = indices may be permuted ? = swapped by pivot  */
/*                                                                            */
/*  m = number of input rows   n= number of input columns                     */
/*  input dimension inputd = n-1 (H-rep) or n (V-rep)                         */
/*  lastdv = inputd-nredundcol  (each redundant column removes a dec. var)    */
/*  working dimension d=lastdv-nlinearity (an input linearity removes a slack) */
/*  obj function in row 0, index 0=B[0]  col 0 has index m+d+1=C[d]           */
/*  H-rep: b-vector in col 0, A matrix in columns 1..n-1                      */
/*  V-rep: col 0 all zero, b-vector in col 1, A matrix in columns 1..n        */
/******************************************************************************/

typedef struct lrs_dic_struct	/* dynamic dictionary data */
{
	lrs_mp_matrix A;
	__int32_t m;			/* A has m+1 rows, row 0 is cost row            */
	__int32_t m_A;           	/* =m or m-d if nonnegative flag set            */
	__int32_t d;			/* A has d+1 columns, col 0 is b-vector         */
	__int32_t d_orig;		/* value of d as A was allocated  (E.G.)        */
	__int32_t lexflag;		/* true if lexmin basis for this vertex         */
	__int32_t depth;			/* depth of basis/vertex in reverse search tree */
	__int32_t i, j;			/* last pivot row and column pivot indices      */
	lrs_mp det;                 /* current determinant of basis                 */
	lrs_mp objnum;		/* objective numerator value                    */
	lrs_mp objden;		/* objective denominator value                  */
	__int32_t *B, *Row;		/* basis, row location indices                  */
	__int32_t *C, *Col;		/* cobasis, column location indices             */
	struct lrs_dic_struct *prev, *next;
} lrs_dic;

typedef struct lrs_dat			/* global problem data   */
{
	lrs_mp_vector Gcd;		/* Gcd of each row of numerators               */
	lrs_mp_vector Lcm;		/* Lcm for each row of input denominators      */
	lrs_mp_vector output;		/* One line of output dimensioned to n         */

	lrs_mp sumdet;		/* sum of determinants                          */
	lrs_mp Nvolume;		/* volume numerator                             */
	lrs_mp Dvolume;		/* volume denominator                           */
	lrs_mp boundn;		/* objective bound numerator                    */
	lrs_mp boundd;		/* objective bound denominator                  */
	__int32_t unbounded;		/* lp unbounded */
	char fname[4096];	/* program name: lrs redund fourier nash        */

	__int32_t *inequality;		/* indices of inequalities corr. to cobasic ind */
	/* initially holds order used to find starting  */
	/* basis, default: m,m-1,...,2,1                */
	__int32_t *facet;		/* cobasic indices for restart in needed        */
	__int32_t *redundcol;		/* holds columns which are redundant            */
	__int32_t *linearity;		/* holds cobasic indices of input linearities   */
	__int32_t *minratio;		/* used for lexicographic ratio test            */
	__int32_t *temparray;		/* for sorting indices, dimensioned to d        */
	__int32_t *isave, *jsave;	/* arrays for estimator, malloc'ed at start     */
	__int32_t inputd;		/* input dimension: n-1 for H-rep, n for V-rep  */

	__int32_t m;      		/* number of rows in input file                 */
	__int32_t n;			/* number of columns in input file              */
	__int32_t lastdv;		/* index of last dec. variable after preproc    */
	/* given by inputd-nredundcol                   */
	__int32_t count[10];		/* count[0]=rays [1]=verts. [2]=base [3]=pivots */
		                /* count[4]=integer vertices                    */

	__int32_t startcount[5];

	__int32_t deepest;		/* max depth ever reached in search             */
	__int32_t nredundcol;		/* number of redundant columns                  */
	__int32_t nlinearity;		/* number of input linearities                  */
	__int32_t totalnodes;		/* count total number of tree nodes evaluated   */
	__int32_t runs;			/* probes for estimate function                 */
	__int32_t seed;			/* seed for random number generator             */
	double cest[10];		/* ests: 0=rays,1=vert,2=bases,3=vol,4=__int32_t vert */
	/**** flags  **********                         */
	__int32_t allbases;		/* TRUE if all bases should be printed          */
	__int32_t bound;                 /* TRUE if upper/lower bound on objective given */
	__int32_t countonly;             /* TRUE if only count totals should be output   */
	__int32_t debug;
	__int32_t dualdeg;		/* TRUE if start dictionary is dual degenerate  */
	__int32_t etrace;		/* turn off debug at basis # strace             */
	__int32_t frequency;		/* frequency to pr__int32_t cobasis indices           */
	__int32_t geometric;		/* TRUE if incident vertex prints after each ray */
	__int32_t getvolume;		/* do volume calculation                        */
	__int32_t givenstart;		/* TRUE if a starting cobasis is given          */
	__int32_t homogeneous;		/* TRUE if all entries in column one are zero   */
	__int32_t hull;			/* do convex hull computation if TRUE           */
	__int32_t incidence;             /* pr__int32_t all tight inequalities (vertices/rays) */
	__int32_t lponly;		/* true if only lp solution wanted              */
	__int32_t maxdepth;	/* max depth to search to in treee              */
	__int32_t maximize;		/* flag for LP maximization                     */
	__int32_t maxoutput;     	/* if positive, maximum number of output lines  */
	__int32_t maxcobases;     	/* if positive, after maxcobasis unexplored subtrees reported */
	__int32_t minimize;		/* flag for LP minimization                     */
	__int32_t mindepth;	/* do not backtrack above mindepth              */
	__int32_t nash;                  /* TRUE for computing nash equilibria           */
	__int32_t nonnegative;		/* TRUE if last d constraints are nonnegativity */
	__int32_t polytope;		/* TRUE for facet computation of a polytope     */
	__int32_t printcobasis;		/* TRUE if all cobasis should be printed        */
	__int32_t printslack;		/* TRUE if indices of slack inequal. printed    */
	__int32_t truncate;              /* TRUE: truncate tree when moving from opt vert*/
	__int32_t verbose;               /* FALSE for minimalist output                  */
	__int32_t restart;		/* TRUE if restarting from some cobasis         */
	__int32_t strace;		/* turn on  debug at basis # strace             */
	__int32_t voronoi;		/* compute voronoi vertices by transformation   */
        __int32_t subtreesize;  /* in estimate mode, iterates if cob_est >= subtreesize */
        __int32_t triangulation;     /* TRUE: the cobases printed triangulate the polytope */
        __int32_t newstart;          /* TRUE: lrs is restarted with new arithmetic         */

	/* Variables for saving/restoring cobasis,  db */

	__int32_t id;			/* numbered sequentially */
	char *name;			/* passed by user */

	__int32_t saved_count[5];	        /* Saves Q->count[*] */
	__int32_t *saved_C;
	lrs_mp saved_det;
	lrs_mp saved_sumdet;
	__int32_t saved_depth;
	__int32_t saved_d;

	__int32_t saved_flag;		/* There is something in the saved cobasis */

	/* Variables for cacheing dictionaries, db */
	lrs_dic *Qhead, *Qtail;

}lrs_dat, lrs_dat_p;

/***************************/
/* mplrs hooks and hacks   */
/***************************/
void lrs_open_outputblock(void); /* prevent mplrs output flushes   */
void lrs_close_outputblock(void);/* re-enable mplrs output flushes */
void lrs_return_unexplored(lrs_dic *P,lrs_dat *Q);  /* send cobasis data for unexplored nodes */

#ifndef LRSLONG
void lrs_overflow(__int32_t i);
void lrs_exit(__int32_t i);
#endif

#ifdef PLRS
#define lrs_printf stream_printf
#else
#define lrs_printf fprintf
#endif

#ifdef PLRS
/****************/
/*      PLRS    */
/****************/
#define plrs_cobasisstring suf(plrs_cobasisstring)
void post_output(const char *, const char *);
void open_outputblock(void);
void close_outputblock(void);
void mplrs_cleanstop(__int32_t checkpoint);
void mplrs_emergencystop(const char *);
__int32_t stream_printf(FILE *str, const char *fmt, ...);
#endif


/*******************************/
/* functions  for external use */
/*******************************/

lrs_dat *lrs_alloc_dat (const char *name);	/* allocate for lrs_dat structure "name"       */
lrs_dic *lrs_alloc_dic (lrs_dat * Q);	/* allocate for lrs_dic structure corr. to Q   */
__int32_t lrs_estimate (lrs_dic * P, lrs_dat * Q);	/* get estimates only and returns est number of cobases in subtree */
__int32_t lrs_read_dat (lrs_dat * Q, __int32_t argc, char *argv[]);	/* read header and set up lrs_dat               */
__int32_t lrs_read_dic (lrs_dic * P, lrs_dat * Q);	/* read input and set up problem and lrs_dic    */
__int32_t lrs_checkbound (lrs_dic *P, lrs_dat * Q);  /* TRUE if current objective value exceeds specified bound */
__int32_t lrs_getfirstbasis (lrs_dic ** P_p, lrs_dat * Q, lrs_mp_matrix * Lin,__int32_t no_output); /* gets first basis, FALSE if none,P may get changed if lin. space Lin found  no_output is TRUE supresses output headers P may get changed if lin. space Lin found    */
void lrs_getinput(lrs_dic *P,lrs_dat *Q,__int32_t *num,__int32_t *den, __int32_t m, __int32_t d); /* reads input matrix b A in lrs/cdd format */
__int32_t lrs_getnextbasis (lrs_dic ** dict_p, lrs_dat * Q, __int32_t prune); /* gets next lrs tree basis, FALSE if none backtrack if prune is TRUE                   */
__int32_t lrs_getsolution (lrs_dic * P, lrs_dat * Q, lrs_mp_vector output, __int32_t col);
__int32_t lrs_getray (lrs_dic * P, lrs_dat * Q, __int32_t col, __int32_t comment, lrs_mp_vector output);
__int32_t lrs_getvertex (lrs_dic * P, lrs_dat * Q, lrs_mp_vector output);
void lrs_close (char *name);	/* close lrs lib program "name"                 */
__int32_t lrs_init (char *name);	/* initialize lrslib and arithmetic package for prog "name" */
void lrs_lpoutput(lrs_dic * P,lrs_dat * Q, lrs_mp_vector output); /* pr__int32_t LP primal and dual solutions */
void lrs_printcobasis (lrs_dic * P, lrs_dat * Q, __int32_t col); /* pr__int32_t cobasis for column col(verted or ray)  */
void lrs_print_header(char *name);
void lrs_printoutput (lrs_dat * Q, lrs_mp_vector output); /* pr__int32_t output array                           */
void lrs_printrow (char name[], lrs_dat * Q, lrs_mp_vector output, __int32_t rowd); /*pr__int32_t row of A matrix in output[0..rowd]      */
void lrs_printsol (lrs_dic * P, lrs_dat * Q, __int32_t col, __int32_t comment);	/* pr__int32_t out solution from col, comment= 0=normal,-1=geometric ray,1..inputd=linearity */
void lrs_printtotals (lrs_dic * P, lrs_dat * Q);/* pr__int32_t final totals for lrs                   */
__int32_t lrs_set_digits (__int32_t dec_digits );  /* set lrsmp digits to equiv. of decimal dec_digits */
__int32_t lrs_solvelp (lrs_dic * P, lrs_dat * Q, __int32_t maximize);/* solve primal feas LP:TRUE bounded else FALSE */

__int32_t lrs_stdin_to_file (char *name);
__int32_t lrs_file_to_cache(FILE *ifp);
__int32_t lrs_cache_to_file(char *name, char *args);


/*******************************/
/* functions  for internal use */
/*******************************/



/*******************************/
/* basic dictionary functions  */
/*******************************/
__int32_t getabasis (lrs_dic * P, lrs_dat * Q, __int32_t order[]); /* Try to find a starting basis  */
void getnextoutput (lrs_dic * P, lrs_dat * Q, __int32_t i, __int32_t col, lrs_mp out);	/* get A[B[i][col] and copy to out */
__int32_t ismin (lrs_dic * P, lrs_dat * Q, __int32_t r, __int32_t s); /* test if A[r][s] is a min ratio for col s */
__int32_t lexmin (lrs_dic * P, lrs_dat * Q, __int32_t col); /* test A to see if current basis is lexmin       */
void pivot (lrs_dic * P, lrs_dat * Q, __int32_t bas, __int32_t cob);	/* Qpivot routine for array A  */
__int32_t primalfeasible (lrs_dic * P, lrs_dat * Q);	/* Do dual pivots to get primal feasibility       */
__int32_t lrs_ratio (lrs_dic * P, lrs_dat * Q, __int32_t col); /* find lex min. ratio  */
__int32_t removecobasicindex (lrs_dic * P, lrs_dat * Q, __int32_t k);	/* remove C[k] from problem  */
__int32_t restartpivots (lrs_dic * P, lrs_dat * Q); /* restart problem from given cobasis   */
__int32_t reverse (lrs_dic * P, lrs_dat * Q, __int32_t *r, __int32_t s); /* TRUE if B[*r] C[s] is a reverse lex-pos pivot  */
__int32_t selectpivot (lrs_dic * P, lrs_dat * Q, __int32_t *r, __int32_t *s);	/* select pivot indices using lexicographic rule  */
__int32_t dan_selectpivot (lrs_dic * P, lrs_dat * Q, __int32_t *r, __int32_t *s); /* select pivot indices using dantzig-lex rule    */
void update (lrs_dic * P, lrs_dat * Q, __int32_t *i, __int32_t *j); /* update the B,C, LOC arrays after a pivot       */
void updatevolume (lrs_dic * P, lrs_dat * Q); /* rescale determinant and update the volume      */


/*******************************/
/* other functions using P,Q   */
/*******************************/
__int32_t lrs_degenerate (lrs_dic * P, lrs_dat * Q);	/* TRUE if the dictionary is primal degenerate    */
void print_basis (FILE * fp, lrs_dat * Q);
void printA (lrs_dic * P, lrs_dat * Q);	/* raw pr__int32_t of dictionary, bases for debugging   */
void pimat (lrs_dic * P, __int32_t r, __int32_t s, lrs_mp Nt, char name[]); /* pr__int32_t the row r col s of A                     */
__int32_t readfacets (lrs_dat * Q, __int32_t facet[]);	/* read and check facet list                      */
__int32_t readlinearity (lrs_dat * Q);	/* read and check linearity list                  */
void rescaledet (lrs_dic * P, lrs_dat * Q, lrs_mp Vnum, lrs_mp Vden);	/* rescale determinant to get its volume */
void rescalevolume (lrs_dic * P, lrs_dat * Q, lrs_mp Vnum, lrs_mp Vden);	/* adjust volume for dimension          */
__int32_t lrs_leaf(lrs_dic *P, lrs_dat *Q);                    /* true if current dictionary is leaf of reverse search tree  */


/***************************************************/
/* Routines for redundancy checking                */
/***************************************************/
__int32_t checkredund (lrs_dic * P, lrs_dat * Q);/* solve primal lp to check redund of obj fun. returns TRUE if redundant, else FALSE          */
__int32_t checkcobasic (lrs_dic * P, lrs_dat * Q, __int32_t index); /* TRUE if index is cobasic and nondegenerate  FALSE if basic, or degen. cobasic, where it will get pivoted out  */
__int32_t checkindex (lrs_dic * P, lrs_dat * Q, __int32_t index); /* index=0 non-red.,1 red., 2 input linearity NOTE: row is returned all zero if redundant!!  */


/***************************************************/
/* Routines for caching and restoring dictionaries */
/***************************************************/
void lrs_free_dic ( lrs_dic *P, lrs_dat *Q);
void lrs_free_dic2 ( lrs_dic *P, lrs_dat *Q);  /* same as lrs_free_dic but no cache*/
void lrs_free_dat ( lrs_dat *Q);
void copy_dict (lrs_dat * global, lrs_dic * dest, lrs_dic * src);
lrs_dic *alloc_memory (lrs_dat * Q);
lrs_dic * lrs_getdic(lrs_dat *Q);
lrs_dic *resize (lrs_dic * P, lrs_dat * Q);
void lrs_free_all_memory(lrs_dic * P, lrs_dat * Q);

/*******************************/
/* utilities                   */
/*******************************/
void lprat (const char *name, __int32_t Num, __int32_t Den);   /* Pr__int32_t Num/Den without reducing  */
__int32_t lreadrat (__int32_t *Num, __int32_t *Den);   /* read a rational string and convert to __int32_t integers            */
void reorder (__int32_t a[], __int32_t range);	/* reorder array in increasing order with one misplaced element   */
void reorder1 (__int32_t a[], __int32_t b[], __int32_t newone, __int32_t range); /* reorder array a in increasing order with misplaced element newone elements of b go a__int32_t for the ride */

/***************************/
/* lp_solve like functions */
/***************************/
__int32_t lrs_solve_lp(lrs_dic *P, lrs_dat *Q);/* solve lp only for given dictionary */
void lrs_set_row(lrs_dic *P, lrs_dat *Q, __int32_t row, __int32_t num[], __int32_t den[], __int32_t ineq);/* load row i of dictionary from num[]/den[] ineq=GE       */ 
void lrs_set_row_mp(lrs_dic *P, lrs_dat *Q, __int32_t row, lrs_mp_vector num, lrs_mp_vector den, __int32_t ineq);/* same as lrs_set_row except num/den is lrs_mp type       */
void lrs_set_obj(lrs_dic *P, lrs_dat *Q, __int32_t num[], __int32_t den[], __int32_t max); /* set up objective function with coeffs num[]/den[] max=MAXIMIZE or MINIMIZE  */
void lrs_set_obj_mp(lrs_dic *P, lrs_dat *Q, lrs_mp_vector num, lrs_mp_vector den, __int32_t max);/* same as lrs_set_obj but num/den has lrs_mp type */

/* plrs related */

void lrs_open_outputblock(void);
void lrs_close_outputblock(void);
void lrs_return_unexplored(lrs_dic *P,lrs_dat *Q);
void lrs_exit(__int32_t i);
