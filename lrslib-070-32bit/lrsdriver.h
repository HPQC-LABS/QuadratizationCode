#ifndef LRS_DRIVER_H
#define LRS_DRIVER_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


__int32_t lrs_main (int argc, char *argv[]);    /* lrs driver, argv[1]=input file, [argc-1]=output file */
__int32_t lrs1_main(int argc, char *argv[],__int32_t overf, char *tmp)/*__attribute__ ((visibility ("default") ))*/;
__int32_t lrs2_main(int argc, char *argv[],__int32_t overf, char *tmp)/*__attribute__ ((visibility ("default") ))*/;
__int32_t lrsgmp_main(int argc, char *argv[],__int32_t overf, char *tmp)/*__attribute__ ((visibility ("default") ))*/;

__int32_t redund_main (int argc, char *argv[]); /* redund driver, argv[1]=input file, [2]=output file */
__int32_t redund1_main(int argc, char *argv[],__int32_t overf,char *tmp)/*__attribute__ ((visibility ("default") ))*/;
__int32_t redund2_main(int argc, char *argv[],__int32_t overf,char *tmp)/*__attribute__ ((visibility ("default") ))*/;
__int32_t redundgmp_main(int argc, char *argv[],__int32_t overf,char *tmp)/*__attribute__ ((visibility ("default") ))*/;

char** makenewargv(int *argc,char** argv,char* tmp);


extern FILE *lrs_cfp;			/* output file for checkpoint information       */
extern FILE *lrs_ifp;			/* input file pointer       */
extern FILE *lrs_ofp;			/* output file pointer      */
#endif
