/* lrslong.c      library code for lrs extended precision arithmetic */
/* Version 4.0, April 13, 2000                                       */
/* Copyright: David Avis 1999, avis@cs.mcgill.ca                     */

/* Derived from prs_single.c ( rational arithmetic for lrs and prs)  */
/* authored by  Ambros Marzetta    Revision 1.2  1998/05/27          */

#ifdef PLRS
#include <sstream>
#include <iostream>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "lrslong.h"

__int32_t lrs_digits;		/* max permitted no. of digits   */
__int32_t lrs_record_digits;		/* this is the biggest acheived so far.     */


#define MAXINPUT 1000		/*max length of any input rational */

void 
gcd (lrs_mp u, lrs_mp v)
     /* Returns u=gcd(u,v) using classic Euclid's algorithm.
        v is destroyed.  Knuth, II, p.320 */

{
#ifndef B128
   __int32_t ul, vl, r;

  ul = labs (*u);
  vl = labs (*v);
#else
  __int32_t ul, vl, r;
  ul = abs32 (*u);
  vl = abs32 (*v);
#endif
  if (ul == 0)
    {
      *u = vl;
      return;
    }

  while (vl != 0)
    {
      r= ul % vl;
      ul = vl;
      vl = r;
    }
  *u = ul;
}				/* gcd */

void 
lcm (lrs_mp a, lrs_mp b)			/* a = least common multiple of a, b; b is preserved */
{
  lrs_mp u, v;
  copy (u, a);
  copy (v, b);
  gcd (u, v);
  exactdivint (a, u, v);		/* v=a/u   a contains remainder = 0 */
  mulint (v, b, a);
}				/* end of lcm */


/***************************************************************/
/*                                                             */
/*     Package of routines for rational arithmetic             */
/*     (Built on top of package for multiprecision arithmetic  */
/*                                                             */
/***************************************************************/

void 
reduce (lrs_mp Na, lrs_mp Da)	/* reduces Na/Da by gcd(Na,Da) */
{
  lrs_mp Nb, Db, Nc, Dc;
  copy (Nb, Na);
  copy (Db, Da);
  storesign (Nb, POS);
  storesign (Db, POS);
  copy (Nc, Na);
  copy (Dc, Da);
  gcd (Nb, Db);			/* Nb is the gcd(Na,Da) */
  exactdivint (Nc, Nb, Na);
  exactdivint (Dc, Nb, Da);
}

void 
reduce__int32_t (lrs_mp Na, lrs_mp Da)	/* divide Na by Da and return */
{
  lrs_mp Temp;
  copy (Temp, Na);
  exactdivint (Temp, Da, Na);
}


__int32_t
comprod (lrs_mp Na, lrs_mp Nb, lrs_mp Nc, lrs_mp Nd)    /* +1 if Na*Nb > Nc*Nd  */
                          /* -1 if Na*Nb < Nc*Nd  */
                          /*  0 if Na*Nb = Nc*Nd  */
{
  lrs_mp mc, md;
  itomp(ZERO,mc); itomp(ZERO,md);   /*just to please some compilers */
  mulint (Na, Nb, mc);
  mulint (Nc, Nd, md);
  if (mp_greater(mc,md))
    return (1);
  if (mp_greater(md,mc))
    return (-1);
  return (0);
}


void 
linrat (lrs_mp Na, lrs_mp Da, __int32_t ka, lrs_mp Nb,  lrs_mp Db, __int32_t kb, lrs_mp Nc, lrs_mp Dc)		
/* computes Nc/Dc = ka*Na/Da  +kb* Nb/Db and reduces answer by gcd(Nc,Dc) */
{
  lrs_mp c;
  itomp(ZERO,c);   /*just to please some compilers */
  mulint (Na, Db, Nc);
  mulint (Da, Nb, c);
  linint (Nc, ka, c, kb);	/* Nc = (ka*Na*Db)+(kb*Da*Nb)  */
  mulint (Da, Db, Dc);		/* Dc =  Da*Db           */
  reduce (Nc, Dc);
}


void 
divrat (lrs_mp Na, lrs_mp Da, lrs_mp Nb, lrs_mp Db, lrs_mp Nc, lrs_mp Dc)	
/* computes Nc/Dc = (Na/Da)  / ( Nb/Db ) and reduces answer by gcd(Nc,Dc) */
{
  mulint (Na, Db, Nc);
  mulint (Da, Nb, Dc);
  reduce (Nc, Dc);
}


void 
mulrat (lrs_mp Na, lrs_mp Da, lrs_mp Nb, lrs_mp Db, lrs_mp Nc, lrs_mp Dc)	
/* computes Nc/Dc = Na/Da  * Nb/Db and reduces answer by gcd(Nc,Dc) */
{
  mulint (Na, Nb, Nc);
  mulint (Da, Db, Dc);
  reduce (Nc, Dc);
}

/***************************************************************/
/*                                                             */
/*     Conversion and I/O functions                            */
/*                                                             */
/***************************************************************/

void 
atomp (const char *s, lrs_mp a)	/*convert string to lrs_mp integer */
     /* based on  atoi KR p.58 */
{
  __int32_t diff, ten, i, sig;
  lrs_mp mpone;
  itomp (ONE, mpone);
  ten = 10;
  for (i = 0; s[i] == ' ' || s[i] == '\n' || s[i] == '\t'; i++);
  /*skip white space */
  sig = POS;
  if (s[i] == '+' || s[i] == '-')	/* sign */
    sig = (s[i++] == '+') ? POS : NEG;
  itomp (0L, a);
  while (s[i] >= '0' && s[i] <= '9')
    {
      diff = s[i] - '0';
      linint (a, ten, mpone, diff);
      i++;
    }
  storesign (a, sig);
  if (s[i])
    {
      fprintf (stderr, "\nIllegal character in number: '%s'\n", s + i);
      exit (1);
    }
}				/* end of atomp */

void 
atoaa (const char *in, char *num, char *den)
     /* convert rational string in to num/den strings */
{
  __int32_t i, j;
  for (i = 0; in[i] != '\0' && in[i] != '/'; i++)
    num[i] = in[i];
  num[i] = '\0';
  den[0] = '\0';
  if (in[i] == '/')
    {
      for (j = 0; in[j + i + 1] != '\0'; j++)
	den[j] = in[i + j + 1];
      den[j] = '\0';
    }
}				/* end of atoaa */

void 
mptodouble (lrs_mp a, double *x)	/* convert lrs_mp to double */
{
  (*x) = (*a);
}

__int32_t
mptoi (lrs_mp a)        /* convert lrs_mp to __int32_t */
{
  return (*a);
}

void 
rattodouble (lrs_mp a, lrs_mp b, double *x)	/* convert lrs_mp rati
						   onal to double */

{
  double y;
  mptodouble (a, &y);
  mptodouble (b, x);
  *x = y / (*x);
}

__int32_t 
readrat (lrs_mp Na, lrs_mp Da)	/* read a rational or integer and convert to lrs_mp */
	       /* returns true if denominator is not one       */
{
  char in[MAXINPUT], num[MAXINPUT], den[MAXINPUT];
  if(fscanf (lrs_ifp, "%s", in)==EOF)
                 {
                   fprintf (lrs_ofp, "\nInvalid input: check you have entered enough data!\n");
                   exit(1);
                 }

  if(!strcmp(in,"end"))          /*premature end of input file */
    {
     return (999);
    }

  atoaa (in, num, den);		/*convert rational to num/dem strings */
  atomp (num, Na);
  if (den[0] == '\0')
    {
      itomp (1L, Da);
      return (FALSE);
    }
  atomp (den, Da);
  return (TRUE);
}

#ifdef PLRS
/* read a rational or integer and convert to lrs_mp with base BASE */
/* returns true if denominator is not one                      */
__int32_t plrs_readrat (lrs_mp Na, lrs_mp Da, const char* rat)
{
  	char in[MAXINPUT], num[MAXINPUT], den[MAXINPUT];
 	strcpy(in, rat);
	atoaa (in, num, den);		/*convert rational to num/dem strings */
	atomp (num, Na);
	if (den[0] == '\0')
	{
		itomp (1L, Da);
		return (FALSE);
	}
	atomp (den, Da);
	return (TRUE);
}

#endif

void 
readmp (lrs_mp a)		/* read an integer and convert to lrs_mp */
{
  __int32_t in;
  if(fscanf (lrs_ifp, "%d", &in)==EOF)
                 {
                   fprintf (lrs_ofp, "\nInvalid integer input");
                   exit(1);
                 }

  itomp (in, a);
}

#ifdef PLRS

string sprat (char name[], lrs_mp Nin, lrs_mp Din)	/*reduce and pr__int32_t Nin/Din  */
{

	//create stream to collect output
	stringstream ss;
	string str;
	
	lrs_mp Nt, Dt;
	copy (Nt, Nin);
	copy (Dt, Din);
	reduce (Nt, Dt);
	if (sign (Nt) != NEG)
		ss<<" ";
#ifndef B128
	ss<<name<<*Nt;
	if (*Dt != 1)
		ss<<"/"<<*Dt;
#else
	__int32_t a = *Nt;
	char buf[64]; /* more than enough for 2^128 */
	buf[63] = '\0';
	__int32_t i;
	if (a < 0)
		ss << '-';
	a = abs32(a);
	for (i=62; a!=0 && i>=0; a=a/10, i--)
		buf[i] = '0' + a%10;
	if (i == 62) /* *Nt == 0 */
		ss << '0';
	else
		ss << buf+i+1;
	if (*Dt != 1)
	{
		ss << '/';
		a = *Dt;
		if (a < 0)
			ss<< '-';
		a = abs32(a);
		for (i=62; a!=0 && i>=0; a=a/10, i--)
			buf[i] = '0' + a%10;
		if (i == 62) /* *Dt == 0, uh oh */
			ss << '0';
		else
			ss << buf+i+1;
	}
#endif
	ss<<" ";
	//pipe stream to single string
	str = ss.str();
	return str;
}

char *cprat (char name[], lrs_mp Nin, lrs_mp Din)
{
        char *ret;
        unsigned __int32_t len;
        __int32_t i, offset=0;
	string s;
        const char *cstr;

	s = sprat(name,Nin,Din);
	cstr = s.c_str();
        len = strlen(cstr);
        ret = (char *)malloc(sizeof(char)*(len+1));

        for (i=0; i+offset<len+1;)
        {
                if (cstr[i+offset]!=' ')
                {
                        ret[i] = cstr[i+offset];
                        i++;
                }
                else /* skip whitespace */
                        offset++;
        }

        return ret;
}

string spmp (char name[], lrs_mp Nt)	/*pr__int32_t the __int32_t precision integer a */
{
	//create stream to collect output
	stringstream ss;
	string str;

	ss<<name;
	if(sign(Nt) != NEG)
		ss<<" ";
#ifndef B128
	ss<<*Nt<<" ";
#else
      {
	__int32_t a = *Nt;
        char buf[64]; /* more than enough for 2^128 */
	buf[63] = '\0';
	__int32_t i;
	if (a < 0)
		ss << '-';
        a = abs32(a);
	for (i=62; a!=0 && i>=0; a=a/10, i--)
		buf[i] = '0' + a%10;
	if (i == 62) /* *Nt == 0 */
		ss << '0';
	else
		ss << buf+i+1;
      }
      ss << " ";
#endif
	//pipe stream to single string
	str = ss.str();
	return str;
}
#endif
void 
pmp (char name[], lrs_mp Nt)
{
  lrs_printf (lrs_ofp, "%s", name);
  if (sign (Nt) != NEG)
    lrs_printf (lrs_ofp, " ");
#ifndef B128
  lrs_printf (lrs_ofp, "%d", *Nt);
#else
  {
    __int32_t lower = *Nt % P10_INT64;
    __int32_t upper = *Nt / P10_INT64;
    if (upper != 0)
      lrs_printf(lrs_ofp, "%d", upper);
    else if (lower < 0)
      lrs_printf(lrs_ofp, "-");
    lrs_printf(lrs_ofp, "%d", abs128(lower));
  }
#endif
  lrs_printf (lrs_ofp, " ");
}

void 
prat (char name[], lrs_mp Nin, lrs_mp Din)
     /*pr__int32_t the __int32_t precision rational Nt/Dt  */
{
  lrs_mp Nt, Dt;
  copy (Nt, Nin);
  copy (Dt, Din);
  reduce (Nt, Dt);
  if (sign (Nt) != NEG)
    lrs_printf (lrs_ofp, " ");
#ifndef B128
  lrs_printf (lrs_ofp, "%s%d", name, *Nt);
  if (*Dt != 1)
    lrs_printf (lrs_ofp, "/%d", *Dt);
#else
  {
    __int32_t lower = *Nt % P10_INT64;
    __int32_t upper = *Nt / P10_INT64;
    lrs_printf(lrs_ofp, "%s", name);
    if (upper != 0)
      lrs_printf(lrs_ofp, "%d", upper);
    else if (lower < 0)
      lrs_printf(lrs_ofp, "-");
    lrs_printf(lrs_ofp, "%d", abs32(lower));
    if (*Dt != 1)
    {
      lower = *Dt % P10_INT64;
      upper = *Dt / P10_INT64;
      lrs_printf(lrs_ofp, "/");
      if (upper != 0)
        lrs_printf(lrs_ofp, "%d", upper);
      if (lower < 0)
        lrs_printf(lrs_ofp, "-");
      lrs_printf(lrs_ofp, "%d", abs32(lower));
    }
  }
#endif
  lrs_printf (lrs_ofp, " ");
}				/* prat */


/***************************************************************/
/*                                                             */
/*     Memory allocation functions                             */
/*                                                             */
/***************************************************************/
lrs_mp_t
lrs_alloc_mp_t ()
 /* dynamic allocation of lrs_mp number */
{
  lrs_mp_t p;
  p=(lrs_mp_t)calloc (1, sizeof (lrs_mp));
  return p;
}


lrs_mp_vector 
lrs_alloc_mp_vector (__int32_t n)
 /* allocate lrs_mp_vector for n+1 lrs_mp numbers */
{
  lrs_mp_vector p;
  __int32_t i;

  p = (lrs_mp_t *) CALLOC ((n + 1), sizeof (lrs_mp_t));
  for (i = 0; i <= n; i++)
    p[i] = (lrs_mp_t) CALLOC (1, sizeof (lrs_mp));

  return p;
}

void
lrs_clear_mp_vector (lrs_mp_vector p, __int32_t n)
/* free space allocated to p */
{
  __int32_t i;
  for (i = 0; i <= n; i++)
    free (p[i]);
  free (p);
}

lrs_mp_matrix 
lrs_alloc_mp_matrix (__int32_t m, __int32_t n)
/* allocate lrs_mp_matrix for m+1 x n+1 lrs_mp numbers */
{
  lrs_mp_matrix a;
  lrs_mp_t araw;
  __int32_t mp_width, row_width;
  __int32_t i, j;

  mp_width = lrs_digits + 1;
  row_width = (n + 1) * mp_width;

  araw = (lrs_mp_t) calloc ((m + 1) * row_width, sizeof (lrs_mp));
  a = (lrs_mp_t **) calloc ((m + 1), sizeof (lrs_mp_vector));

  for (i = 0; i < m + 1; i++)
    {
      a[i] = (lrs_mp_t *) calloc ((n + 1), sizeof (lrs_mp_t));

      for (j = 0; j < n + 1; j++)
	a[i][j] = (araw + i * row_width + j * mp_width);
    }
  return a;
}

void
lrs_clear_mp_matrix (lrs_mp_matrix p, __int32_t m, __int32_t n)
/* free space allocated to lrs_mp_matrix p */
{
  __int32_t i;

/* p[0][0] is araw, the actual matrix storage address */

 free(p[0][0]);

 for (i = 0; i < m + 1; i++)
      free (p[i]);
/* 2015.9.9 memory leak fix */
 free(p);
}

void 
lrs_getdigits (__int32_t *a, __int32_t *b)
{
/* send digit information to user */
  *a = DIG2DEC (ZERO);
  *b = DIG2DEC (ZERO);
  return;
}

void *
xcalloc (__int32_t n, __int32_t s, __int32_t l, char *f)
{
  void *tmp;

  tmp = calloc (n, s);
  if (tmp == 0)
    {
      char buf[200];

      sprintf (buf, "\n\nFatal error on line %d of %s", l, f);
      perror (buf);
      exit (1);
    }
  return tmp;
}

__int32_t 
lrs_mp_init (__int32_t dec_digits, FILE * fpin, FILE * fpout)
/* max number of decimal digits for the computation */
/* __int32_t __int32_t version                                 */
{
#ifndef PLRS
  lrs_ifp = fpin;
  lrs_ofp = fpout;
#endif
  lrs_record_digits = 0;
  lrs_digits =  0;		/* max permitted no. of digits   */
  return TRUE;
}

void 
notimpl (char s[])
{
  fflush (stdout);
  fprintf (stderr, "\nAbnormal Termination  %s\n", s);
  exit (1);
}

/***************************************************************/
/*                                                             */
/*     Misc. functions                                         */
/*                                                             */
/***************************************************************/

void 
reducearray (lrs_mp_vector p, __int32_t n)
/* find largest gcd of p[0]..p[n-1] and divide through */
{
  lrs_mp divisor;
  lrs_mp Temp;
  __int32_t i = 0L;

  while ((i < n) && zero (p[i]))
    i++;
  if (i == n)
    return;

  copy (divisor, p[i]);
  storesign (divisor, POS);
  i++;

  while (i < n)
    {
      if (!zero (p[i]))
	{
	  copy (Temp, p[i]);
	  storesign (Temp, POS);
	  gcd (divisor, Temp);
	}
      i++;
    }

/* reduce by divisor */
  for (i = 0; i < n; i++)
    if (!zero (p[i]))
      reduce__int32_t (p[i], divisor);
}				/* end of reducearray */


__int32_t 
myrandom (__int32_t num, __int32_t nrange)
/* return a random number in range 0..nrange-1 */

{
  __int32_t i;
  i = (num * 401 + 673) % nrange;
  return (i);
}

void 
getfactorial (lrs_mp factorial, __int32_t k)		/* compute k factorial
						   in lrs_mp */
{
  lrs_mp temp;
  __int32_t i;
  itomp (ONE, factorial);
  for (i = 2; i <= k; i++)
    {
      itomp (i, temp);
      mulint (temp, factorial, factorial);
    }
}				/* end of getfactorial */

void 
stringcpy (char *s, char *t)	/*copy t to s pointer version */
{
  while (((*s++) = (*t++)) != '\0');
}
