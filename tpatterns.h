/*******************************************************
 *                                                     *
 *  This code is copyright (C) T J R Cutts 1998        *
 *  tjrc1@mole.bio.cam.ac.uk                           *
 *                                                     *
 *  It may be redistributed and modified, as long as   *
 *  this copyright notice is retained.                 *
 *                                                     *
 *******************************************************/

#include <limits.h>
#include <stdio.h>
#include <mpi.h>
#include "tplib.h"

#ifndef TNAME
#define TNAME "TPatterns"
#endif

/* NVSN enables code which allows nucleic patterns to be compared
   against nucleic databases.  TPatterns is very slow at doing this,
   and the code is not quite right anyway, so turn the capability
   off. */
#ifdef NVSN
#undef NVSN
#endif

/* This defines the number of () subexpressions allowed within the
   pattern.  TPatterns does not report the matches of such
   subexpressions, although sometimes they are necessary when
   specifying alternatives. */
#define NMATCHES 10

/* Some flags used to tell findpattern() what sort of comparison this
   is.  In fact these only affect the way matches are reported. */
#define SQ_NOFRAME -1
#define SQ_PROTEIN 0
#define SQ_FORWARD 0
#define SQ_REVERSE 1

/* The characters used by FASTA to flag protein or nucleic databases */
#define FF_PROTEIN '0'
#define FF_NUCLEIC '1'

/* The characters FASTA uses to define library formats */
#define DB_FASTA '0'
#define DB_GCG '6'
#define DB_MULTIFILE '@'
#define DB_TIM 'T'

/* Definition for linked list of library files to be searched */
typedef struct
{
   void *next;
   char *title;
   char *file;
   char fmt;
   char letter;
   char dbtype;
} libent;

/* Global variables in main.c */
extern char *results;
#ifdef NVSN
extern int nucleic_query;
#endif
extern int findall;
extern int numseqs;
extern int nummatches;
extern libent *firstlib;
extern libent *currlib;

extern int numworkers;
extern int thisworker;
extern MPI_Status mpi_status;
extern int mpi_error;
void mpi_receive_data(pcre *);

/* Functions in options.c */
void usage(char *, char *);
int options(int, char **, int *);

/* Functions in interact.c */
void doprompts(char *, char *);

/* Functions in fastlibs.c */
typedef void (*searchfunc)(FILE *, char dbtype, pcre *);
void searchdb(libent *);
int build_libents(char **, int);
void freelibents(libent *);

/* Functions and variables in compare.c */
pcre *comp_regexp(char *);
int findpattern(char *, char *, int, int, pcre *, int, int);

/* Functions in gcg.c */
pcre *gcg_init(void);
int read_gcg(FILE *);

/* Functions in tim.c */
int read_tim(FILE *);
