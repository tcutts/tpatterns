#include <stdio.h>

#ifndef PATH_MAX
#define PATH_MAX 1024
#endif

#define MAXSEQLEN 350000

/* Buffer size for file IO

   Setting this to zero uses the system default buffers

 */

#ifdef _HIUX_SOURCE
/* The Hitach SR2201 needs a large IO buffer
   for speed */
#define IOBUFSZ 512*1024
#else
#define IOBUFSZ 0
#endif

extern int len[2];
extern char seqbuf[MAXSEQLEN];
extern char seqtitle[1024];

void send_mpi_data(void);

/* Translation matrix */
typedef char txmatrix[4][4][4];
extern txmatrix matrix;

/* Functions in translate.c */
void make_revmatrix(void);
void initbasebits(void);
int compilemx(char *);
void translate(char *, char **, int*);
#ifdef NVSN
void rev_comp(char *, char *);
#endif

/* The type for the sequence reading functions */
typedef int (*readfunc)(FILE *);
typedef void (*initfunc)(void);

/* Functions in fasta.c */

void init_read_fasta(void);
int read_fasta(FILE *);

/* Error stuff in error.c */
typedef enum {
  TPERR_SUCCESS = 0,
  TPERR_MEM,
  TPERR_MXEOF,
  TPERR_REDIR,
  TPERR_SEQLEN
} tperror;

void tp_error(tperror);
