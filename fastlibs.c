/*******************************************************
 *                                                     *
 *  This code is copyright (C) T J R Cutts 1998        *
 *  tjrc1@mole.bio.cam.ac.uk                           *
 *                                                     *
 *  It may be redistributed and modified, as long as   *
 *  this copyright notice is retained.                 *
 *                                                     *
 *******************************************************/

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "tpatterns.h"

/* Pointers for the beginning of the linked list, and
   the current library */
libent *firstlib = NULL;
libent *currlib;

char iobuffer[IOBUFSZ + 8];
static char *bufptr;
static int lensofar;
int numworkers;
int thisworker;
MPI_Status mpi_status;
int mpi_error;
static char dbtype;

char *fp[6];
int lengths[6];
char frm[6][MAXSEQLEN >> 1];
#ifdef NVSN
char revbuf[MAXSEQLEN];
#endif

int build_libents(char *argv[], int libfile)
{
  FILE *l;
  int r;
  char buf[512], dbtype;
  char *libfilename;
  char t[256];
  char f[256];
  char fmt[4];
  char key[4];

  /* If libfile<0, then user did not specify, so get from FASTLIBS */
  if (libfile < 0)
  {
    libfilename = getenv("FASTLIBS");
    if (libfilename == NULL)
    {
      fprintf(stderr,
              "Error: no library file given, or FASTLIBS not defined\n");
      return 1;
    }
  }
  else
    libfilename = argv[libfile];

  l = fopen(libfilename, "r");

  if (l == NULL)
  {
    sprintf(buf, "%s: Could not open file %s",
            argv[0],
            libfilename);
    perror(buf);
    return 1;
  }

  while (fgets(buf, sizeof(buf), l))
  {
    r = sscanf(buf, "%[^$]$%c%[^/]%s %[0-9T]", t, &dbtype, key, f, fmt);

    if (r >= 3)
    {
#ifdef NVSN
      if (!nucleic_query || (dbtype == FF_NUCLEIC))
#endif
      {
        if (firstlib == NULL)
        {
          firstlib = (libent *)malloc(sizeof(libent));
          currlib = firstlib;
        }
        else
        {
          currlib->next = (libent *)malloc(sizeof(libent));
          currlib = currlib->next;
        }

        if (currlib == NULL)
        {
          tp_error(TPERR_MEM);
          return 1;
        }

        currlib->title = strdup(t);
        currlib->file = strdup(f);

        if ((currlib->title == NULL) ||
            (currlib->file == NULL))
        {
          tp_error(TPERR_MEM);
          return 1;
        }

        currlib->letter = key[0];
        currlib->next = NULL;
        currlib->dbtype = dbtype;

        if ((r == 3) || (key[1] == DB_MULTIFILE))
          currlib->fmt = key[1];
        else
          currlib->fmt = fmt[0];
      }
    }
  }

  fclose(l);
  return 0;
}

/* Frees the linked list */
void freelibents(libent *l)
{
  if (l->next)
    freelibents(l->next);
  if (l->title)
    free(l->title);
  if (l->file)
    free(l->file);
}

void multifile(libent *parent)
{
  FILE *list;
  char *p, buf[512];
  libent tmp;

  list = fopen(parent->file, "r");

  if (list == NULL)
  {
    sprintf(buf, "Could not open %s:", parent->file);
    perror(buf);
    return;
  }

  tmp.dbtype = parent->dbtype;

  while (fgets(buf, sizeof(buf), list) != NULL)
  {
    p = strchr(buf, ' ');
    *p = '\0';
    tmp.file = strdup(buf);
    tmp.title = tmp.file;
    tmp.fmt = p[1];
    searchdb(&tmp);
    free(tmp.file);
  }

  fclose(list);
}

/*
void compare_it(char dbtype, int len, pcre *r)
{
#ifdef NVSN
  if (nucleic_query)
    {
      findpattern(seqtitle, seqbuf, len, 0, r,
      SQ_NOFRAME, SQ_FORWARD);
      rev_comp(seqbuf, revbuf);
      findpattern(seqtitle, revbuf, len, 0, r,
      SQ_NOFRAME, SQ_REVERSE);
    }
  else
#endif
    {
      if (dbtype == FF_PROTEIN)
  findpattern(seqtitle, seqbuf, len, 0, r,
        SQ_NOFRAME, SQ_PROTEIN);
      else
  {
    translate(seqbuf, fp, lengths);
    findpattern(seqtitle, fp[0], lengths[0],
          0, r, 0, SQ_FORWARD);
    findpattern(seqtitle, fp[1], lengths[1],
          0, r, 1, SQ_FORWARD);
    findpattern(seqtitle, fp[2], lengths[2],
          0, r, 2, SQ_FORWARD);
    findpattern(seqtitle, fp[3], lengths[3],
          0, r, 0, SQ_FORWARD);
    findpattern(seqtitle, fp[4], lengths[4],
          0, r, 1, SQ_FORWARD);
    findpattern(seqtitle, fp[5], lengths[5],
          0, r, 2, SQ_FORWARD);
  }
    }
}
*/

void mpi_receive_data(pcre *r)
{

  int length;
  char *seqptr;
  char *tptr;

  fp[0] = frm[0];
  fp[1] = frm[1];
  fp[2] = frm[2];
  fp[3] = frm[3];
  fp[4] = frm[4];
  fp[5] = frm[5];

  for (;;)
  {

    MPI_Recv(seqbuf,
             MAXSEQLEN,
             MPI_CHAR,
             0,
             0,
             MPI_COMM_WORLD,
             &mpi_status);

    MPI_Get_count(&mpi_status,
                  MPI_CHAR,
                  &length);

    /* If we are sent less than four bytes, it's all over */

    if (length < 4)
      return;

    for (lensofar = 0, bufptr = seqbuf;
         lensofar < length;)
    {
      memcpy(len, bufptr, 2 * sizeof(int));

      tptr = bufptr + 2 * sizeof(int);

      seqptr = tptr + len[0] + 1;

      /* Translate and search */

      translate(seqptr, fp, lengths);

      findpattern(tptr, fp[0], lengths[0],
                  0, r, 0, SQ_FORWARD);
      findpattern(tptr, fp[1], lengths[1],
                  0, r, 1, SQ_FORWARD);
      findpattern(tptr, fp[2], lengths[2],
                  0, r, 2, SQ_FORWARD);
      findpattern(tptr, fp[3], lengths[3],
                  0, r, 0, SQ_REVERSE);
      findpattern(tptr, fp[4], lengths[4],
                  0, r, 1, SQ_REVERSE);
      findpattern(tptr, fp[5], lengths[5],
                  0, r, 2, SQ_REVERSE);

      /* Move pointer to beginning of next sequence */
      lensofar += (len[0] + len[1] + 2 * sizeof(int) + 1);

      bufptr = seqbuf + lensofar;
    }
  }
}

void send_mpi_data(void)
{
  static unsigned int nextworker;

  nextworker = (nextworker + 1) % (numworkers - 1);

  MPI_Rsend(seqbuf,
            lensofar,
            MPI_CHAR,
            nextworker + 1,
            0,
            MPI_COMM_WORLD);

  lensofar = 0;
  bufptr = seqbuf;
}

void searchdb(libent *lib)
{

  FILE *db;

  int n;

  char buf[256];

  lensofar = 0;
  bufptr = seqbuf;

  printf("Searching %s...\n", lib->title);

  if (lib->fmt != DB_TIM)
  {
    fprintf(stderr, "I don't understand format %d\n", lib->fmt);
    return;
  }

  db = fopen(lib->file, "r");
  if (db == NULL)
  {
    sprintf(buf, "Could not open %s:", lib->file);
    perror(buf);
    return;
  }

#if IOBUFSZ > 0
  if (setvbuf(db,
              iobuffer,
              _IOFBF,
              IOBUFSZ) != 0)
  {
    fprintf(stderr, "Could not set I/O buffer size\n");
    return;
  }
#endif

  if (lib->file != lib->title)
    printf("\tSearching %s...\n", lib->file);

  for (;;)
  {

    n = fread(len, sizeof(int), 2, db);

    if (n < 1)
    {
      if (lensofar > 0)
        send_mpi_data();

      return;
    }

    if (lensofar + len[0] + len[1] + 2 * sizeof(int) + 1 > MAXSEQLEN)
    {

      if (lensofar > 0)
        send_mpi_data();
      else
      {
        tp_error(TPERR_SEQLEN);
        return;
      }
    }

    memcpy(bufptr, len, 2 * sizeof(int));
    bufptr += 2 * sizeof(int);
    fread(bufptr, sizeof(char), len[0], db);
    bufptr += len[0];
    *bufptr = '\0';
    bufptr++;

    fread(bufptr, sizeof(char), len[1], db);
    bufptr += len[1];

    lensofar += 2 * sizeof(int) + 1 + len[0] + len[1];

    numseqs++;
  }
}
