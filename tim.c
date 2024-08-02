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
#include <errno.h>
#include "tplib.h"

int read_tim(FILE *f)
{
  int n;

  n = fread(len, sizeof(int), 2, f);

  if (n < 1)
  {
    if (n < 0)
      perror("Error reading file:");
    return 0;
  }

  if (len[1] > MAXSEQLEN)
  {
    tp_error(TPERR_SEQLEN);
    return 0;
  }

  fread(seqtitle, sizeof(char), len[0], f);
  seqtitle[len[0]] = '\0';

  fread(seqbuf, sizeof(char), len[1], f);
  seqbuf[len[1]] = '\0';

  return len[1];
}
