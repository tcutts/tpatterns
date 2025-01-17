/*******************************************************
 *                                                     *
 *  This code is copyright (C) T J R Cutts 1998        *
 *  tjrc1@mole.bio.cam.ac.uk                           *
 *                                                     *
 *  It may be redistributed and modified, as long as   *
 *  this copyright notice is retained.                 *
 *                                                     *
 *******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tplib.h"

/* The following matrix encodes the Universal genetic code.  During
   execution any matrix supplied via the -m option will overwrite it */

txmatrix matrix = {
    {"KNKN", "TTTT", "RSRS", "IIMI"},
    {"QHQH", "PPPP", "RRRR", "LLLL"},
    {"EDED", "AAAA", "GGGG", "VVVV"},
    {"*Y*Y", "SSSS", "*CWC", "LFLF"}};

txmatrix revmatrix;

/* make_revmatrix() creates the reverse complemented matrix by
   swapping bits 4 and 5 with bits 0 and 1, inverting the bits, and
   then keeping only bits 0-5 */

void make_revmatrix(void)
{
  int n, r;

  for (n = 0; n < 64; n++)
  {
    r = (~((n & 0x0C) |
           ((n & 0x03) << 4) |
           ((n >> 4) & 0x03))) &
        0x3F;

    ((char *)revmatrix)[r] = ((char *)matrix)[n];
  }
}

/* compilemx() reads a file of codon -> amino acid translations and
   places them in matrix, above, for later use by the translate()
   function */

int compilemx(char *filename)
{
  char buf[80];
  char path[PATH_MAX];

  FILE *f;
  int a, b, c, res;
  char ta, tb, tc, tx, *r;

  if ((f = fopen(filename, "r")) == NULL)
  {
    /* Don't bother with the environment if an absolute path was
 given */
    if (filename[0] != '/')
    {
      if ((r = getenv("TP_TXMDIR")) != NULL)
      {
        sprintf(path, "%s/%s", r, filename);
        f = fopen(path, "r");
      }
      if (f == NULL)
      {
        sprintf(path, "%s/%s", TFLIB, filename);
        f = fopen(path, "r");
      }
    }
  }

  if (f == NULL)
  {
    perror("Could not open translation matrix file:");
    return 1;
  }

  for (a = 0; a < 4; a++)
  {
    for (b = 0; b < 4; b++)
    {
      for (c = 0; c < 4; c++)
      {
        do
        {
          r = fgets(buf, sizeof(buf), f);
          if (r == NULL)
          {
            /* Premature end of file encountered */
            tp_error(TPERR_MXEOF);
            return 1;
          }
          res = sscanf(buf,
                       "%c%c%c %c",
                       &ta, &tb, &tc, &tx);
        } while (res != 4);

        matrix[a][b][c] = tx;
      }
    }
  }

  fclose(f);
  return 0;
}

/* initbasebits() initialises the lookup table for the translate()
   function */

int basebits[256];

#ifdef NVSN
int comphash[256];
#endif

void initbasebits(void)
{

  /* First of all the hash for N->P translation */

  memset(basebits, 255, sizeof(int) * 256);
  basebits['A'] = 0;
  basebits['C'] = 1;
  basebits['G'] = 2;
  basebits['T'] = 3;
  basebits['U'] = 3;

#ifdef NVSN
  /* And now the hash for reverse complementation */
  memset(comphash, 'x', sizeof(int) * 256);
  comphash['A'] = 'T';
  comphash['B'] = 'V';
  comphash['C'] = 'G';
  comphash['D'] = 'H';
  comphash['G'] = 'C';
  comphash['H'] = 'D';
  comphash['K'] = 'M';
  comphash['M'] = 'K';
  comphash['N'] = 'N';
  comphash['R'] = 'Y';
  comphash['S'] = 'S';
  comphash['T'] = 'A';
  comphash['U'] = 'A';
  comphash['V'] = 'B';
  comphash['X'] = 'X';
  comphash['Y'] = 'Y';
  comphash['.'] = '.';
  comphash['~'] = '~';
#endif
}

/* translate() a nucleic acid sequence in all six reading frames
   simultaneously.  out should be an array of six (char *) pointers,
   all of which must be large enough to hold the translated sequence.
   No sanity checking is performed here. l is an array of six
   intergers, which will take the lengths of the translated sequences.

   Thanks to Peter Benie for help optimising the code */

void translate(char *in, char **out, int *l)
{
  int n, rem, index;
  char *p, *end, *r3, *r4, *r5;
  char *r0 = out[0];
  char *r1 = out[1];
  char *r2 = out[2];

  n = strlen(in);
  end = &in[n] - 5;

  /* How much of a part codon do we have at the end */
  rem = n % 3;

  /* How many residues long will the frame 0 sequence be? */
  n /= 3;

  l[0] = l[1] = l[2] = l[3] = l[4] = l[5] = n;

  /* Position the reverse reading frame pointers at the
     ends of where they will be */

  r3 = &out[3][n];
  r4 = &out[4][n];
  r5 = &out[5][n];

  /* If there aren't two bases at the end, the translations
     will not all be the same length; one or two of them
     will be one residue shorter */
  if (rem != 2)
  {
    r5--;
    l[5]--;
    l[2]--;
    if (rem == 0)
    {
      r4--;
      l[4]--;
      l[1]--;
    }
  }

  *r3 = *r4 = *r5 = '\0';
  r3--;
  r4--;
  r5--;

  index = (basebits[in[0]] << 2) | basebits[in[1]];

  for (p = in; p <= end; p += 3)
  {

    /* Frame 0/3 */

    index = ((index << 2) & 0x3C) | basebits[p[2]];

    if (index & 0x80)
    {
      *r0++ = 'X';
      *r3-- = 'X';
    }
    else
    {
      *r0++ = ((char *)matrix)[index];
      *r3-- = ((char *)revmatrix)[index];
    }

    /* Frame 1/4 */

    index = ((index << 2) & 0x3C) | basebits[p[3]];

    if (index & 0x80)
    {
      *r1++ = 'X';
      *r4-- = 'X';
    }
    else
    {
      *r1++ = ((char *)matrix)[index];
      *r4-- = ((char *)revmatrix)[index];
    }

    /* Frame 2/5 */

    index = ((index << 2) & 0x3C) | basebits[p[4]];

    if (index & 0x80)
    {
      *r2++ = 'X';
      *r5-- = 'X';
    }
    else
    {
      *r2++ = ((char *)matrix)[index];
      *r5-- = ((char *)revmatrix)[index];
    }
  }

  /* We may, at this point, have one or two codons still untranslated,
     so we need to deal with those */

  if (rem < 2)
  {
    index = ((index << 2) & 0x3C) | basebits[p[2]];

    if (index & 0x80)
    {
      *r0++ = 'X';
      *r3-- = 'X';
    }
    else
    {
      *r0++ = ((char *)matrix)[index];
      *r3-- = ((char *)revmatrix)[index];
    }

    if (rem == 1)
    {
      index = ((index << 2) & 0x3C) | basebits[p[3]];

      if (index & 0x80)
      {
        *r1++ = 'X';
        *r4-- = 'X';
      }
      else
      {
        *r1++ = ((char *)matrix)[index];
        *r4-- = ((char *)revmatrix)[index];
      }
    }
  }

  *r0 = *r1 = *r2 = '\0';
}

#ifdef NVSN
/* rev_comp() produces the reverse complement of a nucleic acid
   sequence.  It is only used when comparing a nucleotide pattern
   against a nucleotide db */

void rev_comp(char *in, char *out)
{
  register char *p = in;
  register char *q;

  q = &out[strlen(p)];
  *q-- = '\0';

  while (q >= out)
    *q-- = comphash[*p++];
}
#endif
