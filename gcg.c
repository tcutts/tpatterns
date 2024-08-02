/*******************************************************
 *                                                     *
 *  This code is copyright (C) T J R Cutts 1998        *
 *  tjrc1@mole.bio.cam.ac.uk                           *
 *                                                     *
 *  It may be redistributed and modified, as long as   *
 *  this copyright notice is retained.                 *
 *                                                     *
 *******************************************************/

#include <stdlib.h>
#include <string.h>
#include "tplib.h"
#include <pcre.h>

/* Regular expression for the information line in GCG format files */
static pcre *inf;
static pcre_extra *infhints;

/* The regular expression to use */
#define INFOFORMAT "^>+[A-Z]+.*(2BIT|ASCII).*Len: *([0-9]+)"

/* One greater than the number of subexpressions in this */
#define INFOPARTS 3

/* Lookup table from GCG 2BIT bytes to groups of four nucleotides */
const char twobit[256][4] = {
  "CCCC", "CCCT", "CCCA", "CCCG", "CCTC", "CCTT", "CCTA", "CCTG",
  "CCAC", "CCAT", "CCAA", "CCAG", "CCGC", "CCGT", "CCGA", "CCGG",
  "CTCC", "CTCT", "CTCA", "CTCG", "CTTC", "CTTT", "CTTA", "CTTG",
  "CTAC", "CTAT", "CTAA", "CTAG", "CTGC", "CTGT", "CTGA", "CTGG",
  "CACC", "CACT", "CACA", "CACG", "CATC", "CATT", "CATA", "CATG",
  "CAAC", "CAAT", "CAAA", "CAAG", "CAGC", "CAGT", "CAGA", "CAGG",
  "CGCC", "CGCT", "CGCA", "CGCG", "CGTC", "CGTT", "CGTA", "CGTG",
  "CGAC", "CGAT", "CGAA", "CGAG", "CGGC", "CGGT", "CGGA", "CGGG",
  "TCCC", "TCCT", "TCCA", "TCCG", "TCTC", "TCTT", "TCTA", "TCTG",
  "TCAC", "TCAT", "TCAA", "TCAG", "TCGC", "TCGT", "TCGA", "TCGG",
  "TTCC", "TTCT", "TTCA", "TTCG", "TTTC", "TTTT", "TTTA", "TTTG",
  "TTAC", "TTAT", "TTAA", "TTAG", "TTGC", "TTGT", "TTGA", "TTGG",
  "TACC", "TACT", "TACA", "TACG", "TATC", "TATT", "TATA", "TATG",
  "TAAC", "TAAT", "TAAA", "TAAG", "TAGC", "TAGT", "TAGA", "TAGG",
  "TGCC", "TGCT", "TGCA", "TGCG", "TGTC", "TGTT", "TGTA", "TGTG",
  "TGAC", "TGAT", "TGAA", "TGAG", "TGGC", "TGGT", "TGGA", "TGGG",
  "ACCC", "ACCT", "ACCA", "ACCG", "ACTC", "ACTT", "ACTA", "ACTG",
  "ACAC", "ACAT", "ACAA", "ACAG", "ACGC", "ACGT", "ACGA", "ACGG",
  "ATCC", "ATCT", "ATCA", "ATCG", "ATTC", "ATTT", "ATTA", "ATTG",
  "ATAC", "ATAT", "ATAA", "ATAG", "ATGC", "ATGT", "ATGA", "ATGG",
  "AACC", "AACT", "AACA", "AACG", "AATC", "AATT", "AATA", "AATG",
  "AAAC", "AAAT", "AAAA", "AAAG", "AAGC", "AAGT", "AAGA", "AAGG",
  "AGCC", "AGCT", "AGCA", "AGCG", "AGTC", "AGTT", "AGTA", "AGTG",
  "AGAC", "AGAT", "AGAA", "AGAG", "AGGC", "AGGT", "AGGA", "AGGG",
  "GCCC", "GCCT", "GCCA", "GCCG", "GCTC", "GCTT", "GCTA", "GCTG",
  "GCAC", "GCAT", "GCAA", "GCAG", "GCGC", "GCGT", "GCGA", "GCGG",
  "GTCC", "GTCT", "GTCA", "GTCG", "GTTC", "GTTT", "GTTA", "GTTG",
  "GTAC", "GTAT", "GTAA", "GTAG", "GTGC", "GTGT", "GTGA", "GTGG",
  "GACC", "GACT", "GACA", "GACG", "GATC", "GATT", "GATA", "GATG",
  "GAAC", "GAAT", "GAAA", "GAAG", "GAGC", "GAGT", "GAGA", "GAGG",
  "GGCC", "GGCT", "GGCA", "GGCG", "GGTC", "GGTT", "GGTA", "GGTG",
  "GGAC", "GGAT", "GGAA", "GGAG", "GGGC", "GGGT", "GGGA", "GGGG"};

/* tx2bit() translates GCG 2bit into ASCII

   Note that no sanity checking is done; the target buffer
   must be large enough to accomodate the ASCII form of the
   sequence */

void tx2bit(char *bytes, char *target, int blen)
{
  register char *p;
  register char *q;
  unsigned int n;
  
  for (n = 0, p=bytes, q=target;
       n<blen;
       n++, p++, q += 4)
    {
      memcpy(q, twobit[(unsigned char)*p], 4);
    } 
}

/* The regular expression for extracting length and format
   information from each GCG format sequence is going to
   be applied many thousands of times, so it's best to compile
   it just once in advance */

pcre *gcg_init(void)
{
  const char *err;
  int errptr;

  inf = pcre_compile(INFOFORMAT, 0, &err, &errptr);

  /* In PCRE 1.09, this study is ineffectual.  However, later
     versions may add some optimisations, so this doesn't hurt. */

  if (inf != NULL)
    infhints = pcre_study(inf, 0, &err);

  return inf;
}

int read_gcg(FILE *f, int fd)
{
  static char tmpbuf[MAXSEQLEN>>2];
  static char info[1024];
  char *s, *p;
  int match[INFOPARTS<<1];
  int n, length, bytes, isascii;
  char c;

  while (!feof(f))
    {
    /* Get the info line first */
      p = fgets(info, sizeof(info), f);
      if (p)
	{
	  /* Got info line, title is next */
	  fgets(seqtitle, sizeof(seqtitle), f);

	  /* Interpret the info line */
	  n = pcre_exec(inf, infhints, info,
			strlen(info), 0, match, INFOPARTS*2);

	  if (n>=0)
	    {
	      /* Is this ASCII or 2 Bit format? */
	      isascii = (info[match[2]] == 'A') ? 1 : 0;
	    
	      /* How many bytes do we need to get? */
	      c = info[match[5]];
	      info[match[5]] = '\0';
	      length = atoi(&info[match[4]]);
	      info[match[5]] = c;

	      if (isascii)
		{
		  s = seqbuf;
		  bytes = length + 1;
		}
	      else
		{
		  s = tmpbuf;
		  bytes = (length >> 2) + 2;
		}

	      if (length>MAXSEQLEN)
		{
		  tp_error(TPERR_SEQLEN);
		  return 0;
		}

	      /* Read the actual sequence */
	      n = fread(s, sizeof(char), bytes, f);
	      if (n==0)
		{
		  perror("File read error: ");
		  return 0;
		}
	    
	      if (isascii)
		s[bytes-1] = '\0';
	      else
		{
		  tx2bit(s, seqbuf, bytes);
		  seqbuf[length]='\0';
		}
	      titlen = strlen(seqtitle);
	      seqlen = length;
	      return length;
	    }
	  else
	    {
	      /* This should never happen, unless GCG change their file
		 format, or your .seq file does not have full INFO lines
		 (.seq files created by dataset do not, for example)
		 
		 In these cases you should create the GCG database using
		 dataset -fasta, and then flag it in FASTLIBS as a FASTA
		 format database, not GCG. */
	      
	      fprintf(stderr,
		      "Failed to understand info line:\n%s",
		      info);
	      return 0;
	    }
	}
      
    }
  return 0;
}

