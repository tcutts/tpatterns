#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include "tplib.h"

/* This trivial program saves a FASTA format database from STDIN
   as a binary file with the lengths of the title and sequences binary
   encoded at the beginning of each record, so that programs can use
   read(2) to load the data from disk very much more quickly */

int main(int argc, char *argv[])
{
  int tlen, slen;
  int maxslen = 0;

  init_read_fasta();
  while ((slen = read_fasta(stdin)) > 0)
  {
    if (slen > maxslen)
      maxslen = slen;

    tlen = strlen(seqtitle);

    if (tlen > 257)
      tlen = 257;

    memset(&seqtitle[tlen], '\0', sizeof(char) * (258 - tlen));

    tlen--;

    write(STDOUT_FILENO, &tlen, sizeof(int));
    write(STDOUT_FILENO, &slen, sizeof(int));
    write(STDOUT_FILENO, &seqtitle[1], tlen);
    write(STDOUT_FILENO, &seqbuf, slen);
  }

  fprintf(stderr, "Maximum sequence length: %d\n", maxslen);

  return 0;
}
