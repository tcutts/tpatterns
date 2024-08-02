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

char buf[1024];
char tmptitle[1024];

void init_read_fasta(void)
{
	tmptitle[0] = '\0';
}

int read_fasta(FILE *f)
{
	int length, n;
	char *p;

	length = 0;
	seqbuf[0] = '\0';

	while (!feof(f))
	{
		p = fgets(buf, sizeof(buf), f);
		if (p)
		{
			if (*p == '>')
			{
				/* We have found a comment line */
				p[strlen(p) - 1] = '\0';

				if (tmptitle[0] == '\0')
				{
					/* This is our first sequence, so carry on reading */
					memcpy(tmptitle, buf, sizeof(buf));
				}
				else
				{
					memcpy(seqtitle, tmptitle, sizeof(tmptitle));
					len[0] = strlen(seqtitle);
					memcpy(tmptitle, buf, sizeof(buf));
					if (len[0] != 0)
					{
						seqbuf[length] = '\0';
						len[1] = length;
						return length;
					}
					return 0;
				}
			}
			else
			{
				/* Data line */
				n = strlen(p) - 1;
				if ((length + n + 3) > MAXSEQLEN)
				{
					tp_error(TPERR_SEQLEN);
					return 0;
				}
				memcpy(&seqbuf[length], buf, n);
				length += n;
			}
		}
	}
	if (length != 0)
	{
		seqbuf[length] = '\0';
		len[1] = length;
		memcpy(seqtitle, tmptitle, sizeof(seqtitle));
		len[0] = strlen(seqtitle);
	}
	return length;
}
