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
#include <string.h>
#include "tpatterns.h"

#ifdef NVSN
#define NUCOPT "[-n] "
#else
#define NUCOPT ""
#endif

/* usage() prints a usage message and returns a non-zero
   exit status */

void usage(char *a, char *b)
{
	fprintf(stderr,
			"%s\n\nUsage:\n\t%s %s[-1] [-l libfile] [-m translation] [-o outfile] [[pattern] [database]]\n",
			b, a, NUCOPT);
}

/* options() processes the command line arguments,
   returning the index of the first argument that
   is not an option */

int options(int argc, char *argv[], int *libfile)
{
#ifdef NVSN
	char *opts = "1hl:m:no:";
#else
	char *opts = "1hl:m:o:";
#endif
	extern char *optarg;
	extern int optind;
	int c;

	while ((c = getopt(argc, argv, opts)) != -1)
	{
		switch (c)
		{
		case '1':
			findall = 0;
			break;
#ifdef NVSN
		case 'n':
			nucleic_query++;
			break;
#endif
		case 'h':
			usage(argv[0], "");
			break;
		case 'l':
			libfile[0] = optind - 1;
			break;
		case 'm':
			if (compilemx(argv[optind - 1]))
				return -1;
			break;
		case 'o':
		{
			if ((results = strdup(argv[optind - 1])) == NULL)
				tp_error(TPERR_MEM);
			break;
		}
		default:
			usage(argv[0], "Unknown option");
		}
	}

	return optind;
}
