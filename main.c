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
#include <errno.h>
#include <stdlib.h>
#include "tpatterns.h"

int numseqs = 0;
int nummatches = 0;
int findall = 1;
char *results = NULL;
char pattern[1024];
char srchlibs[256];

#ifdef NVSN
int nucleic_query = 0;
#endif

int main(int argc, char *argv[])
{
	pcre *r;
	int offset;
	int libfile = -1;
	int n;
	libent *lib;

	mpi_error = MPI_Init(&argc, &argv);

	if (mpi_error != MPI_SUCCESS)
	{
		printf("MPI initialisation error\n");
		return 1;
	}

	MPI_Comm_size(MPI_COMM_WORLD, &numworkers);
	MPI_Comm_rank(MPI_COMM_WORLD, &thisworker);

	offset = options(argc, argv, &libfile);

	/* Build linked list of available libraries */
	if ((offset > 0) && (build_libents(argv, libfile) == 0))
	{
		if (thisworker == 0)
		{

			if ((argc - offset) == 1)
			{
				usage(argv[0], "");
				MPI_Finalize();
				return 1;
			}

			printf("%s %s (c) T J R Cutts 1998\n\n",
				   TNAME, TFVER);

			if (argc == offset)
			{
				usage(argv[0], "MPI version does not support interactive use");
				MPI_Finalize();
				return 1;
			}
			else
				strcpy(srchlibs, argv[offset + 1]);
		}
		else
		{
			/* Worker-specific initialisation */

			initbasebits();
			make_revmatrix();

			strcpy(pattern, argv[offset]);

			if ((r = comp_regexp(pattern)) == NULL)
			{
				if (thisworker == 1)
					fprintf(stderr, "Could not understand your pattern\n");

				MPI_Finalize();
				return 1;
			}
		}

		/* Now the main work of the instances */

		if (thisworker == 0)
		{
			for (lib = firstlib; lib != NULL; lib = lib->next)
			{
				if (strchr(srchlibs, lib->letter) != NULL)
				{
					searchdb(lib);

					fprintf(stderr, "Finished searching %c\n",
							lib->letter);
					fflush(stderr);
				}
			}

			/* Send a tiny message to tell all the workers
			   the data is all gone */

			for (n = 1; n < numworkers; n++)
				MPI_Send(seqbuf,
						 1,
						 MPI_CHAR,
						 n,
						 0,
						 MPI_COMM_WORLD);
		}
		else
		{
#ifdef DEBUG
			fprintf(stderr, "I am worker %d, with pid %d\n", thisworker, getpid());
			sleep(20);
#endif
			mpi_receive_data(r);
			pcre_free(r);
		}
	}

	MPI_Finalize();

	freelibents(firstlib);

	return 0;
}
