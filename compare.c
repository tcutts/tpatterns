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
#include <pcre.h>
#include "tpatterns.h"

static pcre_extra *hints;

int findpattern(char *title,
				char *sequence,
				int length,
				int offset,
				pcre *query,
				int frame,
				int reversed)
{
	int result;
	char *srchfrom, c;
	int match[NMATCHES * 2];

	srchfrom = &sequence[offset];

	result = pcre_exec(query,
					   hints,
					   srchfrom,
					   length - offset,
					   0,
					   PCRE_NOTEMPTY,
					   match,
					   NMATCHES * 2);

	if (result >= 0)
	{
		c = srchfrom[match[1]];
		srchfrom[match[1]] = '\0';
		nummatches++;

		if (frame == SQ_NOFRAME)
		{
#ifdef NVSN
			if (nucleic_query)
			{
				printf("\n%s\n%sNucleotides: %d-%d\tMatch: %s\n",
					   title,
					   reversed == SQ_REVERSE ? "Reversed, " : "",
					   offset + match[0] + 1,
					   offset + match[1],
					   &srchfrom[match[0]]);
			}
			else
#endif
				printf("\n%s\nResidues: %d-%d\tMatch: %s\n",
					   title,
					   offset + match[0] + 1,
					   offset + match[1],
					   &srchfrom[match[0]]);
		}
		else
		{
			printf("\n%s\n%sReading Frame: %d\tResidues: %d-%d\tMatch: %s\n",
				   title,
				   reversed == SQ_REVERSE ? "Reverse " : "",
				   frame,
				   offset + match[0] + 1,
				   offset + match[1],
				   &srchfrom[match[0]]);
		}

		srchfrom[match[1]] = c;

		/* Start another search from after where this one started.
	   I love recursion!
		 */

		if (findall)
			findpattern(title,
						sequence,
						length,
						offset + match[0] + 1,
						query,
						frame,
						reversed);

		return match[0];
	}
	else
	{
#ifdef DEBUG
		if (result != PCRE_ERROR_NOMATCH)
			fprintf(stderr, "\n%spcre_exec returned %d\n", title, result);
#endif
		return -1;
	}
}

/* Protein ambiguity expansion */

const char Xexp[] = "[^*]";

#ifdef NVSN
/* Nucleotide ambiguity expansions */
const char Aexp[] = "[ADHMNRVWX]";
const char Bexp[] = "[CGTBU]";
const char Cexp[] = "[BCHMNSVXY]";
const char Dexp[] = "[AGTDU]";
const char Gexp[] = "[BDGKNRSVX]";
const char Hexp[] = "[ACTHU]";
const char Kexp[] = "[GTKU]";
const char Mexp[] = "[ACM]";
const char Nexp[] = "[ACGTYRNXU]";
const char Rexp[] = "[AGR]";
const char Sexp[] = "[CGS]";
const char Texp[] = "[BDHKTWXYU]";
const char Vexp[] = "[ACGV]";
const char Wexp[] = "[ATWU]";
const char Yexp[] = "[CTYU]";
#endif

#define CHAREXPAND(x, y)             \
	case x:                          \
		memcpy(q, y, sizeof(y) - 1); \
		q += sizeof(y) - 1;          \
		break

/* This function performs some mangling of the input regexp before
   compiling it.  In particular it allows for X to specify any residue */

pcre *comp_regexp(char *raw)
{
	char *munged, *p, *q;
	int in_square = 0;
	int in_brace = 0;
	int escaped = 0;
	pcre *r;
	const char *err;
	int errptr;

#ifdef DEBUG
	fprintf(stderr, "comp_regexp(): raw = %s\n", raw);
#endif

	munged = (char *)malloc(strlen(raw) << 4);

	if (!munged)
		return NULL;

	p = raw;
	q = munged;

	while (*p)
	{
		if (escaped)
			escaped = 0;
		else
		{
			switch (*p)
			{
			case '\\':
				escaped = 1;
				break;
			case '[':
				in_square = 1;
				break;
			case ']':
				in_square = 0;
				break;
			case '{':
				in_brace = 1;
				break;
			case '}':
				in_brace = 0;
				break;
			}
		}

		/* Check whether the user has escaped an ambiguity symbol */
		if (escaped)
		{
#ifdef NVSN
			if ((nucleic_query && strchr("BDHKMNRSVWX", p[1])) ||
				((!nucleic_query) && (p[1] == 'X')))
			{
				p++;
			}
#else
			if (p[1] == 'X')
			{
				p++;
			}
#endif /* NVSN else */
		}

		if (escaped || in_square || in_brace)
			*q++ = *p;
		else
		{
			/* Doing this test for every character is extremely wasteful,
			   but since this procedure is only executed once at the
			   beginning of a run it doesn't really matter too much */
#ifdef NVSN
			if (nucleic_query)
			{
				switch (*p)
				{
					CHAREXPAND('A', Aexp);
					CHAREXPAND('B', Bexp);
					CHAREXPAND('C', Cexp);
					CHAREXPAND('D', Dexp);
					CHAREXPAND('G', Gexp);
					CHAREXPAND('H', Hexp);
					CHAREXPAND('K', Kexp);
					CHAREXPAND('M', Mexp);
					CHAREXPAND('R', Rexp);
					CHAREXPAND('S', Sexp);
				case 'T':
					CHAREXPAND('U', Texp);
					CHAREXPAND('V', Vexp);
					CHAREXPAND('W', Wexp);
				case 'X':
					CHAREXPAND('N', Nexp);
				default:
					*q++ = *p;
				}
			}
			else
#endif /* NVSN */
			{
				switch (*p)
				{
					CHAREXPAND('X', Xexp);
				case '+':
				case '*':
					if (p[1] == '?')
						fprintf(stderr,
								"WARNING:  Use of %c with ? can result in very long matches.  Consider using the {m, n} construct instead.\n\n",
								*p);
				default:
					*q++ = *p;
				}
			}
		}

		p++;
	}

	*q = '\0';

#ifdef DEBUG
	fprintf(stderr, "comp_regexp(): munged = \"%s\"\n", munged);
#endif

	/* Note on PCRE_UNGREEDY:

	   PCRE_UNGREEDY causes + and * to make the smallest possible match.
	   Note that this breaks exact compatibility with Perl regular
	   expressions, but is more useful behaviour for our purposes. */

	r = pcre_compile(munged, PCRE_UNGREEDY, &err, &errptr, NULL);

	if (r == NULL)
	{
		fprintf(stderr,
				"tpatterns: error in pcre at offset %d: %s\n",
				errptr, err);
	}
	else
	{

		/* Using pcre_study can make certain patterns twice as fast; namely
	   those beginning with an ambiguity */

		hints = pcre_study(r, 0, &err);

		if (err != NULL)
		{
			fprintf(stderr,
					"tpatterns: error while studing pcre: %s\n",
					err);
			return NULL;
		}
	}

	free(munged);

	return r;
}
