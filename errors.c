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
#include "tplib.h"

const char *messages[] = {
    "Success",                      /* TF_SUCCESS */
    "Memory allocation failure",    /* TF_MEM */
    "Premature end of matrix file", /* TF_MXEOF */
    "Could not open file for output",
    "Sequence length exceeds MAXSEQLEN"};

void tp_error(tperror code)
{
  fputs(messages[code], stderr);
}
