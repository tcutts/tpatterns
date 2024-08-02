/*******************************************************
 *                                                     *
 *  This code is copyright (C) T J R Cutts 1998        *
 *  tjrc1@mole.bio.cam.ac.uk                           *
 *                                                     *
 *  It may be redistributed and modified, as long as   *
 *  this copyright notice is retained.                 *
 *                                                     *
 *******************************************************/

#include <string.h>
#include "tpatterns.h"

void doprompts(char *pattern, char *dbs)
{
  libent *p;
  char possibles[256];
  int m;
  int n = 0;
  
  /* Print the linked list of databases, creating
     a string of the possible letters the user can
     choose from (possibles) */

  printf("\nAvaliable databases:\n\n");
  
  for (p = firstlib;
       p != NULL;
       p = p->next)
    {
      printf("\t%c %s\n", p->letter, p->title);
      possibles[n++] = p->letter;
    }
  
  possibles[n++] = '\n';
  possibles[n] = '\0';
  
  /* Get the list of databases the user wants to search */

  for (m=0,n=-1;n!=m;)
    {
      printf("\nWhich would you like to search (e.g. AO)? ");
      
      fgets(dbs, sizeof(possibles), stdin);

      m = strlen(dbs);
      
      if (m<2)
	printf("\nNo databases given\n");
      else
	{
	  /* Did the user ask for any databases not in the list */
	  n = strspn(dbs, possibles);

	  /* If so, tell them */
	  if (n!=m)
	    printf("\nUnknown database specified: %c\n", dbs[n]);
	}
    }
  
  /* Now get the pattern to search with */
  
  for (n=0; n<2;)
    {
      printf("\nEnter the amino acid pattern to search for: ");
      
      fgets(pattern, 1024, stdin);
      n = strlen(pattern);
      /* Remove the carriage return */
      pattern[n-1] = '\0';
    }

  /* Get a filename for the output */
  
  printf("\nEnter a filename for output: ");
  fgets(possibles, sizeof(possibles), stdin);

  if ((n=strlen(possibles))<2)
    printf("\nResults will be printed to the terminal\n\n");
  else
    {
      possibles[n-1] = '\0';
      results = strdup(possibles);
      if (results==NULL)
	tp_error(TPERR_MEM);
    }

}
