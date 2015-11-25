/* -------------------------------------------------------------------------------- */
/*                                    ascii2dump                                    */
/* -------------------------------------------------------------------------------- */
/*   (c) University of Munich, Database Group                                       */
/*   written by M. Breunig, breunig@dbs.informatik.uni-muenchen.de                  */
/* -------------------------------------------------------------------------------- */
/* Converts an ascii file into a binary file.                                       */
/* -------------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>

void main(int argc, char **argv) {
  int a=0, dim;
  char *infile, *outfile;
  float fval;
  FILE *afp, *bfp;
  
  printf("\nConvert ASCII-data into binary-data\n");
  if(argc!=4) { 
    printf("  ** Usage: %s <ascii-file> <dimension> <binary-file>\n\n",argv[0]); exit(-2); 
  }
  
  infile = argv[1];
  dim = atoi(argv[2]);
  outfile = argv[3];
  
  afp = fopen(infile, "rb");
  if (!afp) { perror("cannot open input file\n"); exit(-1); }
  
  bfp = fopen(outfile, "wb");
  if (!bfp) { perror("cannot open output file\n"); exit(-1); }
  
  printf("Input file opened - converting...\n");
  
  while ( fscanf(afp, "%f", &fval) > 0 ) {
    fwrite(&fval, sizeof(float), 1, bfp);
    a++;
	  if ((a%(1000*dim)) == 0) {printf("\n%d ", (int) (a/dim)); fflush(stdout); }
  }
  printf("\r%d   \n\n", (int) (a/dim)); fflush(stdout);

  fclose(afp);
  fclose(bfp);
}
