/*------------------------- begin pgm_io.c-------------------------*/
/*******************************************************************************
* FILE: pgm_io.c
* This code was written by Mike Heath. heath@csee.usf.edu (in 1995).
*******************************************************************************/
/*******************************************************************************
* FILE: pgm_io.c
* This code was written by Mike Heath. (in 1995).
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/******************************************************************************
* Function: read_pgm_image
* Purpose: Reads an image in PGM format.
******************************************************************************/
int read_pgm_image(char *infilename, unsigned char **image, int *rows, int *cols) {
   FILE *fp;
   char buf[71];

   if (infilename == NULL) fp = stdin;
   else {
      if ((fp = fopen(infilename, "r")) == NULL) {
         fprintf(stderr, "Error reading the file %s in read_pgm_image().\n", infilename);
         return 0;
      }
   }

   fgets(buf, 70, fp);
   if (strncmp(buf, "P5", 2) != 0) {
      fprintf(stderr, "The file %s is not in PGM format.\n", infilename);
      if (fp != stdin) fclose(fp);
      return 0;
   }
   do { fgets(buf, 70, fp); } while (buf[0] == '#');
   sscanf(buf, "%d %d", cols, rows);
   do { fgets(buf, 70, fp); } while (buf[0] == '#');

   if (((*image) = (unsigned char *) malloc((*rows) * (*cols))) == NULL) {
      fprintf(stderr, "Memory allocation failure in read_pgm_image().\n");
      if (fp != stdin) fclose(fp);
      return 0;
   }
   if ((*rows) != fread((*image), (*cols), (*rows), fp)) {
      fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
      if (fp != stdin) fclose(fp);
      free((*image));
      return 0;
   }

   if (fp != stdin) fclose(fp);
   return 1;
}

/******************************************************************************
* Function: write_pgm_image
* Purpose: Writes an image in PGM format.
******************************************************************************/
int write_pgm_image(char *outfilename, unsigned char *image, int rows, int cols, char *comment, int maxval) {
   FILE *fp;

   if (outfilename == NULL) fp = stdout;
   else {
      if ((fp = fopen(outfilename, "w")) == NULL) {
         fprintf(stderr, "Error writing the file %s in write_pgm_image().\n", outfilename);
         return 0;
      }
   }

   fprintf(fp, "P5\n%d %d\n", cols, rows);
   if (comment != NULL)
      if (strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
   fprintf(fp, "%d\n", maxval);

   if (rows != fwrite(image, cols, rows, fp)) {
      fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
      if (fp != stdout) fclose(fp);
      return 0;
   }

   if (fp != stdout) fclose(fp);
   return 1;
}

/*------------------------- end pgm_io.c -------------------------*/
