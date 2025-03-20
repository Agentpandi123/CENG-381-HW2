/*"Canny" edge detector code:
---------------------------

This text file contains the source code for a "Canny" edge detector. It
was written by Mike Heath (heath@csee.usf.edu) using some pieces of a
Canny edge detector originally written by someone at Michigan State
University.

There are three 'C' source code files in this text file. They are named
"canny_edge.c", "hysteresis.c" and "pgm_io.c". They were written and compiled
under SunOS 4.1.3. Since then they have also been compiled under Solaris.
To make an executable program: (1) Separate this file into three files with
the previously specified names, and then (2) compile the code using

  gcc -o canny canny.c hysteresis.c pgm_io.c -lm
  (Note: You can also use optimization such as -O3)

The resulting program, canny_edge, will process images in the PGM format.
Parameter selection is left up to the user. A broad range of parameters to
use as a starting point are: sigma 0.60-2.40, tlow 0.20-0.50 and,
thigh 0.60-0.90.

If you are using a Unix system, PGM file format conversion tools can be found
at ftp://wuarchive.wustl.edu/graphics/graphics/packages/pbmplus/.
Otherwise, it would be easy for anyone to rewrite the image I/O procedures
because they are listed in the separate file pgm_io.c.

If you want to check your compiled code, you can download grey-scale and edge
images from http://marathon.csee.usf.edu/edge/edge_detection.html. You can use
the parameters given in the edge filenames and check whether the edges that
are output from your program match the edge images posted at that address.

Mike Heath
(10/29/96)*/

/*------------------------- begin canny_edge.c -------------------------*/
/*******************************************************************************
* --------------------------------------------
*(c) 2001 University of South Florida, Tampa
* Use, or copying without permission prohibited.
* PERMISSION TO USE
* In transmitting this software, permission to use for research and
* educational purposes is hereby granted.  This software may be copied for
* archival and backup purposes only.  This software may not be transmitted
* to a third party without prior permission of the copyright holder. This
* permission may be granted only by Mike Heath or Prof. Sudeep Sarkar of
* University of South Florida (sarkar@csee.usf.edu). Acknowledgment as
* appropriate is respectfully requested.
* 
*  Heath, M., Sarkar, S., Sanocki, T., and Bowyer, K. Comparison of edge
*    detectors: a methodology and initial study, Computer Vision and Image
*    Understanding 69 (1), 38-54, January 1998.
*  Heath, M., Sarkar, S., Sanocki, T. and Bowyer, K.W. A Robust Visual
*    Method for Assessing the Relative Performance of Edge Detection
*    Algorithms, IEEE Transactions on Pattern Analysis and Machine
*    Intelligence 19 (12),  1338-1359, December 1997.
*  ------------------------------------------------------
*
* PROGRAM: canny_edge
* PURPOSE: This program implements a "Canny" edge detector. The processing
* steps are as follows:
*
*   1) Convolve the image with a separable gaussian filter.
*   2) Take the dx and dy the first derivatives using [-1,0,1] and [1,0,-1]'.
*   3) Compute the magnitude: sqrt(dx*dx+dy*dy).
*   4) Perform non-maximal suppression.
*   5) Perform hysteresis.
*
* The user must input three parameters. These are as follows:
*
*   sigma = The standard deviation of the gaussian smoothing filter.
*   tlow  = Specifies the low value to use in hysteresis. This is a 
*           fraction (0-1) of the computed high threshold edge strength value.
*   thigh = Specifies the high value to use in hysteresis. This fraction (0-1)
*           specifies the percentage point in a histogram of the gradient of
*           the magnitude. Magnitude values of zero are not counted in the
*           histogram.
*
* NAME: Mike Heath
*       Computer Vision Laboratory
*       University of South Floeida
*       heath@csee.usf.edu
*
* DATE: 2/15/96
*
* Modified: 5/17/96 - To write out a floating point RAW headerless file of
*                     the edge gradient "up the edge" where the angle is
*                     defined in radians counterclockwise from the x direction.
*                     (Mike Heath)
*******************************************************************************/
/*******************************************************************************
* PROGRAM: canny_edge
* PURPOSE: This program implements a "Canny" edge detector. 
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define VERBOSE 0
#define BOOSTBLURFACTOR 90.0

int read_pgm_image(char *infilename, unsigned char **image, int *rows, int *cols);
int write_pgm_image(char *outfilename, unsigned char *image, int rows, int cols, char *comment, int maxval);

void canny(unsigned char *image, int rows, int cols, float sigma,
           float tlow, float thigh, unsigned char **edge, char *fname);
void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma, short int **smoothedim);
void make_gaussian_kernel(float sigma, float **kernel, int *windowsize);
void derivative_x_y(short int *smoothedim, int rows, int cols, short int **delta_x, short int **delta_y);
void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols, short int **magnitude);
void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols, float tlow, float thigh, unsigned char *edge);
void radian_direction(short int *delta_x, short int *delta_y, int rows, int cols, float **dir_radians, int xdirtag, int ydirtag);
double angle_radians(double x, double y);
void non_max_supp(short int *mag, short int *gradx, short int *grady, int rows, int cols, unsigned char *result);

int main(int argc, char *argv[]) {
    char *infilename = NULL;
    char *dirfilename = NULL;
    char outfilename[128];
    char composedfname[128];
    unsigned char *image;
    unsigned char *edge;
    int rows, cols;
    float sigma, tlow, thigh;

    if (argc < 5) {
        fprintf(stderr, "\n<USAGE> %s image sigma tlow thigh [writedirim]\n", argv[0]);
        exit(1);
    }

    infilename = argv[1];
    sigma = atof(argv[2]);
    tlow = atof(argv[3]);
    thigh = atof(argv[4]);

    if (argc == 6) dirfilename = infilename;
    else dirfilename = NULL;

    if (read_pgm_image(infilename, &image, &rows, &cols) == 0) {
        fprintf(stderr, "Error reading the input image, %s.\n", infilename);
        exit(1);
    }

    if (dirfilename != NULL) {
        sprintf(composedfname, "%s_s_%3.2f_l_%3.2f_h_%3.2f.fim", infilename, sigma, tlow, thigh);
        dirfilename = composedfname;
    }
    canny(image, rows, cols, sigma, tlow, thigh, &edge, dirfilename);

    sprintf(outfilename, "%s_s_%3.2f_l_%3.2f_h_%3.2f.pgm", infilename, sigma, tlow, thigh);
    if (write_pgm_image(outfilename, edge, rows, cols, "", 255) == 0) {
        fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
        exit(1);
    }

    free(image);
    free(edge);
    return 0;
}

void canny(unsigned char *image, int rows, int cols, float sigma,
           float tlow, float thigh, unsigned char **edge, char *fname) {
    FILE *fpdir = NULL;
    unsigned char *nms;
    short int *smoothedim, *delta_x, *delta_y, *magnitude;
    float *dir_radians = NULL;

    gaussian_smooth(image, rows, cols, sigma, &smoothedim);

    derivative_x_y(smoothedim, rows, cols, &delta_x, &delta_y);

    if (fname != NULL) {
        radian_direction(delta_x, delta_y, rows, cols, &dir_radians, -1, -1);
        if ((fpdir = fopen(fname, "wb")) == NULL) {
            fprintf(stderr, "Error opening the file %s for writing.\n", fname);
            exit(1);
        }
        fwrite(dir_radians, sizeof(float), rows * cols, fpdir);
        fclose(fpdir);
        free(dir_radians);
    }

    magnitude_x_y(delta_x, delta_y, rows, cols, &magnitude);

    if ((nms = (unsigned char *) calloc(rows * cols, sizeof(unsigned char))) == NULL) {
        fprintf(stderr, "Error allocating the nms image.\n");
        exit(1);
    }
    non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);

    if ((*edge = (unsigned char *) calloc(rows * cols, sizeof(unsigned char))) == NULL) {
        fprintf(stderr, "Error allocating the edge image.\n");
        exit(1);
    }
    apply_hysteresis(magnitude, nms, rows, cols, tlow, thigh, *edge);

    free(smoothedim);
    free(delta_x);
    free(delta_y);
    free(magnitude);
    free(nms);
}

void derivative_x_y(short int *smoothedim, int rows, int cols, short int **delta_x, short int **delta_y) {
    int r, c, pos;

    if (((*delta_x) = (short int *) calloc(rows * cols, sizeof(short int))) == NULL) {
        fprintf(stderr, "Error allocating the delta_x image.\n");
        exit(1);
    }
    if (((*delta_y) = (short int *) calloc(rows * cols, sizeof(short int))) == NULL) {
        fprintf(stderr, "Error allocating the delta_y image.\n");
        exit(1);
    }

    for (r = 0; r < rows; r++) {
        pos = r * cols;
        (*delta_x)[pos] = smoothedim[pos + 1] - smoothedim[pos];
        pos++;
        for (c = 1; c < cols - 1; c++, pos++) {
            (*delta_x)[pos] = smoothedim[pos + 1] - smoothedim[pos - 1];
        }
        (*delta_x)[pos] = smoothedim[pos] - smoothedim[pos - 1];
    }

    for (c = 0; c < cols; c++) {
        pos = c;
        (*delta_y)[pos] = smoothedim[pos + cols] - smoothedim[pos];
        pos += cols;
        for (r = 1; r < rows - 1; r++, pos += cols) {
            (*delta_y)[pos] = smoothedim[pos + cols] - smoothedim[pos - cols];
        }
        (*delta_y)[pos] = smoothedim[pos] - smoothedim[pos - cols];
    }
}

void radian_direction(short int *delta_x, short int *delta_y, int rows, int cols, float **dir_radians, int xdirtag, int ydirtag) {
    int r, c, pos;
    float dx, dy;

    if (((*dir_radians) = (float *) calloc(rows * cols, sizeof(float))) == NULL) {
        fprintf(stderr, "Error allocating the direction image.\n");
        exit(1);
    }

    for (r = 0; r < rows; r++) {
        for (c = 0; c < cols; c++) {
            pos = r * cols + c;
            dx = (float)delta_x[pos];
            dy = (float)delta_y[pos];
            (*dir_radians)[pos] = atan2(dy, dx);
        }
    }
}

void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols, short int **magnitude) {
    int r, c, pos;

    if (((*magnitude) = (short int *) calloc(rows * cols, sizeof(short int))) == NULL) {
        fprintf(stderr, "Error allocating the magnitude image.\n");
        exit(1);
    }

    for (r = 0; r < rows; r++) {
        for (c = 0; c < cols; c++) {
            pos = r * cols + c;
            (*magnitude)[pos] = (short int)(sqrt((double)(delta_x[pos] * delta_x[pos] + delta_y[pos] * delta_y[pos])));
        }
    }
}

void non_max_supp(short int *mag, short int *gradx, short int *grady, int rows, int cols, unsigned char *result) {
    int r, c, pos;
    float angle, mag1, mag2;
    const float PI = 3.14159265358979323846;

    for (r = 1; r < rows - 1; r++) {
        for (c = 1; c < cols - 1; c++) {
            pos = r * cols + c;

            if (gradx[pos] == 0) {
                angle = 90.0;
            } else {
                angle = atan2(grady[pos], gradx[pos]) * 180.0 / PI;
            }

            if ((angle >= -22.5 && angle < 22.5) || (angle >= 157.5 || angle < -157.5)) {
                mag1 = mag[pos - 1];
                mag2 = mag[pos + 1];
            } else if ((angle >= 22.5 && angle < 67.5) || (angle >= -157.5 && angle < -112.5)) {
                mag1 = mag[pos - cols + 1];
                mag2 = mag[pos + cols - 1];
            } else if ((angle >= 67.5 && angle < 112.5) || (angle >= -112.5 && angle < -67.5)) {
                mag1 = mag[pos - cols];
                mag2 = mag[pos + cols];
            } else {
                mag1 = mag[pos - cols - 1];
                mag2 = mag[pos + cols + 1];
            }

            if (mag[pos] >= mag1 && mag[pos] >= mag2) {
                result[pos] = (unsigned char)mag[pos];
            } else {
                result[pos] = 0;
            }
        }
    }
}

void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma, short int **smoothedim) {
    int r, c, rr, cc, windowsize, center;
    float *tempim, *kernel, dot, sum;

    make_gaussian_kernel(sigma, &kernel, &windowsize);
    center = windowsize / 2;

    if ((tempim = (float *) calloc(rows * cols, sizeof(float))) == NULL) {
        fprintf(stderr, "Error allocating the buffer image.\n");
        exit(1);
    }
    if (((*smoothedim) = (short int *) calloc(rows * cols, sizeof(short int))) == NULL) {
        fprintf(stderr, "Error allocating the smoothed image.\n");
        exit(1);
    }

    for (r = 0; r < rows; r++) {
        for (c = 0; c < cols; c++) {
            dot = 0.0;
            sum = 0.0;
            for (cc = -center; cc <= center; cc++) {
                if ((c + cc) >= 0 && (c + cc) < cols) {
                    dot += (float)image[r * cols + (c + cc)] * kernel[center + cc];
                    sum += kernel[center + cc];
                }
            }
            tempim[r * cols + c] = dot / sum;
        }
    }

    for (c = 0; c < cols; c++) {
        for (r = 0; r < rows; r++) {
            dot = 0.0;
            sum = 0.0;
            for (rr = -center; rr <= center; rr++) {
                if ((r + rr) >= 0 && (r + rr) < rows) {
                    dot += tempim[(r + rr) * cols + c] * kernel[center + rr];
                    sum += kernel[center + rr];
                }
            }
            (*smoothedim)[r * cols + c] = (short int)(dot * BOOSTBLURFACTOR / sum + 0.5);
        }
    }

    free(tempim);
    free(kernel);
}

void make_gaussian_kernel(float sigma, float **kernel, int *windowsize) {
    int i, center;
    float x, fx, sum = 0.0;

    *windowsize = 1 + 2 * ceil(2.5 * sigma);
    center = (*windowsize) / 2;

    if ((*kernel = (float *) calloc((*windowsize), sizeof(float))) == NULL) {
        fprintf(stderr, "Error allocating the Gaussian kernel array.\n");
        exit(1);
    }

    for (i = 0; i < (*windowsize); i++) {
        x = (float)(i - center);
        fx = exp(-0.5 * x * x / (sigma * sigma)) / (sigma * sqrt(6.2831853));
        (*kernel)[i] = fx;
        sum += fx;
    }

    for (i = 0; i < (*windowsize); i++) (*kernel)[i] /= sum;
}


// Implement the rest of the functions (canny, gaussian_smooth, etc.)...

/*------------------------- end canny_edge.c -------------------------*/
