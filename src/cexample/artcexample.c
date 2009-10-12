/*
 * =====================================================================================
 *
 *       Filename:  artcexample.c
 *
 *    Description:  This program is an example code to demonstrate how to use
 *    		    the Artimagen library in C.
 *
 *        Version:  1.0
 *        Created:  09/29/2009 12:41:55 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dr. Petr Cizmar (pc), petr.cizmar@nist.gov
 *        Company:  National Institute od Standards and Technology
 *
 * =====================================================================================
 */

/* 
    Artificial SEM Image Generator (ArtImaGen)
    2009  Petr Cizmar @ National Institute of Standards and Technology
    E-mail: petr.cizmar@nist.gov

    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to Gnu.org public domain is compatible with GPL.

 */

/*
 * COMPILATION NOTE:
 *
 * If this program may be compiled with a gcc (C compiler) only when linking
 * against shared version of artimagen library. Static linking is possible only
 * using the g++ (C++ compiler), although, this program users plain C only.
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "../libartimagen/artimagen.h"
#include "../../config.h"


/* missing functions */
#ifndef HAVE_SRANDOM
 #ifdef WIN32
  #define srandom srand
 #else
  #define srandom srand48
 #endif
#endif

#ifndef HAVE_RANDOM
 #ifdef WIN32
  #define random (long)rand
 #else
  #define random (long)lrand48
 #endif
#endif

int main(){
   printf ("\n*********\nArtimagen in C example.\nby Petr Cizmar @ NIST, September 2009\n*********\n\n");

   DIST_TYPE size = 512;

   t_gc_definition gc_def; 
   t_std_image_def_struct im_def;

   srandom(time(0));

// define the sample
   gc_def.sizex = size;
   gc_def.sizey = size;
   gc_def.ee_coefficient = 0.2;
   gc_def.ee_top_above_basic = 0.3;
   gc_def.basic_level = 0.5;
   gc_def.basic_level_variation = 0.1;
   gc_def.grain_min_size = 5;
   gc_def.grain_max_size = 15;
   gc_def.number_of_grains = size*size/2000;
   gc_def.rotation = 0;
   gc_def.fs_density = 7e-3;
   gc_def.fs_min_r = 6;
   gc_def.fs_max_r = 10;
   gc_def.fs_min_coe = 0.95;
   gc_def.fs_max_coe = 1.05;

//define the image
   im_def.sizex = size;
   im_def.sizey = size;
   im_def.bg_min_gl = 0.1;
   im_def.bg_max_gl = 0.3;
   im_def.bg_dens_x = 5;
   im_def.bg_dens_y = 5;
   im_def.beam_sigma = 1;
   im_def.beam_astig_ratio = 1;
   im_def.beam_astig_angle = 30;
   im_def.shift_x = 0;
   im_def.shift_y = 0;
   im_def.vib_pixel_dwell_time = 10000;
   im_def.vib_min_frequency = 3;
   im_def.vib_max_frequency = 200;
   im_def.vib_max_amplitude = 0.2;
   im_def.vib_number_of_frequencies = 10;
   im_def.vib_pixel_dead_time = 100;
   im_def.vib_line_dead_time = 100000;
   im_def.noise_sigma = 0.08;

   void *sample = generate_gc_sample(&gc_def);
   printf("sample OK\n");
   void *image = generate_standard_image(sample, &im_def);
   printf("image OK\n");

   save_image(image, (char *)"cex.tiff", (char *)"generated using Artimagen library by Petr Cizmar @ NIST.");

   destroy_gc_sample(sample);
   destroy_image(image);
   return 0;
}

//vim:set cindent
