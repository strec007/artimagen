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

   fill_gc_sample_definition_structure(
	 &gc_def,                  //t_gc_definition *def_str,
	 size,                  //int sizex,
	 size,                  //int sizey,
	 0.2,                //float ee_coefficient,
	 0.3,                  //IM_STORE_TYPE ee_top_above_basic,
	 0,                  //reserved,
	 0.5,                  //IM_STORE_TYPE basic_level,
	 0.1,                  //IM_STORE_TYPE basic_level_variation,
	 5,                  //DIST_TYPE grain_min_size,
	 15,                  //DIST_TYPE grain_max_size,
	 size*size/2000,                  //int number_of_grains,
	 0,                  //float rotation,
	 7e-3,                  //float fs_density,
	 6,                  //DIST_TYPE fs_min_r,
	 10,                  //DIST_TYPE fs_max_r,
	 0.95,                  //float fs_min_coe,
	 1.05                  //float fs_max_coe
	 );
   
   fill_image_definition_structure(
	 &im_def,	//      t_std_image_def_struct *def_str,
	 size,		//      DIST_TYPE sizex,
	 size,		//      DIST_TYPE sizey,
	 0.1,		//      IM_STORE_TYPE bg_min_gl,
	 0.3,		//      IM_STORE_TYPE bg_max_gl,
	 5,		//      int bg_dens_x,
	 5,		//      int bg_dens_y,
	 1,		//      float beam_sigma,
	 1,		//      float beam_astig_ratio,
	 30,		//      float beam_astig_angle,
	 0,		//      float shift_x,
	 0,		//      float shift_y,
	 10000,		//      unsigned int vib_pixel_dwell_time,
	 3,		//      float vib_min_frequency,
	 200,		//      float vib_max_frequency,
	 0.2,		//      float vib_max_amplitude,
	 10,		//      int vib_number_of_frequencies,
	 100,		//      unsigned int vib_pixel_dead_time,
	 100000,	//      unsigned int vib_line_dead_time,
	 0.08		//      double noise_sigma);
	 );

   void *sample = generate_gc_sample(&gc_def);
   printf("sample OK\n");
   void *image = generate_standard_image(sample, &im_def);
   printf("image OK\n");

   save_image(image, "cex.tiff", "generated using Artimagen library by Petr Cizmar @ NIST.");

   destroy_gc_sample(sample);
   destroy_image(image);
   return 0;
}

//vim:set cindent
