/*
 * =====================================================================================
 *
 *       Filename:  wrapper.cxx
 *
 *    Description:  This file is a part of the artimagen library. It serves as a
 *    		    wrapper, so some basic artimagen methods could be called
 *    		    from C programs.  The wrapper functions are available to C++
 *    		    applications.
 *
 *        Version:  1.0
 *        Created:  09/28/2009 01:31:00 PM
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
 

#include "artimagen.h"
#include "artimagen_i.h"
#include "../../config.h"


extern "C" void *fill_image_definition_structure(
      t_std_image_def_struct *def_str,
      IM_STORE_TYPE bg_min_gl,
      IM_STORE_TYPE bg_max_gl,
      int bg_dens_x,
      int bg_dens_y,
      float beam_sigma,
      float beam_astig_ratio,
      float beam_astig_angle,
      float shift_x,
      float shift_y,
      unsigned int vib_pixel_dwell_time,
      float vib_min_frequency,
      float vib_max_frequency,
      float vib_max_amplitude,
      int vib_number_of_frequencies,
      unsigned int vib_pixel_dead_time,
      unsigned int vib_line_dead_time,
      double noise_sigma){
   def_str->bg_min_gl = bg_min_gl;
   def_str->bg_max_gl = bg_max_gl;
   def_str->bg_dens_x = bg_dens_x;
   def_str->bg_dens_y = bg_dens_y;
   def_str->beam_sigma = beam_sigma;
   def_str->beam_astig_ratio = beam_astig_ratio;
   def_str->beam_astig_angle = beam_astig_angle;
   def_str->shift_x = shift_x;
   def_str->shift_y = shift_y;
   def_str->vib_pixel_dwell_time = vib_pixel_dwell_time;
   def_str->vib_min_frequency = vib_min_frequency;
   def_str->vib_max_frequency = vib_max_frequency;
   def_str->vib_max_amplitude = vib_max_amplitude;
   def_str->vib_number_of_frequencies = vib_number_of_frequencies;
   def_str->vib_pixel_dead_time = vib_pixel_dead_time;
   def_str->vib_line_dead_time = vib_line_dead_time;
   def_str->noise_sigma = noise_sigma;
}

extern "C" void *fill_gc_sample_definition_structure(
      t_gc_definition *def_str,
      int sizex,
      int sizey,
      float ee_coefficient,
      IM_STORE_TYPE ee_top_above_basic,
      DIST_TYPE ee_thickness,
      IM_STORE_TYPE basic_level,
      IM_STORE_TYPE basic_level_variation,
      DIST_TYPE grain_min_size,
      DIST_TYPE grain_max_size,
      int number_of_grains,
      float rotation,
      float fs_density,
      DIST_TYPE fs_min_r,
      DIST_TYPE fs_max_r,
      float fs_min_coe,
      float fs_max_coe
){
   def_str->sizex = sizex;
   def_str->sizey = sizey;
   def_str->ee_coefficient = ee_coefficient;
   def_str->ee_top_above_basic = ee_top_above_basic;
   def_str->ee_thickness = ee_thickness;
   def_str->basic_level = basic_level;
   def_str->basic_level_variation = basic_level_variation;
   def_str->grain_min_size = grain_min_size;
   def_str->grain_max_size = grain_max_size;
   def_str->number_of_grains = number_of_grains;
   def_str->rotation = rotation;
   def_str->fs_density = fs_density;
   def_str->fs_min_r = fs_min_r;
   def_str->fs_max_r = fs_max_r;
   def_str->fs_min_coe = fs_min_coe;
   def_str->fs_max_coe = fs_max_coe;
}
 // vim: cindent
