/*
 * =====================================================================================
 *
 *       Filename:  artimagen.h
 *
 *    Description:  This is a header file to be used by the C, C++, or other
 *    		    applications.
 *
 *        Version:  1.0
 *        Created:  09/28/2009 01:39:00 PM
 *       Revision:  none
 *       Compiler:  gcc / g++
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
 * IMPORTANT NOTE: include ../../config.h before this file, or define HAVE_LUA,
 * if you are compiling with LUA support.
 */

#ifndef ARTIMAGEN_H
#define ARTIMAGEN_H

#ifdef HAVE_LUA
extern "C" {
#include <lua5.1/lualib.h>
#include <lua5.1/lauxlib.h>
}
#endif 

// definition of distance types
#define DIST_TYPE double
#define GE_INDEX_TYPE unsigned long int

//definition of image bit depth flags
#define BI_8B 0
#define BI_16B 1

//definition of image maximal values and types
#define TYPE_8B uint8_t
#define TYPE_16B uint16_t

// definition of storage type and type for image coordinates
#define IM_COORD_TYPE int
#define IM_STORE_TYPE double


// t_std_image_def_struct 
typedef struct{/*{{{*/
   DIST_TYPE sizex;
   DIST_TYPE sizey;
   IM_STORE_TYPE bg_min_gl;
   IM_STORE_TYPE bg_max_gl;
   int bg_dens_x;
   int bg_dens_y;
   float beam_sigma;
   float beam_astig_ratio;
   float beam_astig_angle;
   float shift_x;
   float shift_y;
   unsigned int vib_pixel_dwell_time;
   float vib_min_frequency;
   float vib_max_frequency;
   float vib_max_amplitude;
   int vib_number_of_frequencies;
   unsigned int vib_pixel_dead_time;
   unsigned int vib_line_dead_time;
   double noise_sigma;
} t_std_image_def_struct;/*}}}*/

// t_gc_definition structure
typedef struct {/*{{{*/
   int sizex;
   int sizey;
   float ee_coefficient;
   IM_STORE_TYPE ee_top_above_base;
   DIST_TYPE ee_thickness;
   IM_STORE_TYPE base_level;
   IM_STORE_TYPE base_level_variation;
   DIST_TYPE grain_min_size;
   DIST_TYPE grain_max_size;
   int number_of_grains;
   float rotation;
   float fs_density;
   DIST_TYPE fs_min_r;
   DIST_TYPE fs_max_r;
   float fs_min_coe;
   float fs_max_coe;
} t_gc_definition;/*}}}*/

// t_cor_definition structure
typedef struct {/*{{{*/
   int sizex;
   int sizey;
   float ee_coefficient;
   IM_STORE_TYPE ee_top_above_base;
   DIST_TYPE ee_thickness;
   IM_STORE_TYPE base_level;
   IM_STORE_TYPE base_level_variation;
   DIST_TYPE lsize;
   DIST_TYPE tsize;
   DIST_TYPE distance; 
   double rotation;
   float fs_density;
   DIST_TYPE fs_min_r;
   DIST_TYPE fs_max_r;
   float fs_min_coe;
   float fs_max_coe;
} t_cor_definition;/*}}}*/

// t_rct_definition structure
typedef struct {/*{{{*/
   int sizex;
   int sizey;
   float ee_coefficient;
   IM_STORE_TYPE ee_top_above_base;
   DIST_TYPE ee_thickness;
   IM_STORE_TYPE base_level;
   IM_STORE_TYPE base_level_variation;
   DIST_TYPE lsize;
   DIST_TYPE tsize;
   float rotation;
   float fs_density;
   DIST_TYPE fs_min_r;
   DIST_TYPE fs_max_r;
   float fs_min_coe;
   float fs_max_coe;
} t_rct_definition;/*}}}*/

// t_crs_definition structure
typedef struct {/*{{{*/
   int sizex;
   int sizey;
   float ee_coefficient;
   IM_STORE_TYPE ee_top_above_base;
   DIST_TYPE ee_thickness;
   IM_STORE_TYPE base_level;
   IM_STORE_TYPE base_level_variation;
   DIST_TYPE lsize;
   DIST_TYPE tsize;
   DIST_TYPE distance; 
   float rotation;
   float fs_density;
   DIST_TYPE fs_min_r;
   DIST_TYPE fs_max_r;
   float fs_min_coe;
   float fs_max_coe;
} t_crs_definition;/*}}}*/

enum{
   AIG_MSG_PRINT_COMMENT=1000, // For debugging, these messages will mostly be supressed.
   AIG_MSG_CREATING,
   AIG_MSG_APPLYING,
   AIG_MSG_PAINTING,
   AIG_MSG_SUCCESS,
   AIG_MSG_OOPS,
   AIG_MSG_FATAL_ERROR,
   AIG_MSG_LUA_ERROR,
   AIG_MSG_SAVING
};

enum{ // object IDs.
   AIG_ID_OBJECT = 9000,
   AIG_ID_VECTOR,
   AIG_ID_LINE,
   AIG_ID_TRIANGLE,
   AIG_ID_CURVE,
   AIG_ID_STRAIGHTLINE,
   AIG_ID_BEZIER,
   AIG_ID_VORONOI,
   AIG_ID_FEATURE,
   AIG_ID_GENERICFEATURE,
   AIG_ID_GOLDENGRAIN,
   AIG_ID_RECTANGLE,
   AIG_ID_SNAKE,
   AIG_ID_CORNER,
   AIG_ID_CROSS,
   AIG_ID_EFFECT,
   AIG_ID_EDGEEFFECT,
   AIG_ID_FINESTRUCTUREEFFECT,
   AIG_ID_SAMPLE,
   AIG_ID_GOLDONCARBONSAMPLE,
   AIG_ID_PERIODICCORNERSAMPLE,
   AIG_ID_SINGLERECTANGLESAMPLE,
   AIG_ID_SNAKESAMPLE,
   AIG_ID_PERIODICCROSSSAMPLE,
   AIG_ID_IMAGE,
   AIG_ID_IMAGEEFFECT,
   AIG_ID_PSF,
   AIG_ID_GAUSSIANPSF,
   AIG_ID_VIBRATION,
   AIG_ID_NOISE,
   AIG_ID_POISSONNOISE,
   AIG_ID_GAUSSIANNOISE,
   AIG_ID_BACKGROUND,
   AIG_ID_EVENBACKGROUND,
   AIG_ID_WAVYBACKGROUND
};

enum {
   AIG_CHECK_ID_NO_MATCH = 0,
   AIG_CHECK_ID_INHERITED,
   AIG_CHECK_ID_EXACT
};

// t_message definition
typedef struct { // This is a type for messages sent by the individual classes to the application./*{{{*/
   	  // It's good for progress bars and other progress or debug related
	  // meters.
   const char *sender_id;
   int message;
   const char *comment;
} t_message;/*}}}*/


#ifdef __cplusplus
extern "C" {
#endif
void *generate_gc_sample(t_gc_definition *def);
void destroy_gc_sample(void *poi);
void *generate_standard_image(void *sample, t_std_image_def_struct *def);
void destroy_image(void *poi);
void save_image(void *image, char* filename, char *comment);
#ifdef __cplusplus
}
#endif

#ifdef HAVE_LUA
int exec_lua_file(const char *fn);
#endif

/////////////// Footer - do not write below this line //////////////
#endif
// vim: cindent

