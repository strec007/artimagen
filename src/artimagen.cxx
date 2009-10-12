/* 
    Artificial SEM Image Generator (ArtImaGen)
    2009  Petr Cizmar @ National Institute of Standards and Technology
    E-mail: petr.cizmar@nist.gov

    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to Gnu.org public domain is compatible with GPL.

 */

#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "libartimagen/artimagen_i.h"
#include "libartimagen/artimagen.h"
#include "../config.h"

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
 
#define SIZEX 1124
#define SIZEY 868

using namespace std;
using namespace artimagen;

void print_version_info(){
   cout << "\n **** ArtImaGen v." << PACKAGE_VERSION << " ****\n * by Petr Cizmar @ NIST\n * 2008-2009\n\n";
}

void initialize(){
   srandom(time(0));
   print_version_info();
}

void generate_gold_on_carbon_image(CSample **sam){
   t_gc_definition gcdef;

   gcdef.sizex=SIZEX;
   gcdef.sizey=SIZEY;
   gcdef.ee_coefficient=0.2;
   gcdef.ee_top_above_basic=0.3;
   gcdef.ee_thickness=10;
   gcdef.basic_level=0.5;
   gcdef.basic_level_variation=0.1;
   gcdef.grain_min_size=5;
   gcdef.grain_max_size=20;
   gcdef.number_of_grains=SIZEX*SIZEY/(50*50)/2;
   gcdef.fs_density = 7e-3;
   gcdef.fs_min_r = 6;
   gcdef.fs_max_r = 10;
   gcdef.fs_min_coe = 0.9;
   gcdef.fs_max_coe = 1.1;
   *sam = new CGoldOnCarbonSample(&gcdef);
}

void generate_corner_structure_image(CSample **sam){
   t_cor_definition cordef;

   cordef.sizex=SIZEX;
   cordef.sizey=SIZEY;
   cordef.ee_coefficient=0.01;
   cordef.ee_top_above_basic=0.1;
   cordef.ee_thickness=20;
   cordef.basic_level=0.5;
   cordef.basic_level_variation=0.2;
   cordef.lsize=50;
   cordef.tsize=150;
   cordef.distance=200;
   cordef.fs_density = 7e-4;
   cordef.fs_min_r = 6;
   cordef.fs_max_r = 10;
   cordef.fs_min_coe = 0.95;
   cordef.fs_max_coe = 1.05;
   *sam = new CPeriodicCornerSample(&cordef);
}

void generate_cross_structure_image(CSample **sam){
   t_crs_definition crsdef;

   crsdef.sizex=SIZEX;
   crsdef.sizey=SIZEY;
   crsdef.ee_coefficient=0.2;
   crsdef.ee_top_above_basic=0.5;
   crsdef.ee_thickness=10;
   crsdef.basic_level=0.3;
   crsdef.basic_level_variation=0.2;
   crsdef.lsize=70;
   crsdef.tsize=20;
   crsdef.distance=200;
   crsdef.rotation=0;
   crsdef.fs_density = 2e-2;
   crsdef.fs_min_r = 2;
   crsdef.fs_max_r = 6;
   crsdef.fs_min_coe = 0.9;
   crsdef.fs_max_coe = 1.1;
   try {
      CSample *sample;
      sample = new CPeriodicCrossSample(&crsdef);
      sample->move_by(CVector(70,70));
      *sam = sample;
   }
   catch (int exc){
      if (exc == AIG_EX_FEATURE_OVERLAP) 
	 cout << "Bad configuration, features overlap" << endl;
   }
}

void generate_rectangle_structure_image(CSample **sam){
   t_rct_definition rctdef;

   rctdef.sizex=SIZEX;
   rctdef.sizey=SIZEY;
   rctdef.ee_coefficient=0.2;
   rctdef.ee_top_above_basic=0.5;
   rctdef.ee_thickness=10;
   rctdef.basic_level=0.3;
   rctdef.basic_level_variation=0.2;
   rctdef.lsize=300;
   rctdef.tsize=200;
   rctdef.rotation=30*3.1416/180;
   rctdef.fs_density = 2e-2;
   rctdef.fs_min_r = 2;
   rctdef.fs_max_r = 6;
   rctdef.fs_min_coe = 0.9;
   rctdef.fs_max_coe = 1.1;
   *sam = new CSingleRectangleSample(&rctdef);
}


void generate_snake_structure_image(CSample **sam){
   t_rct_definition rctdef;

   rctdef.sizex=SIZEX;
   rctdef.sizey=SIZEY;
   rctdef.ee_coefficient=0.2;
   rctdef.ee_top_above_basic=0.5;
   rctdef.ee_thickness=10;
   rctdef.basic_level=0.3;
   rctdef.basic_level_variation=0.02;
   rctdef.lsize=1200;
   rctdef.tsize=50;
   rctdef.rotation=0*3.1415927/180;
   rctdef.fs_density = 2e-2;
   rctdef.fs_min_r = 2;
   rctdef.fs_max_r = 6;
   rctdef.fs_min_coe = 0.9;
   rctdef.fs_max_coe = 1.1;
   *sam = new CSnakeSample(&rctdef);
}

int main(int argc, char **argv){


   CSample *sam = NULL;
   initialize();
   //try{ 
      CImage im(SIZEX,SIZEY);

      //CEvenBackgroud back(0.2);
      CWavyBackgroud back(0.2,0.4,5,5);
      back.apply(&im);

      generate_snake_structure_image(&sam);
      //generate_rectangle_structure_image(&sam);
      //generate_corner_structure_image(&sam);
      //generate_cross_structure_image(&sam);
      sam->paint(&im);

      im.set_ifft_blocking(IM_IFFT_BLOCK);
      CGaussianPsf psf(SIZEX,SIZEY,0.7,2,30); //sizex, sizey, sigma, astig. ratio, ast. angle.
      psf.apply(&im);
//      im.shift(0.5,0.5);
      im.calculate_ifft();

      // create the drift-distortion function and apply it
      t_vib_definition vibdef;
      vibdef.pixel_dwell_time = 100;
      vibdef.min_frequency = 1;
      vibdef.max_frequency = 100;
      vibdef.max_amplitude = 0.2;
      vibdef.number_of_frequencies = 10;
      vibdef.pixel_dead_time = 25; //estimate values
      vibdef.line_dead_time = 200; //estimate

      CVibration vib(&vibdef);
      vib.apply(&im, 0);
      im.crop(50,50,SIZEX-50,SIZEY-50);
      CGaussianNoise gn(0.07);
      gn.apply(&im);

      im.tiff_write("test.tiff","Generated by ArtImaGen software by Petr Cizmar @ NIST", BI_8B);
}

// vim: cindent
