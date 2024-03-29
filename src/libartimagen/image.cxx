/* 
    Artificial Charged-Particle-Microscope Image Generator (ARTIMAGEN)
    2009  Petr Cizmar @ National Institute of Standards and Technology
    E-mail: petr.cizmar@nist.gov

    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to GNU.org public domain is compatible with GPL.

 */
#include <sstream>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <assert.h>
#include "strings.h"
#include "artimagen_i.h"
#include <config.h>

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

using namespace std;
using namespace artimagen;

///////////////// CImage ////////////////////////////

CImage::CImage(IM_COORD_TYPE sizex, IM_COORD_TYPE sizey){ /*{{{*/
   sender_id = "CImage";
   ident(AIG_ID_IMAGE);
   send_message(AIG_MSG_CREATING,"Empty image");

   //this is a constructor which allocates the buffer space and stores the x-
   //and y- sizes of the image.

   this->buffer = new IM_STORE_TYPE[sizex * sizey];
   this->sizex = sizex;
   this->sizey = sizey;
   fft_valid = 0;
   fft_data = NULL;
   plan_initialized = 0;
   block_ifft = IM_NO_IFFT_BLOCK;
}/*}}}*/

CImage::CImage(CImage *im){/*{{{*/
   sender_id = "CImage";
   ident(AIG_ID_IMAGE);
   send_message(AIG_MSG_CREATING,"Copying image");

   this->sizex = im->give_sizex();
   this->sizey = im->give_sizey();
   this->fft_data = NULL;
   this->buffer = new IM_STORE_TYPE[sizex * sizey];
   
   memcpy(buffer, im->give_buffer(), sizex * sizey * sizeof(IM_STORE_TYPE));

   plan_initialized = 0;
   fft_valid = 0;
   fft_data = NULL;
   block_ifft = IM_NO_IFFT_BLOCK;

   if (im->is_fft_valid()) {
      fft_valid = 1;
      make_fft_plan();
      memcpy(fft_data, im->give_fft_data(), sizex * sizey * sizeof(fftw_complex));
   }

}/*}}}*/

CImage::CImage(CImage *im, IM_COORD_TYPE sizex, IM_COORD_TYPE sizey){/*{{{*/
   sender_id = "CImage";
   send_message(AIG_MSG_CREATING,"Resizing image");
   // reszizing constructor
   this->sizex = sizex;
   this->sizey = sizey;
   this->fft_data = NULL;
   this->buffer = new IM_STORE_TYPE[sizex * sizey];

   fft_valid = 0;
   fft_data = NULL;
   plan_initialized = 0;
   block_ifft = IM_NO_IFFT_BLOCK;

   make_fft_plan();
   memset(fft_data, 0, sizex * sizey * sizeof(fftw_complex));

   im->calculate_fft();
   fftw_complex *im_fft_data = im->give_fft_data();

   memset(buffer, 0, sizex * sizey * sizeof(IM_STORE_TYPE));
   
   /* FFT image copying cycle */
   IM_COORD_TYPE msizex = (sizex < im->sizex)? sizex : im->sizex; //minimum x-size
   IM_COORD_TYPE msizey = (sizey < im->sizey)? sizey : im->sizey; //minimum y-size

   for (IM_COORD_TYPE j=0; j<msizey/4; j++)
      for (IM_COORD_TYPE i=0; i<msizex/4; i++){
	 fft_data[index(i,j,sizex)][0] = im_fft_data[index(i,j,im->sizex)][0];
	 fft_data[index(i,j,sizex)][1] = im_fft_data[index(i,j,im->sizex)][1];

	 fft_data[index(sizex-i-1,j,sizex)][0] = im_fft_data[index(im->sizex-i-1,j,im->sizex)][0];
	 fft_data[index(sizex-i-1,j,sizex)][1] = im_fft_data[index(im->sizex-i-1,j,im->sizex)][1];

	 fft_data[index(i,sizey-j-1,sizex)][0] = im_fft_data[index(i,im->sizey-j-1,im->sizex)][0];
	 fft_data[index(i,sizey-j-1,sizex)][1] = im_fft_data[index(i,im->sizey-j-1,im->sizex)][1];
	 
	 fft_data[index(sizex-i-1,sizey-j-1,sizex)][0] = im_fft_data[index(im->sizex-i-1,im->sizey-j-1,im->sizex)][0];
	 fft_data[index(sizex-i-1,sizey-j-1,sizex)][1] = im_fft_data[index(im->sizex-i-1,im->sizey-j-1,im->sizex)][1];
      }
      
   calculate_ifft();

   IM_STORE_TYPE des_mean = im->give_mean_value();
   IM_STORE_TYPE mean = give_mean_value();
   for (long int i=0; i<sizex*sizey; i++) buffer[i] *= des_mean / mean;
}
/*}}}*/

CImage::CImage(const string filename){ /*{{{*/
   sender_id = "CImage";
   ident(AIG_ID_IMAGE);
   send_message(AIG_MSG_CREATING,"Reading image from TIFF.");

   uint32_t i, j, sizex_h, sizey_h;

   unsigned char *helper_buffer;
   TIFF *tif = TIFFOpen(filename.data(), "r");
   if (tif) {
      TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &sizex_h);
      TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &sizey_h);

      this->sizex = sizex_h;
      this->sizey = sizey_h;
      this->buffer = new IM_STORE_TYPE[sizex * sizey];
      
      assert(helper_buffer = (unsigned char*) malloc(sizeof(uint32_t) * sizex_h * sizey_h)); // 8-bit images only
      TIFFReadRGBAImage(tif, sizex_h, sizey_h, (uint32_t *)helper_buffer, 1);
      for (j=0; j < sizey_h; j++) 
	 for (i=0; i < sizex_h; i++) {
	    uint32_t v_source;
	    v_source = 4 * (sizex_h * j + i);
	    uint32_t v_dest;
	    v_dest = sizex_h * (sizey_h - 1 - j) + i;
	    buffer[v_dest]= (IM_STORE_TYPE) helper_buffer[v_source] / (uint8_t)(-1); //copy every 4. valye to the grayscale buffer
	 }
      free(helper_buffer);

   } else throw(AIG_EX_CANNOT_READ_TIFF);
    
   TIFFClose(tif);

   fft_valid = 0;
   fft_data = NULL;
   plan_initialized = 0;
   block_ifft = IM_NO_IFFT_BLOCK;
}/*}}}*/

CImage::~CImage(){/*{{{*/
   // This is the destructor, disposes the buffers and plans.
   destroy_fft_plan();
   delete [] this->buffer;
}/*}}}*/

void CImage::zero(){/*{{{*/

   // this is to zero the buffer

   //IM_COORD_TYPE x,y;

   memset(this->buffer, 0, this->sizex * this->sizey * sizeof(IM_STORE_TYPE));
}/*}}}*/

void CImage::tiff_write(const string filename, const string tiff_comment, int bits){/*{{{*/
   switch (bits){
      case BI_8B:
	 this->tiff_write<TYPE_8B>(filename, tiff_comment);
	 break;
      case BI_16B:
	 this->tiff_write<TYPE_16B>(filename, tiff_comment);
	 break;
   }

}/*}}}*/

template <class save_type>
void CImage::tiff_write(const string filename, const string tiff_comment){/*{{{*/

   CSaveBuffer<save_type> sb(this);
   // this function writes the image to the tiff file including the provided
   // comment.

   TIFF *output;

   if((output = TIFFOpen(filename.data(), "w")) == NULL){
      cerr << "Could not write the image" << endl;
      throw(AIG_EX_CANNOT_WRITE_TIFF);
   }

   // Write the tiff tags to the file
   TIFFSetField(output, TIFFTAG_IMAGEWIDTH, this->sizex);
   TIFFSetField(output, TIFFTAG_IMAGELENGTH, this->sizey);
   TIFFSetField(output, TIFFTAG_COPYRIGHT, tiff_comment.data());

   //TIFFSetField(output, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
   TIFFSetField(output, TIFFTAG_COMPRESSION, COMPRESSION_NONE); //uncompr. images
   
   TIFFSetField(output, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
   TIFFSetField(output, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
   TIFFSetField(output, TIFFTAG_BITSPERSAMPLE, sb.give_element_size()*8);
   TIFFSetField(output, TIFFTAG_SAMPLESPERPIXEL, 1);

   // Actually write the image
   if(TIFFWriteEncodedStrip(output, 0,(void *) sb.give_save_buffer(),
	    this->sizex * this->sizey * sb.give_element_size()) == 0){
      cerr << "Could not write image." << endl;
      throw(AIG_EX_CANNOT_WRITE_TIFF);
   }
   TIFFClose(output);

}/*}}}*/

IM_COORD_TYPE CImage::give_sizex(){/*{{{*/
   return this->sizex;
}/*}}}*/

IM_COORD_TYPE CImage::give_sizey(){/*{{{*/
   return this->sizey;
}/*}}}*/

IM_STORE_TYPE *CImage::give_buffer(){/*{{{*/
   return this->buffer;
}/*}}}*/

void CImage::dot(IM_COORD_TYPE x, IM_COORD_TYPE y, IM_STORE_TYPE color){/*{{{*/
   if ((x>=0) && (y>=0) && (x<sizex) && (y<sizey))
	buffer[y * sizex + x] = color;
   fft_valid = 0;
}/*}}}*/

int CImage::make_fft_plan(){/*{{{*/
   fft_data = (fftw_complex *) fftw_malloc(sizex * sizey * sizeof(fftw_complex));
   assert(fft_data);

#ifdef HAVE_FFTW3_THREADS
   int nthreads = ::AIGApp->get_number_of_threads();
   fftw_plan_with_nthreads(nthreads);
#endif
   plan_im = fftw_plan_dft_2d(sizey, sizex, fft_data, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);

#ifdef HAVE_FFTW3_THREADS
   fftw_plan_with_nthreads(nthreads);
#endif
   plan_imi = fftw_plan_dft_2d(sizey, sizex, fft_data, fft_data, FFTW_BACKWARD, FFTW_ESTIMATE);

   plan_initialized = 1;

   return 1;

}/*}}}*/

int CImage::destroy_fft_plan(){/*{{{*/
   if (plan_initialized) {
      fftw_destroy_plan(plan_imi);
      fftw_destroy_plan(plan_im);
   }
   if (fft_data) {
      fftw_free(fft_data);
      fft_data = NULL;
   }
   plan_initialized = 0;
   return 1;
}/*}}}*/

void CImage::calculate_fft(){/*{{{*/
   if (!plan_initialized) make_fft_plan();
   if (!fft_valid) {
      for (int i=0; i<sizex*sizey; i++) {
	 fft_data[i][0] = buffer[i];
	 fft_data[i][1] = 0;
      }
      fftw_execute(plan_im);
      fft_valid = 1;
   }
}/*}}}*/

void CImage::calculate_ifft(){/*{{{*/
   fftw_execute(plan_imi);
   for (int i=0; i<sizex*sizey; i++) buffer[i] = fft_data[i][0] / (sizex*sizey);
   fft_valid = 0;
}/*}}}*/

int CImage::is_fft_valid(){/*{{{*/
   return fft_valid;
}/*}}}*/

fftw_complex *CImage::give_fft_data(){/*{{{*/
   return fft_data;
}/*}}}*/

void CImage::set_ifft_blocking(unsigned char blk){/*{{{*/
   block_ifft = blk;
}/*}}}*/

unsigned char CImage::give_ifft_blocking(){/*{{{*/
   return block_ifft;
}/*}}}*/

void CImage::shift(double shx, double shy){/*{{{*/

   calculate_fft();
   
   for (IM_COORD_TYPE j=0; j<sizey; j++) {
      double y = j;
      if (y > sizey/2) y -= sizey;
      for (IM_COORD_TYPE i=0; i<sizex; i++) {
	 double x = i;
	 if (x > sizex/2) x -= sizex;
	 double q1 = -2*M_PI*(x*shx/sizex+y*shy/sizey);
	 GE_INDEX_TYPE ind = index(i,j,sizex);
	 double si = sin(q1);
	 double co = cos(q1);
	 fftw_complex h;
	 h[0] = fft_data[ind][0] * co - fft_data[ind][1] * si;
	 h[1] = fft_data[ind][1] * co + fft_data[ind][0] * si;
	 fft_data[ind][0] = h[0];
	 fft_data[ind][1] = h[1];
      }
   }
   if (!block_ifft) calculate_ifft();
}/*}}}*/

void CImage::crop(IM_COORD_TYPE x1, IM_COORD_TYPE y1, IM_COORD_TYPE x2, IM_COORD_TYPE y2){/*{{{*/
   for (IM_COORD_TYPE j=0; j < (y2 - y1); j++)
      for (IM_COORD_TYPE i=0; i < (x2 - x1); i++){
	 buffer[index(i,j,x2-x1)] = buffer[index(i+x1,j+y1,sizex)];
      }
   sizex = x2-x1;
   sizey = y2-y1;
   fft_valid = 0;
   destroy_fft_plan();
   
} /*}}}*/

IM_STORE_TYPE CImage::give_mean_value(){/*{{{*/
   IM_STORE_TYPE mean = 0;
   for (long int i=0; i<(sizex*sizey);i++) mean +=buffer[i];
   mean /= sizex*sizey;
   return mean;
}/*}}}*/

//////////////// CImageEffect ///////////////////////////

int CImageEffect::apply(CImage *im){/*{{{*/
   sender_id = "CImageEffect";
   ident(AIG_ID_IMAGEEFFECT);
   return 0;
}/*}}}*/

//////////////// CPsf ///////////////////////////////////

CPsf::CPsf(IM_COORD_TYPE sizex, IM_COORD_TYPE sizey){/*{{{*/
   sender_id = "CPsf";
   ident(AIG_ID_PSF);
   psf_data = NULL;
   image_data = NULL;
   image_fft_data = NULL;

   this->sizex = sizex;
   this->sizey = sizey;

   psf_data = new double[sizex * sizey];
   memset(psf_data, 0, sizex * sizey * sizeof(psf_data));

}/*}}}*/

int CPsf::calculate_fft_psf(){/*{{{*/
   send_message(AIG_MSG_CREATING,"Calculating PSF");
   fft_psf_data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *sizex * sizey);
   for (int i=0; i<sizex*sizey; i++){
      fft_psf_data[i][0] = psf_data[i];	
      fft_psf_data[i][1] = 0;	
   }

   fftw_plan pl = fftw_plan_dft_2d(sizey, sizex, fft_psf_data, fft_psf_data, FFTW_FORWARD, FFTW_ESTIMATE);

   fftw_execute(pl);
   fftw_destroy_plan(pl);
   return 1;
}/*}}}*/

int CPsf::convolve(CImage *im){/*{{{*/
   if ((im->give_sizex() != sizex) || (im->give_sizey() != sizey)){
      send_message(AIG_MSG_FATAL_ERROR, "Wrong image size - cannot convolve, make a new PSF!!!");
      return -1;
   }

   im->calculate_fft();
   image_fft_data = im->give_fft_data();
   
   for (IM_COORD_TYPE j=0; j<sizey; j++)
      for (IM_COORD_TYPE i=0; i<sizex; i++){
	 fftw_complex h;
	 h[0] = image_fft_data[index(i,j,sizex)][0] * fft_psf_data[index(i,j,sizex)][0] - 
	    image_fft_data[index(i,j,sizex)][1] * fft_psf_data[index(i,j,sizex)][1];
	 h[1] = image_fft_data[index(i,j,sizex)][0] * fft_psf_data[index(i,j,sizex)][1] + 
	    image_fft_data[index(i,j,sizex)][1] * fft_psf_data[index(i,j,sizex)][0];
	 image_fft_data[index(i,j,sizex)][0] = h[0];
	 image_fft_data[index(i,j,sizex)][1] = h[1];
      }
   if (im->give_ifft_blocking() == IM_NO_IFFT_BLOCK) im->calculate_ifft();
   return 1;
}/*}}}*/

int CPsf::apply(CImage *im){/*{{{*/
   send_message(AIG_MSG_APPLYING,"Convolution");
   return convolve(im);
}/*}}}*/

CPsf::~CPsf(){/*{{{*/
   if (fft_psf_data) fftw_free(fft_psf_data);
   if (psf_data) delete [] psf_data;
}/*}}}*/

//////////////// CGaussianPsf ///////////////////////////////////

CGaussianPsf::CGaussianPsf(IM_COORD_TYPE sizex, IM_COORD_TYPE sizey, float sigma, float astigmatism_ratio, float astigmatism_angle):CPsf(sizex, sizey){/*{{{*/

   sender_id = "CGaussianPsf";
   ident(AIG_ID_GAUSSIANPSF);
   IM_COORD_TYPE psfsize = sigma * astigmatism_ratio*3;

   for (IM_COORD_TYPE j=-psfsize; j<psfsize; j++)
      for (IM_COORD_TYPE i=-psfsize; i<psfsize; i++){
	 
	 IM_COORD_TYPE cx = i;
	 IM_COORD_TYPE cy = j;
	 
	 if (cx<0) cx += sizex;
	 if (cy<0) cy += sizey;
	 if (cx>=sizex) cx -= sizex;
	 if (cy>=sizey) cy -= sizey;

	 double xx=(i*cos(astigmatism_angle)+j*sin(astigmatism_angle))*astigmatism_ratio;
	 double yy=(-i*sin(astigmatism_angle)+j*cos(astigmatism_angle))/astigmatism_ratio;

	 psf_data[index(cx,cy,sizex)] = (double)1/2/sqrt(2*M_PI)/pow(sigma,2)*exp(-(pow(xx,2)+pow(yy,2))/2/pow(sigma,2));
      }

   calculate_fft_psf();
}/*}}}*/

//////////////// CVibration /////////////////////////////////////

CVibration::CVibration(t_vib_definition *def){/*{{{*/
   const_init();
   memcpy(&this->def, def, sizeof(t_vib_definition)); // stores the definitions
   int nf = def->number_of_frequencies;
   frequencies = new float[nf];
   phases = new float[nf];
   amplitudes_x = new float[nf];
   amplitudes_y = new float[nf];

   for (int i=0; i<nf; i++){
      frequencies[i] = rand() * (def->max_frequency - def->min_frequency) / RAND_MAX + def->min_frequency;
      phases[i] = rand() * 2* M_PI / RAND_MAX;
      double amp = rand() * def->max_amplitude / RAND_MAX;
      double ang = rand() * 2 * M_PI / RAND_MAX;
      amplitudes_x[i] = amp * cos(ang);
      amplitudes_y[i] = amp * sin(ang);
   }
}/*}}}*/

CVibration::CVibration(t_vib_definition *def, vector <float*> vibs){/*{{{*/
   const_init();
   memcpy(&this->def, def, sizeof(t_vib_definition)); // stores the definitions
   int nf = vibs.size();
   this->def.number_of_frequencies = nf;

   frequencies = new float[nf];
   phases = new float[nf];
   amplitudes_x = new float[nf];
   amplitudes_y = new float[nf];

   for (int i=0; i<nf; i++){
      frequencies[i] = vibs[i][2];
      phases[i] = vibs[i][3];
      amplitudes_x[i] = vibs[i][0];
      amplitudes_y[i] = vibs[i][1];
   }
}/*}}}*/

void CVibration::const_init(){/*{{{*/
   sender_id = "CVibration";
   ident(AIG_ID_VIBRATION);
   send_message(AIG_MSG_CREATING,"Calculating vibration functions");
   motion_x = NULL;
   motion_y = NULL;
}/*}}}*/

CVibration::~CVibration(){/*{{{*/
   delete [] frequencies;
   delete [] phases;
   delete [] amplitudes_x;
   delete [] amplitudes_y;
   if (motion_x) delete [] motion_x;
   if (motion_y) delete [] motion_y;
}/*}}}*/

void CVibration::generate_curves(const int sizex, const int sizey, unsigned long time_shift){/*{{{*/
   assert(!motion_x); // if the curves were already generated, this would cause a memory leak.
   assert(!motion_y); // if you want to regenerate curves, delete the original  first.
   motion_x = new double[(int) (sizex * sizey)];
   motion_y = new double[(int) (sizex * sizey)];

   unsigned long t = (double) time_shift/def.pixel_dwell_time;

   for (int j=0; j< sizey; j++){
      t += def.line_dead_time;
      for (int i=0; i< sizex; i++){
	 t += def.pixel_dwell_time+def.pixel_dead_time;
	 motion_x[index(i,j,sizex)] = 0;
	 motion_y[index(i,j,sizex)] = 0;
	 for (int f = 0; f<def.number_of_frequencies; f++){
	    motion_x[index(i,j,sizex)] += amplitudes_x[f] * sin(2 * M_PI * frequencies[f] * t * 1e-9 + phases[f]);
	    motion_y[index(i,j,sizex)] += amplitudes_y[f] * sin(2 * M_PI * frequencies[f] * t * 1e-9 + phases[f]);
	 }
	 // Part of debug output // debfile << motion_x[index(i,j,sizex)] << ", " << motion_y[index(i,j,sizex)] << endl;
      }
   }
   // Part of debug output // debfile.close();
}/*}}}*/

void CVibration::destroy_curves(){/*{{{*/
   if (motion_x) {
      delete [] motion_x;
      motion_x = NULL;
   }
   if (motion_y) {
      delete [] motion_y;
      motion_y = NULL;
   }
}/*}}}*/

int CVibration::apply(CImage *im, unsigned long time_shift){/*{{{*/
   send_message(AIG_MSG_APPLYING,"Applying vibration");
   const IM_COORD_TYPE sx = im->give_sizex();
   const IM_COORD_TYPE sy = im->give_sizey();

   IM_STORE_TYPE *helper_array = new IM_STORE_TYPE[sx*sy];

   if ((motion_x) || (motion_y)) destroy_curves(); // if the vib has been used, delete it
   generate_curves(sx,sy,time_shift);

   memcpy(helper_array, im->give_buffer(),sx*sy*sizeof(IM_STORE_TYPE));

   for (IM_COORD_TYPE j=0; j<sy-1; j++)
      for (IM_COORD_TYPE i=0; i<sx-1; i++){
	 GE_INDEX_TYPE ind = index(i, j, sx);
	 if ((i+motion_x[ind] >= 0) && (i+motion_x[ind] < sx) && (j+motion_y[ind] >= 0) && (j+motion_y[ind] < sy)){
	    double x = motion_x[ind]; 
	    double y = motion_y[ind]; 
	    double wx = x - floor(x); // x-weight (for interpolation) of the vibration vector ( = (motion_x[ind], motion_y[ind]))
	    double wy = y - floor(y); // y-weight 
	    double fq1 = (helper_array[index(i+floor(x),j+floor(y),sx)]*(1-wy)+helper_array[index(i+floor(x),j+ceil(y),sx)]*wy);
	    double fq2 = (helper_array[index(i+ceil(x),j+floor(y),sx)]*(1-wy)+helper_array[index(i+ceil(x),j+ceil(y),sx)]*wy);
	    double fp = (fq1 * (1-wx) + fq2 * wx);
	    im->dot(i, j, fp);
	 }
      }
   delete [] helper_array;
   return 0;
}/*}}}*/

/////////////// CNoise /////////////////////////////////////////

CNoise::CNoise(){/*{{{*/
   sender_id = "CNoise";
   ident(AIG_ID_NOISE);
}/*}}}*/

int CNoise::apply(CImage *im){/*{{{*/
   send_message(AIG_MSG_APPLYING,"Applying");
   const IM_COORD_TYPE sx = im->give_sizex();
   const IM_COORD_TYPE sy = im->give_sizey();

   for (IM_COORD_TYPE j=0; j<sy; j++)
      for (IM_COORD_TYPE i=0; i<sx; i++){
	 IM_STORE_TYPE v = im->give_buffer()[index(i,j,sx)];
	 v = noise_value(v);
	 if (v>1) v=1;  // clipping is better than overflow
	 if (v<0) v=0; 
	 im->dot(i, j, v);
      }
   return 1;
}/*}}}*/

double CNoise::noise_value(double n){/*{{{*/
   return 0;
}/*}}}*/

//////////////// CPoissonNoise /////////////////////////////////

CPoissonNoise::CPoissonNoise(double particles_per_unit){/*{{{*/
   sender_id = "CPoissonNoise";
   ident(AIG_ID_POISSONNOISE);
   if (particles_per_unit > 1000) throw (int) AIG_EX_POISSON_TOO_HIGH_LAMBDA;
   this->particles_per_unit = particles_per_unit;
}/*}}}*/

long int factorial(long int n){
   if (n > 1) return n*factorial(n-1);
   return 1;
}

double CPoissonNoise::noise_value(double n){/*{{{*/
   double lambda = n*particles_per_unit;

   if (n < 30) {
      double L = exp(-lambda);
      int k = 0;
      double p = 1;
      do {
	 k++;
	 p *= (double)rand()/RAND_MAX;
      } while (p >= L);
      return (k-1)/particles_per_unit;
   } else {


      double c = 0.767 - 3.36/lambda;
      double beta = 3.14159265/sqrt(3.0*lambda);
      double alpha = beta*lambda;
      double k = log(c) - lambda - log(beta);

      double x,y,u,v,lhs,rhs;
      int n;
      do {
	 u = rand()/RAND_MAX;
	 x = (alpha - log((1.0 - u)/u))/beta;
	 n = floor(x + 0.5);
	 if (n >= 0){
	    v = rand()/RAND_MAX;
	    y = alpha - beta*x;
	    lhs = y + log(v/pow((1.0 + exp(y)),2));
	    rhs = k + n*log(lambda) - log(factorial(n));

	 } 
      } while (lhs > rhs);
      return n;
   }
}/*}}}*/

//////////////// CGaussianNoise ///////////////////////////////

CGaussianNoise::CGaussianNoise(double sigma){/*{{{*/
   sender_id = "CGaussianNoise";
   ident(AIG_ID_GAUSSIANNOISE);
   this->sigma = sigma;
}/*}}}*/

double CGaussianNoise::noise_value(double n){/*{{{*/
   double u1;
   do{
      u1 = (double) rand() / RAND_MAX;
   } while (u1 == 0);
   double u2 = (double) rand() / RAND_MAX;

   return n + sigma * sqrt(-2 * log(u1)) * cos(2*M_PI*u2);
}/*}}}*/


//////////////// CBackgroud //////////////////////////////////

CBackgroud::CBackgroud(){/*{{{*/
   sender_id = "CBackgroud";
   ident(AIG_ID_BACKGROUND);
}/*}}}*/

int CBackgroud::apply(CImage *im){ /*{{{*/
   return 0;
}/*}}}*/


//////////////// CEvenBackgroud //////////////////////////////

CEvenBackgroud::CEvenBackgroud(IM_STORE_TYPE gl){/*{{{*/
   sender_id = "CEvenBackgroud";
   ident(AIG_ID_EVENBACKGROUND);
   send_message(AIG_MSG_CREATING,"Creating even background");
   bg_greylevel = gl;
}/*}}}*/

int CEvenBackgroud::apply(CImage *im){ /*{{{*/
   send_message(AIG_MSG_APPLYING,"Applying");
   const IM_COORD_TYPE sx = im->give_sizex();
   const IM_COORD_TYPE sy = im->give_sizey();

   for (IM_COORD_TYPE j=0; j<sy; j++)
      for (IM_COORD_TYPE i=0; i<sx; i++){
	 im->dot(i, j, bg_greylevel);
      }

   return 1;
}/*}}}*/

IM_STORE_TYPE CEvenBackgroud::give_bg_greylevel(){/*{{{*/
   return bg_greylevel;
}/*}}}*/

//////////////// CWavyBackgroud //////////////////////////////

CWavyBackgroud::CWavyBackgroud(IM_STORE_TYPE min_gl, IM_STORE_TYPE max_gl, int densx, int densy){/*{{{*/
   const_init(densx, densy);
   this->min_gl = min_gl;
   this->max_gl = max_gl;
   for (int i=0; i< densx*densy; i++) bg_im->give_buffer()[i] = min_gl + rand() * (max_gl - min_gl) / RAND_MAX;
}/*}}}*/

CWavyBackgroud::CWavyBackgroud(int densx, int densy, vector <IM_STORE_TYPE> values){/*{{{*/
   this->const_init(densx, densy);
   int j = 0;
   for (int i=0; i< densx*densy; i++) {
      bg_im->give_buffer()[i] = values[j++];
   }
}/*}}}*/

void CWavyBackgroud::const_init(int densx, int densy){/*{{{*/
   sender_id = "CWavyBackgroud";
   ident(AIG_ID_WAVYBACKGROUND);
   send_message(AIG_MSG_CREATING,"Creating wavy background");
   this->densx = densx;
   this->densy = densy;
   bg_im = new CImage(densx,densy);
}/*}}}*/

CWavyBackgroud::~CWavyBackgroud(){/*{{{*/
   delete bg_im;
}/*}}}*/

int CWavyBackgroud::apply(CImage *im){/*{{{*/
   send_message(AIG_MSG_APPLYING,"Applying");
   CImage *bg_im_big = new CImage(bg_im, im->give_sizex(), im->give_sizey());
   for (long int i=0; i< (im->give_sizex() * im->give_sizey()); i++) im->give_buffer()[i] = bg_im_big->give_buffer()[i]; 
   delete bg_im_big;
   return 0;
}/*}}}*/

// vim: cindent
