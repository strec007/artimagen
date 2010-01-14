/* 
    Artificial SEM Image Generator (ArtImaGen)
    2009  Petr Cizmar @ National Institute of Standards and Technology
    E-mail: petr.cizmar@nist.gov

    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to Gnu.org public domain is compatible with GPL.

 */

#ifndef ARTIMAGEN_I_H
#define ARTIMAGEN_I_H

#include <cstring>
#include <string>
#include <fftw3.h>
#include <tiffio.h>
#include <cmath>
#include <vector>
#include "artimagen.h"

#define IM_NO_IFFT_BLOCK 0
#define IM_IFFT_BLOCK 1

enum {
CU_SEGMENT_NOT_HIT, 
CU_SEGMENT_HIT, 
CU_SEGMENT_HIT_ERR,
CU_SEGMENT_HIT_VERTEX
};

enum {
   AIG_FE_ARRAY_FULL,
   AIG_FE_CURVE_INSERTED_OK
};


#define SA_ADD_OVERLAP 0
#define SA_ADD_OK 1
#define SA_ADD_FULL 2

enum {
   AIG_EX_ZERO_LINE,
   AIG_EX_READD_EFFECT_CHAIN_ATTEMPT,
   AIG_EX_FEATURE_OVERLAP,
   AIG_EX_GOLDONCARBON_TOO_MANY_FAILS
};


namespace artimagen {
using namespace std;

/* object */

class CApp{/*{{{*/
   public:
      CApp();
      void set_message_call_back(void (*f)(t_message *));
      void broadcast_message(t_message message);
   private:
      void (*call_back)(t_message *);
};/*}}}*/

#ifndef AIGAPP_DECLARATION
extern CApp *AIGApp;
#endif

class CObject{/*{{{*/
   public:
      CObject();
      void disable_messages();
      void enable_messages();
      int check_id(int id);
   protected:
      virtual void send_message(int message, const char *comment);
      const char* sender_id;
      char messaging_enabled;
      std::vector<int> inherited_ids;
      void ident(int); // make yourself an ID;
   private:
      int object_id;
};/*}}}*/

/* geometry */

class CCurve;
class CVector;
class CImage;

inline GE_INDEX_TYPE index(unsigned int x, unsigned int y, unsigned int sizex){
   return y * sizex + x;
}

template <class T>
T sq(T a){
   return (a*a);
}

DIST_TYPE cross(CVector a, CVector b);
DIST_TYPE dot(CVector a, CVector b);

class CVector{/*{{{*/
   public:
      CVector();
      CVector(DIST_TYPE x, DIST_TYPE y);
      DIST_TYPE x;
      DIST_TYPE y;
      CVector operator + (CVector a);
      CVector operator - (CVector a);
      CVector operator = (CVector a);
      CVector operator * (DIST_TYPE a);
      DIST_TYPE sq_length();
      int is_colinear(CVector a);
      string p(); // print vector - fore debugging
      void paint(CImage *im);
      void rotate(CVector center, double angle);
};/*}}}*/

class CLine{/*{{{*/
   public:
      CLine();
      CLine(CVector p0, CVector p1);
      CLine operator = (CLine a); // copy operator
      bool operator || (CLine a); // parallel relation
      bool operator ^ (CLine a); // overlap relation
      double operator % (CLine a); // intersection relation
      int intersection(CLine a, double *tt, double *uu); // intersection calculation function
      DIST_TYPE sq_distance_from(CVector v);
      DIST_TYPE sq_length();
      string p(); // print line coordinates - debugging
      void paint(CImage *im);
      CVector parametrize(double t);
   private:
      CVector p0;
      CVector p1;

};/*}}}*/

class CTriangle:public CObject{/*{{{*/
   public:
      CTriangle(CVector va, CVector vb, CVector vc);

   private:
      CVector va;
      CVector vb;
      CVector vc;
      CLine ea;
      CLine eb;
      CLine ec;
};/*}}}*/

class CCurve:public CObject{/*{{{*/
   public:
      CCurve();
      ~CCurve();
      DIST_TYPE sq_distance_from(CVector v);
      void calculate_bounding_box();
      void give_bounding_box(CVector *tl, CVector *br);
      int give_number_of_vertices();
      CVector *give_vertices();
      void paint(CImage *im);
      virtual int is_hit(CVector x, CVector d); // tells number of hits
      int hits_line(CLine l); // tells number of hits, -1 if error.
      virtual void move_by(CVector mv);
      virtual void rotate(CVector center, double angle);
   protected:
      virtual void init();
      void approximate_by_vertices();
      int number_of_vertices;
      virtual CVector parametrize(DIST_TYPE t); //t should be from 0 to 1
      CVector *vertices;
      CVector bounding_box_tl;
      CVector bounding_box_br;
      int is_segment_hit(CVector p1, CVector p2, CVector x, CVector d); // p1 and p2 are defining the segment and x is the starting point, d is the direction of test
   private:

};/*}}}*/

class CStraightLine:public CCurve{/*{{{*/
   public:
      CStraightLine(CLine l);
   protected:
      virtual CVector parametrize(DIST_TYPE t);
   private:
      CLine line;
};/*}}}*/

class CBezier:public CCurve{/*{{{*/
   public:
      CBezier(CVector p1, CVector p2, CVector p3, CVector p4);
   protected:
      virtual CVector parametrize(DIST_TYPE t);
   private:
      CVector control_points[4];
};/*}}}*/

class CPolygon:public CCurve{/*{{{*/
   public:
      CPolygon();
      ~CPolygon();
      virtual int is_hit(CVector x, CVector d); // tells number of hits
      int is_inside(CVector p);
      void set_map_division(int md);
      bool overlaps(CPolygon *fe);
   protected:
      void create_map();
      void destroy_map();

   private:
      char *map;
      unsigned int map_division;
      int is_inside_with_hits(CVector p);

};/*}}}*/

class CVoronoi:public CObject{/*{{{*/
   public:
      CVoronoi(CPolygon *polygon, CVector *free_points, int num_of_free);
      ~CVoronoi();
      void calculate_map(int hsplit, int vsplit);
      void lloyd();
   private:
      int vsplit, hsplit;  // number of splits in vertical and horizontal directions
      CVector map2vec(int x, int y);
      int *map; // actual voronoi map saying which pixel is closest to the region.
      CPolygon *polygon;
      CVector *free_points;
      int number_of_free;
      CVector bounding_box_tl;
      CVector bounding_box_br;
};/*}}}*/

/* sample */

class CEffect;

class CFeature:public CPolygon{/*{{{*/
   public:
      CFeature();
      CFeature(vector <CCurve *> curve_vec, vector<CEffect *> effect_vec);
      virtual void init();
      ~CFeature();
      virtual void paint(CImage *im);
      void add_effect_chain(CEffect **effects, int number_of_effects);
      IM_STORE_TYPE give_base_gray_level();
      void set_base_gray_level(IM_STORE_TYPE level);
   protected:
      int number_of_curves;
      int actual_number_of_curves;
      int number_of_effects;
      CCurve **curves;
      CEffect **effects;
      IM_STORE_TYPE base_gray_level;
      void build_vertices();
   private:
      int add_curve(CCurve *cu);
      void const_init();
};/*}}}*/

class CGoldenGrain:public CFeature{/*{{{*/
   public:
      CGoldenGrain(const CVector position, const DIST_TYPE size);
   protected:
};/*}}}*/

class CRectangle:public CFeature{/*{{{*/
   public:
      CRectangle(const DIST_TYPE lsize, const DIST_TYPE tsize, double rotation);
};/*}}}*/

class CRectangleRounded:public CFeature{/*{{{*/
   public:
      CRectangleRounded(const DIST_TYPE lsize, const DIST_TYPE tsize, double rotation, const DIST_TYPE rounding);
};/*}}}*/

class CSnakeRounded:public CFeature{/*{{{*/
   public:
      CSnakeRounded(const DIST_TYPE w, const DIST_TYPE a, const DIST_TYPE b, const DIST_TYPE c, double rotation);
};/*}}}*/

class CSnake:public CFeature{/*{{{*/
   public:
      CSnake(const DIST_TYPE w, const DIST_TYPE a, const DIST_TYPE b, const DIST_TYPE c, double rotation);
};/*}}}*/

class CCorner:public CFeature{/*{{{*/
   public:
      CCorner(const DIST_TYPE lsize, const DIST_TYPE tsize, double rotation);
};/*}}}*/

class CCross:public CFeature{/*{{{*/
   public:
      CCross(const DIST_TYPE lsize, const DIST_TYPE tsize, double rotation);
};/*}}}*/

class CEffect:public CObject{/*{{{*/
   public:
      CEffect();
      float give_amplification(CFeature *fe, CVector v); // effect amplifies the original value depenging on coordinates
   protected:
   private:
      virtual double fun(CFeature *fe, CVector v) = 0;
};/*}}}*/

class CEdgeEffect: public CEffect{/*{{{*/
   public:
      CEdgeEffect(float coefficient, IM_STORE_TYPE top_edge_value_above_base);
   protected:
   private:
      double fun(CFeature *fe, CVector v);
      float coefficient;
      IM_STORE_TYPE top_edge_value_above_base;
      IM_STORE_TYPE center_feature_value;

};/*}}}*/

typedef struct {/*{{{*/
	DIST_TYPE r;
	CVector c;
	float coe;
} finestruct_spot;/*}}}*/

class CFineStructureEffect: public CEffect{/*{{{*/
   public:
      CFineStructureEffect(float density, DIST_TYPE min_r, DIST_TYPE max_r, float min_coe, float max_coe);
      ~CFineStructureEffect();
   private:
      float density; // in spots per square pixel
      DIST_TYPE min_r; // min and max radius
      DIST_TYPE max_r; 
      float min_coe; // min and max amplification coefficient
      float max_coe;
      DIST_TYPE sizex, sizey; // size
      finestruct_spot *spots; // array of spots
      int nos; // number of spots;

      void generate_spots(DIST_TYPE sizex, DIST_TYPE sizey);
      virtual double fun(CFeature *fe, CVector v);
};/*}}}*/

class CSample:public CObject{/*{{{*/
   public:
      CSample(DIST_TYPE sizex, DIST_TYPE sizey);
      CSample(DIST_TYPE sizex, DIST_TYPE sizey, vector<CFeature *> feature_vec);
      ~CSample();
      void paint(CImage *im);
      void move_by(CVector mv);
      int add_feature(CFeature *fe); //adds feature if it does not conflict with others 
   protected:
      CFeature **features; // list of features
      int number_of_features;
      CEffect **effects; // list of effects applied
      int number_of_effects;

      DIST_TYPE sizex, sizey; // size of the sample
      void create_map();
      void destroy_map();
      void occupies_map(CFeature *fe, int *map); // fills 1 into the map for occupied elemnts, otherwise fills in 0
      void add_feature_to_map(int i);
      int overlaps(CFeature *fe);
      int actual_number_of_features;
      virtual void rotate(CVector center, double angle);
   private:
      int map_division;
      unsigned int **map; // map of layout - to speed up the overlapping detection 
      int *map_counts; // count of features liying in each map segment
      void add_feature_chain(CFeature **features, int number_of_features);//debug only
      void const_init(DIST_TYPE sizex, DIST_TYPE sizey);
};/*}}}*/

class CGoldOnCarbonSample:public CSample{/*{{{*/
   public:
   CGoldOnCarbonSample(t_gc_definition *def);
   private:
};/*}}}*/

class CPeriodicCornerSample:public CSample{/*{{{*/
   public:
      CPeriodicCornerSample(t_cor_definition *def);
   private:
};/*}}}*/

class CSingleRectangleSample:public CSample{/*{{{*/
   public:
      CSingleRectangleSample(t_rct_definition *def);
   private:
};/*}}}*/

class CPeriodicCrossSample:public CSample{/*{{{*/
   public:
      CPeriodicCrossSample(t_crs_definition *def);
   private:
};/*}}}*/

class CSnakeSample:public CSample{/*{{{*/
   public:
      CSnakeSample(t_rct_definition *def);
   private:
};/*}}}*/

/* image */


class CImage:public CObject{/*{{{*/
   public:
      CImage(IM_COORD_TYPE sizex, IM_COORD_TYPE sizey); //create empty image of a certain size
      CImage(CImage *im); // copy image
      CImage(CImage *im, IM_COORD_TYPE sizex, IM_COORD_TYPE sizey); // copy and resize image
      ~CImage(); // desroy image
      void tiff_write(const string filename, const string tiff_comment, int bits); // write image file in a tiff file
      IM_COORD_TYPE give_sizex();
      IM_COORD_TYPE give_sizey();
      IM_STORE_TYPE *give_buffer();
      void dot(IM_COORD_TYPE x, IM_COORD_TYPE y, IM_STORE_TYPE color);
      fftw_complex *give_fft_data();
      void calculate_fft();
      void calculate_ifft();
      int is_fft_valid();
      void set_ifft_blocking(unsigned char blk);  // WARNING: With ifft_blocking the image and 
                                                //its FFT representation may not be consistent!!!
                                                // user is responsible for make the ifft in the end of 
                                                // fftw-based operations
      unsigned char give_ifft_blocking();
      void shift(double shx, double shy);
      void crop(IM_COORD_TYPE x1, IM_COORD_TYPE y1, IM_COORD_TYPE x2, IM_COORD_TYPE y2);

      IM_STORE_TYPE give_mean_value();

   private:
      int make_fft_plan();
      int destroy_fft_plan();
      template <class save_type>
      void tiff_write(const string filename, const string tiff_comment); 

      fftw_plan plan_im, plan_imi;

      IM_COORD_TYPE sizex, sizey;
      IM_STORE_TYPE *buffer;
      fftw_complex *fft_data;

      unsigned char plan_initialized, fft_valid;
      unsigned char block_ifft;                 // ifft is not automatically calculated. good for 
                                                //saving machine-time when the following steps all work with fft
      
      void zero();

};/*}}}*/

template <class save_type>
class CSaveBuffer:public CObject{/*{{{*/
      public:
      CSaveBuffer(CImage *im);
      ~CSaveBuffer();
      save_type convert_to_save_type(IM_STORE_TYPE x);
      size_t give_element_size();
      void *give_save_buffer();
      private:
      save_type *save_buffer;
};

///////////////// implementation ////////////////////////////
//must be in the header file, since it's a template

template <class save_type>
CSaveBuffer<save_type>::CSaveBuffer(CImage *im){/*{{{*/
   this->save_buffer = new save_type[im->give_sizex() * im->give_sizey()];

   IM_COORD_TYPE i,j;
   IM_STORE_TYPE *buf = im->give_buffer();
   for (j=0; j<im->give_sizey(); j++)
      for (i=0; i<im->give_sizex(); i++){
	 this->save_buffer[j * im->give_sizex() + i] = this->convert_to_save_type(buf[j * im->give_sizex() + i]);
      }
}/*}}}*/

template <class save_type>
CSaveBuffer<save_type>::~CSaveBuffer(){/*{{{*/
   delete [] this->save_buffer;
}/*}}}*/

template <class save_type>
size_t CSaveBuffer<save_type>::give_element_size(){/*{{{*/
   return sizeof(save_type);
}/*}}}*/

template <class save_type>
void *CSaveBuffer<save_type>::give_save_buffer(){/*{{{*/
   return (void *) this->save_buffer;
}/*}}}*/

template <class save_type>
save_type CSaveBuffer<save_type>::convert_to_save_type(IM_STORE_TYPE x){/*{{{*/
	save_type max=-1;
	return (save_type) (x*max);
}/*}}}*/

/*}}}*/

class CImageEffect:public CObject{/*{{{*/
   public:
	virtual int apply(CImage *im);
};/*}}}*/

class CPsf:public CImageEffect{/*{{{*/
   public:
      CPsf(IM_COORD_TYPE sizex, IM_COORD_TYPE sizey);
      ~CPsf();
      virtual int apply(CImage *im);
   protected:
      int convolve(CImage *im);
      double *psf_data;
      int calculate_fft_psf();

      int sizex;
      int sizey;
   
   private:

      fftw_complex *fft_psf_data;
      fftw_complex *image_fft_data;
      double *image_data;
      
};/*}}}*/

class CGaussianPsf:public CPsf{/*{{{*/
   public:
      CGaussianPsf(IM_COORD_TYPE sizex, IM_COORD_TYPE sizey, float sigma, float astigmatism_ratio, float astigmatism_angle);
};/*}}}*/

//Vibration definition structure
typedef struct {/*{{{*/
   unsigned int pixel_dwell_time;
   float min_frequency;
   float max_frequency;
   float max_amplitude;
   int number_of_frequencies;
   unsigned int pixel_dead_time;
   unsigned int line_dead_time;
} t_vib_definition;/*}}}*/

class CVibration:public CImageEffect{/*{{{*/
   public:
      CVibration(t_vib_definition *def);
      CVibration(t_vib_definition *def, vector <float*> vibs);
      ~CVibration();
      virtual int apply(CImage *im, unsigned long time_shift);
   protected:
      void generate_curves(const int sizex, const int sizey, unsigned long time_shift);
      void destroy_curves();
   private:
      float *frequencies;     //[Hz]
      float *phases;     	//[rad]
      float *amplitudes_x;      //[nm]
      float *amplitudes_y;      //[nm]
      double *motion_x;
      double *motion_y;
      t_vib_definition def;
      void const_init();

};/*}}}*/

class CNoise:public CImageEffect{/*{{{*/
   public:
      CNoise();
      virtual int apply(CImage *im);
   protected:
      virtual double noise_value(double n);
};/*}}}*/

class CPoissonNoise:public CNoise{/*{{{*/
   public:
      CPoissonNoise(double particles_per_unit);
   protected:
      virtual double noise_value(double n);
      double particles_per_unit;
};/*}}}*/

class CGaussianNoise:public CNoise{/*{{{*/
   public:
      CGaussianNoise(double sigma);
   protected:
      virtual double noise_value(double n);
   private:
      double sigma;
};/*}}}*/

class CBackgroud:public CImageEffect{/*{{{*/
   public:
      CBackgroud();
      virtual int apply(CImage *im);
};/*}}}*/

class CEvenBackgroud:public CBackgroud{/*{{{*/
   public:
      CEvenBackgroud(IM_STORE_TYPE gl);
      virtual int apply(CImage *im);
      IM_STORE_TYPE give_bg_greylevel();
   private:
      IM_STORE_TYPE bg_greylevel;

};/*}}}*/

class CWavyBackgroud:public CBackgroud{/*{{{*/
   public:
      CWavyBackgroud(IM_STORE_TYPE min_gl, IM_STORE_TYPE max_gl, int densx, int densy);
      CWavyBackgroud(int densx, int densy, vector <IM_STORE_TYPE> values);
      ~CWavyBackgroud();
      virtual int apply(CImage *im);
   private:
      int densx, densy; // density values for both axes
      IM_STORE_TYPE min_gl, max_gl; // min and max gray levels
      CImage *bg_im;
      void const_init(int densx, int densy);
};/*}}}*/

}

#endif

// vim: cindent
