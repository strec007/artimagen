/* 
    Artificial SEM Image Generator (ArtImaGen)
    2009  Petr Cizmar @ National Institute of Standards and Technology
    E-mail: petr.cizmar@nist.gov

    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to Gnu.org public domain is compatible with GPL.

 */

using namespace std;
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include "artimagen_i.h"
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
 
using namespace artimagen;

////////////////////CFeature /////////////////
/* This is the general feature class, all specific classes 
 * are derived from this one */

CFeature::CFeature(){/*{{{*/
   sender_id = "CFeature";
   number_of_curves = 0;
   curves = NULL;
   number_of_effects = 0;
   effects = NULL;
   base_gray_level = 0.6; // default value
}/*}}}*/

// global initialization function
void CFeature::init(){/*{{{*/
   build_vertices();
   calculate_bounding_box();

   CVector bbtl, bbbr;
   give_bounding_box(&bbtl, &bbbr);
   set_map_division(1 + (int)((bbbr.x - bbtl.x)/4)); // 4 pixels per map segment
   create_map();
}/*}}}*/

/* copy vertices from included curves */
void CFeature::build_vertices(){/*{{{*/
   // count all vertices of all included curves.
   number_of_vertices = 0;
   // the last vertex of each curve is cut out, to prevent duplicate 
   // vertices
   for (int i=0; i< number_of_curves; i++) this->number_of_vertices += curves[i]->give_number_of_vertices() -1; 
   number_of_vertices++; // Last point must be added to close the polygon.
   // copy the vertices from curves to the local list of vertices.
   vertices = new CVector[this->number_of_vertices];
   int counter = 0;
   for (int i=0; i< number_of_curves; i++) 
      for (int j=0; j < curves[i]->give_number_of_vertices()-1; j++){
	 vertices[counter++] = curves[i]->give_vertices()[j];
      }
   vertices[number_of_vertices-1] = vertices[0]; // ...closing the polygon.
}/*}}}*/

/* paint feature into an image */
void CFeature::paint(CImage *im){/*{{{*/
   int i;
   IM_COORD_TYPE x, y;
   CVector bbtl, bbbr; // bounding box
   IM_STORE_TYPE val;

   give_bounding_box(&bbtl, &bbbr);

   for (y=bbtl.y; y<=bbbr.y; y++) // go thru the bounding box
      for (x=bbtl.x; x<=bbbr.x; x++){
	 CVector vec(x,y);
	 if (is_inside(vec)) { // if the vec is outside the feature, no paint
	    val = base_gray_level;
	    for (i=0; i<number_of_effects; i++){
	       CVector pixel(x,y);
	       // feature effects are mutliplicative, they return 
	       // an amplification coefficient here applied.
	       val *= effects[i]->give_amplification(this, pixel);
	    }
	    im->dot(x,y,val); // paint the pixel
	 }
      }

}/*}}}*/

IM_STORE_TYPE CFeature::give_base_gray_level(){/*{{{*/
   return base_gray_level;
}/*}}}*/

void CFeature::set_base_gray_level(IM_STORE_TYPE level){/*{{{*/
   base_gray_level = level;
}/*}}}*/

void CFeature::add_effect_chain(CEffect **effects, int number_of_effects){/*{{{*/
   this->effects = effects;
   this->number_of_effects = number_of_effects;
}/*}}}*/

CFeature::~CFeature(){/*{{{*/
   if (curves) {
      for (int i=0; i < number_of_curves; i++){
	 if (curves[i]) delete curves[i];
	 curves[i] = NULL;
      }
      delete [] curves;
      curves = NULL;
   }
}/*}}}*/

//////////////////// CGoldenGrain /////////////////

CGoldenGrain::CGoldenGrain(const CVector position, const DIST_TYPE size){/*{{{*/
   sender_id = "CGoldenGrain";
   const int NOC = 7;
   number_of_curves = NOC;
   curves = new CCurve*[number_of_curves];
   for (int i=0; i<number_of_curves; i++) curves[i]=NULL; // important if curve-constructor throws
   CVector P[NOC]; // end control point
   CVector L[NOC]; // lever control point

   const DIST_TYPE displ_rand = size/3; // pixel displacement randomization

   const float l_radius_factor = 0.2;

   float phi = rand()*2*M_PI/RAND_MAX;

   for (int i=0; i<number_of_curves; i++){ // distribution of basic control points
      // P-points are distributed evenly along a circle + ramdomly displaced
      // (displ_rand)
      P[i].x = position.x + size * cos(phi+i*2*M_PI/number_of_curves) + rand()*displ_rand/RAND_MAX;
      P[i].y = position.y + size * sin(phi+i*2*M_PI/number_of_curves) + rand()*displ_rand/RAND_MAX;

   }

   for (int i=0; i<number_of_curves; i++){ // distribution of "lever" control points
      int im = i-1;
      if (im < 0) im += number_of_curves;
      int ip = i+1;
      if (ip >= number_of_curves) ip -= number_of_curves;
      // the lever L-point is placed has the same direction as
      // P[i+1]-P[i-1] (previous and next P points),
      // distance from the current P poibt is determined by l_radius_factor.
      L[i] = P[i]+(P[im]+P[ip]*-1)*-l_radius_factor;
   }

   for (int i=0; i<number_of_curves; i++){  // creating the Beziers.
      int ip = i+1;
      if (ip >= number_of_curves) ip -= number_of_curves;

      // the second lever point (third control point) must be oposite the 
      // next L-point... due to continuity 
      //
      curves[i] = new CBezier(P[i],L[i],P[ip]*2+L[ip]*-1,P[ip]); 
   }
   init();

}/*}}}*/

///////////////////// CRectangle ////////////////

CRectangle::CRectangle(const DIST_TYPE lsize, const DIST_TYPE tsize, double rotation){/*{{{*/
   sender_id = "CRectangle";
   number_of_curves = 4;
   curves = new CCurve*[number_of_curves];
   for (int i=0; i<number_of_curves; i++) curves[i]=NULL; // important if curve-constructor throws

   
   curves[0] = new CStraightLine(CLine(CVector(0,0), CVector(tsize, 0)));
   curves[1] = new CStraightLine(CLine(CVector(tsize,0), CVector(tsize, lsize)));
   curves[2] = new CStraightLine(CLine(CVector(tsize,lsize) , CVector(0, lsize)));
   curves[3] = new CStraightLine(CLine(CVector(0,lsize), CVector(0,0)));

   for (int i=0; i< number_of_curves; i++) curves[i]->rotate(CVector(0,0), rotation);

   init();
}/*}}}*/

///////////////////// CSnake ////////////////

CSnake::CSnake(const DIST_TYPE w, const DIST_TYPE a, const DIST_TYPE b, const DIST_TYPE c, double rotation){/*{{{*/
   sender_id = "CSnake";
   number_of_curves = 8;
   curves = new CCurve*[number_of_curves];
   for (int i=0; i<number_of_curves; i++) curves[i]=NULL; // important if curve-constructor throws

   curves[0] = new CStraightLine(CLine(CVector(0, 0), CVector(w, 0)));
   curves[1] = new CStraightLine(CLine(CVector(w, 0), CVector(w, a-w)));
   curves[2] = new CStraightLine(CLine(CVector(w, a-w), CVector(b+w, a-w)));
   curves[3] = new CStraightLine(CLine(CVector(b+w, a-w), CVector(b+w, a+c-w)));
   curves[4] = new CStraightLine(CLine(CVector(b+w, a+c-w), CVector(b, a+c-w)));
   curves[5] = new CStraightLine(CLine(CVector(b, a+c-w), CVector(b, a)));
   curves[6] = new CStraightLine(CLine(CVector(b, a), CVector(0, a)));
   curves[7] = new CStraightLine(CLine(CVector(0, a), CVector(0, 0)));

   for (int i=0; i< number_of_curves; i++) curves[i]->rotate(CVector(0,0), rotation);

   init();
}/*}}}*/

///////////////////// CCorner ////////////////

CCorner::CCorner(const DIST_TYPE lsize, const DIST_TYPE tsize, double rotation){/*{{{*/
   sender_id = "CCorner";
   number_of_curves = 6;
   curves = new CCurve*[number_of_curves];
   for (int i=0; i<number_of_curves; i++) curves[i]=NULL; // important if curve-constructor throws

   
   curves[0] = new CStraightLine(CLine(CVector(0,0), CVector(tsize, 0)));
   curves[1] = new CStraightLine(CLine(CVector(tsize,0), CVector(tsize, lsize)));
   curves[2] = new CStraightLine(CLine(CVector(tsize,lsize) , CVector(lsize, lsize)));
   curves[3] = new CStraightLine(CLine(CVector(lsize,lsize) , CVector(lsize, tsize)));
   curves[4] = new CStraightLine(CLine(CVector(lsize,tsize) , CVector(0, tsize)));
   curves[5] = new CStraightLine(CLine(CVector(0,tsize), CVector(0,0)));

   init();
}/*}}}*/

///////////////////// CCross ////////////////

CCross::CCross(const DIST_TYPE lsize, const DIST_TYPE tsize, double rotation){/*{{{*/
   sender_id = "CCross";
   number_of_curves = 12;
   curves = new CCurve*[number_of_curves];
   for (int i=0; i<number_of_curves; i++) curves[i]=NULL; // important if curve-constructor throws

   
   curves[0] = new CStraightLine(CLine(CVector(tsize,-lsize), CVector(tsize,-tsize)));
   curves[1] = new CStraightLine(CLine(CVector(tsize,-tsize), CVector(lsize,-tsize)));
   curves[2] = new CStraightLine(CLine(CVector(lsize,-tsize), CVector(lsize,tsize)));
   curves[3] = new CStraightLine(CLine(CVector(lsize,tsize), CVector(tsize,tsize)));
   curves[4] = new CStraightLine(CLine(CVector(tsize,tsize), CVector(tsize,lsize)));
   curves[5] = new CStraightLine(CLine(CVector(tsize,lsize), CVector(-tsize,lsize)));
   curves[6] = new CStraightLine(CLine(CVector(-tsize,lsize), CVector(-tsize,tsize)));
   curves[7] = new CStraightLine(CLine(CVector(-tsize,tsize), CVector(-lsize,tsize)));
   curves[8] = new CStraightLine(CLine(CVector(-lsize,tsize), CVector(-lsize,-tsize)));
   curves[9] = new CStraightLine(CLine(CVector(-lsize,-tsize), CVector(-tsize,-tsize)));
   curves[10] = new CStraightLine(CLine(CVector(-tsize,-tsize), CVector(-tsize,-lsize)));
   curves[11] = new CStraightLine(CLine(CVector(-tsize,-lsize), CVector(tsize,-lsize)));

   for (int i=0; i< number_of_curves; i++) curves[i]->rotate(CVector(0,0), rotation);
   for (int i=0; i< number_of_curves; i++) curves[i]->move_by(CVector(lsize,lsize));

   init();
}/*}}}*/

//////////////////// CEffect /////////////////
/* general feature-effect class */

CEffect::CEffect(){/*{{{*/
   sender_id = "CEffect";
}/*}}}*/

float CEffect::give_amplification(CFeature *fe, CVector v){/*{{{*/
   return fun(fe, v);
}/*}}}*/


//////////////////// CEdgeEffect /////////////////

CEdgeEffect::CEdgeEffect(float coefficient, IM_STORE_TYPE top_edge_value_above_base, DIST_TYPE thickness):CEffect(){/*{{{*/
   sender_id = "CEdgeEffect";
   this->coefficient = coefficient;
   this->top_edge_value_above_base = top_edge_value_above_base;
   this->thickness = thickness;
}/*}}}*/

double CEdgeEffect::fun(CFeature *fe, CVector v){/*{{{*/
   double a, b, c;

   if (!fe->is_inside(v)) return 1;

   center_feature_value = fe->give_base_gray_level();

   b = coefficient;
   a = (double) top_edge_value_above_base;
   c = (double) center_feature_value;

   DIST_TYPE dist = sqrt(fe->sq_distance_from(v));
   return (a * exp(-b * dist) + c) / (double) center_feature_value;
}/*}}}*/


//////////////////// CFineStructureEffect /////////////////
/* makes little round spots on top of the features */

CFineStructureEffect::CFineStructureEffect(float density, DIST_TYPE min_r, DIST_TYPE max_r, float min_coe, float max_coe):CEffect(){/*{{{*/
   sender_id = "CFineStructureEffect";
   this->density = density;
   this->min_r = min_r; //minimum radius
   this->max_r = max_r; // maximum r
   this->min_coe = min_coe; // minimum aplification coefficient
   this->max_coe = max_coe; // max...
   spots = NULL;
}/*}}}*/

CFineStructureEffect::~CFineStructureEffect(){/*{{{*/
   if (spots) delete [] spots;
   spots = NULL; // unnecessary, but this may change later. :-)
}/*}}}*/

void CFineStructureEffect::generate_spots(DIST_TYPE sizex, DIST_TYPE sizey){/*{{{*/
   nos = sizex * sizey * density;
   if (!spots) {
      // generation of spots  - random within limits
      spots = new finestruct_spot[nos];
      for (int i=0; i<nos; i++) {
	 spots[i].r = min_r+rand()*(max_r-min_r)/RAND_MAX;
	 spots[i].c = CVector(rand()*sizex/RAND_MAX, rand()*sizey/RAND_MAX);
	 spots[i].coe = min_coe + rand()*(max_coe - min_coe)/RAND_MAX;
      }
   }
}/*}}}*/

double CFineStructureEffect::fun(CFeature *fe, CVector v){/*{{{*/
   CVector tl, br;
   fe->give_bounding_box(&tl,&br);
   if (!spots) {
      sizex = br.x-tl.x;
      sizey = br.y-tl.y;
      generate_spots(sizex, sizey);
   }

   double fin_coe = 1; // final coefficient defaults to 1 - no amplification

   for (int i=0; i<nos; i++) {
      // if the vector v is inside the ith spot...
      if ((v-tl-spots[i].c).sq_length() < pow(spots[i].r,2)) {
	 // coefficient applies
	 fin_coe = spots[i].coe;
      }
   }
   return fin_coe;
}/*}}}*/

//////////////////// CSample /////////////////
// general sample class

CSample::CSample(DIST_TYPE sizex, DIST_TYPE sizey){/*{{{*/
   sender_id = "CSample";
   map_division = 11;
   features = NULL;
   number_of_features = 0;
   effects = NULL;
   number_of_effects = 0;
   actual_number_of_features = 0;
   map = NULL;
   map_counts = NULL;
   this->sizex = sizex;
   this->sizey = sizey;
   create_map();
}/*}}}*/

void CSample::add_feature_to_map(int i){/*{{{*/

   int cnt = map_division * map_division;
   int occ_map[cnt]; // allocate allocation map

   occupies_map(features[i], occ_map);

   for (int j = 0 ; j < map_division; j++)
      for (int k = 0; k < map_division; k++)
	 if(occ_map[j * map_division + k]) {
	    int cnt = map_counts[j * map_division + k];

	    unsigned int *orig_seg = map[j * map_division + k]; //get the original segment array
	    unsigned int *new_seg = new unsigned int[cnt+1]; // new segment will be one number bigger
	    for (int q = 0; q<cnt; q++) new_seg[q] = orig_seg[q]; //copy previous segment into the new one
	    new_seg[cnt] = i; // add the new feature index
	    map_counts[j * map_division + k]++;  // increase the number 
	    if (cnt) delete [] orig_seg; // delete only if previously defined
	    map[j * map_division + k] = new_seg; // place the new segment to the map
	 }
}/*}}}*/

void CSample::create_map(){/*{{{*/
   map = new unsigned int*[map_division * map_division];
   map_counts = new int[map_division * map_division];

   for (int i=0; i < map_division * map_division; i++) map_counts[i] = 0; // zero all map counts

}/*}}}*/

void CSample::occupies_map(CFeature *fe, int *map){/*{{{*/

   for (int i=0; i< map_division * map_division; i++) map[i] = 0;


   CVector bbtl, bbbr; // bounding box vectors
   fe->give_bounding_box(&bbtl, &bbbr);

   int miny = floor(bbtl.y /sizey * map_division);
   int maxy = ceil(bbbr.y /sizey * map_division);
   int minx = floor(bbtl.x /sizex * map_division);
   int maxx = ceil(bbbr.x /sizex * map_division);

   if (minx < 0) minx = 0; // must not get out of the map
   if (miny < 0) miny = 0;
   if (maxx > map_division) maxx = map_division;
   if (maxy > map_division) maxy = map_division;

   for (int j = miny; j < maxy; j++)
      for (int k = minx; k < maxx; k++){
	 map[j * map_division + k] = 1;
      }
}/*}}}*/

void CSample::paint(CImage *im){/*{{{*/
   for (int i=0; i<actual_number_of_features; i++){
      features[i]->paint(im);
   }
}/*}}}*/

void CSample::rotate(CVector center, double angle){/*{{{*/
   for (int i=0; i<actual_number_of_features; i++){
      features[i]->rotate(center, angle);
   }
}/*}}}*/

void CSample::destroy_map(){/*{{{*/
   if  ((map) && (map_counts)){ // map initialized and not yet destroyed
      for (int i=0; i < map_division * map_division; i++) if (map_counts[i] > 0) delete [] map[i]; // zero all map counts

      delete [] map;
      map = NULL;
      delete [] map_counts;
      map_counts = NULL;
   }
}/*}}}*/

int CSample::overlaps(CFeature *fe){/*{{{*/
   int cnt = map_division * map_division;
   int occ_map[cnt]; // allocate allocation map
   occupies_map(fe, occ_map); // fill in allocation for the ith feature

   for (int i=0; i<cnt; i++){
      if (occ_map[i]){
	 for(int k=0; k<map_counts[i]; k++){
	    if (map_counts[i] > 0) 
	    if (fe->overlaps(features[map[i][k]])) return SA_ADD_OVERLAP;
	 }
      }
   }

   return SA_ADD_OK;

}/*}}}*/

int CSample::add_feature(CFeature *fe){/*{{{*/
   if (number_of_features == actual_number_of_features) return SA_ADD_FULL;
   if (overlaps(fe) == SA_ADD_OVERLAP) return SA_ADD_OVERLAP;

   features[actual_number_of_features] = fe;

   add_feature_to_map(actual_number_of_features);

   actual_number_of_features++;
   return SA_ADD_OK;

}/*}}}*/

void CSample::move_by(CVector mv){/*{{{*/
   for (int i=0; i<actual_number_of_features; i++) ((CFeature *)features[i])->move_by(mv);
}/*}}}*/

CSample::~CSample(){/*{{{*/
   if (features) {
      for (int i=0; i<number_of_features; i++) if (features[i]) delete features[i];
      delete [] features;
      number_of_features = 0;
      features = NULL; // unnecessary - this is a destructor, but if it's moved...
   }
   if (effects) {
      for (int i=0; i<number_of_effects; i++) if (effects[i]) delete effects[i];
      delete [] effects;
      number_of_effects = 0;
      effects = NULL; // unnecessary - this is a destructor, but if it's moved... 
   }
   destroy_map();

}/*}}}*/

//////////////// CGoldOnCarbonSample /////////////////////

CGoldOnCarbonSample::CGoldOnCarbonSample(t_gc_definition *def):CSample(def->sizex, def->sizey){/*{{{*/

   sender_id = "CGoldOnCarbonSample";
   
   CGoldenGrain *gg; 

   number_of_features = def->number_of_grains;
   features = new CFeature*[number_of_features];
   actual_number_of_features = 0;
   for (int i=0; i<number_of_features; i++) features[i] = NULL; // important if feature constructor throws

    //cout << "Placing Gold Grains" << endl << "[";
   for (int i=0; i<number_of_features; i++){
      //cout << "*";
      //flush(cout);
      int adding_res;

      int size = def->grain_max_size-i*(def->grain_max_size - def->grain_min_size)/number_of_features;
      CVector pos(rand()*sizex/RAND_MAX, rand()*sizey/RAND_MAX); // make a random position
      gg = new CGoldenGrain(pos,size);
      gg->set_base_gray_level(def->base_level*(1+rand()*def->base_level_variation/RAND_MAX));


      effects = new CEffect*[2];
      number_of_effects = 2;
      for (int i=0; i<number_of_effects; i++) effects[i] = NULL; // important if effect constructor throws
      effects[0] = new CEdgeEffect(def->ee_coefficient, def->ee_top_above_base, def->ee_thickness);
      effects[1] = new CFineStructureEffect(def->fs_density, def->fs_min_r, def->fs_max_r, def->fs_min_coe, def->fs_max_coe);

      gg->add_effect_chain(effects, number_of_effects);
      int fails = 0; // placement-fail counter
      do {
	 adding_res = add_feature(gg);
	 if (adding_res == SA_ADD_OVERLAP) { // if the new structure overlaps, 
	    CVector npos(rand()*sizex/RAND_MAX, rand()*sizey/RAND_MAX); // make a new random position
	    gg->move_by(npos-pos); // move it there
	    pos = npos;
	    fails++;
	 }
	 if (fails > 1000) throw (int) AIG_EX_GOLDONCARBON_TOO_MANY_FAILS;
      } while (adding_res == SA_ADD_OVERLAP); // repeat until the placement succeeds
      if (adding_res == SA_ADD_FULL) {
	 cout << "Feature array full - You should never see this message" << cout;
      }
   }
   //IF_CO cout << "]" << endl;
}/*}}}*/

/////////////////// CPeriodicCornerSample /////////////////

CPeriodicCornerSample::CPeriodicCornerSample(t_cor_definition *def):CSample(def->sizex, def->sizey){/*{{{*/
   sender_id = "CPeriodicCornerSample";
   CCorner *cc;

   int countx = def->sizex / def->distance;
   int county = def->sizey / def->distance;

   number_of_features = countx * county;

   features = new CFeature*[number_of_features];
   for (int i=0; i<number_of_features; i++) features[i] = NULL; // important if feature constructor throws

   //IF_CO cout << "Placing Corners" << endl << "[";
   for (int j=0; j<county; j++)
      for (int i=0; i<countx; i++){
	 cc = new CCorner(def->lsize, def->tsize, 0);
	 cc->move_by(CVector(i*def->distance,j*def->distance));
	 cc->set_base_gray_level(def->base_level*(1+rand()*def->base_level_variation/RAND_MAX));

	 effects = new CEffect*[2];
	 number_of_effects = 2;
	 for (int i=0; i<number_of_effects; i++) effects[i] = NULL; // important if effect constructor throws
	 effects[0] = new CEdgeEffect(def->ee_coefficient, def->ee_top_above_base, def->ee_thickness);
	 effects[1] = new CFineStructureEffect(def->fs_density, def->fs_min_r, def->fs_max_r, def->fs_min_coe, def->fs_max_coe);

	 cc->add_effect_chain(effects,2);
	 if (add_feature(cc) == SA_ADD_OVERLAP) throw (int) AIG_EX_FEATURE_OVERLAP;
	 //IF_CO cout << "@";
      }
   //IF_CO cout << "]" << endl;
}/*}}}*/

/////////////////// CSingleRectangleSample /////////////////

CSingleRectangleSample::CSingleRectangleSample(t_rct_definition *def):CSample(def->sizex, def->sizey){/*{{{*/
   sender_id = "CSingleRectangleSample";
   CRectangle *cc;

   effects = new CEffect*[2];
   number_of_effects = 2;
   for (int i=0; i<number_of_effects; i++) effects[i] = NULL; // important if effect constructor throws
   effects[0] = new CEdgeEffect(def->ee_coefficient, def->ee_top_above_base, def->ee_thickness);
   effects[1] = new CFineStructureEffect(def->fs_density, def->fs_min_r, def->fs_max_r, def->fs_min_coe, def->fs_max_coe);

   number_of_features = 1;
   features = new CFeature*[number_of_features];
   for (int i=0; i<number_of_features; i++) features[i] = NULL; // important if feature constructor throws

   //IF_CO cout << "Placing Rectangle" << endl << "[";
   cc = new CRectangle(def->lsize, def->tsize, def->rotation);
   cc->move_by(CVector(def->sizex/2-def->tsize/2,def->sizey/2-def->lsize/2));
   cc->set_base_gray_level(def->base_level*(1+rand()*def->base_level_variation/RAND_MAX));
   cc->add_effect_chain(effects,2);
   if (add_feature(cc) == SA_ADD_OVERLAP) throw (int) AIG_EX_FEATURE_OVERLAP;
   //IF_CO cout << "@";
   //IF_CO cout << "]" << endl;
}/*}}}*/

/////////////////// CSnakeSample /////////////////

CSnakeSample::CSnakeSample(t_rct_definition *def):CSample(def->sizex, def->sizey){/*{{{*/
   sender_id = "CSnakeSample";

   effects = new CEffect*[1];
   number_of_effects = 1;
   for (int i=0; i<number_of_effects; i++) effects[i] = NULL; // important if effect constructor throws
   effects[0] = new CEdgeEffect(def->ee_coefficient, def->ee_top_above_base, def->ee_thickness);
   //effects[1] = new CFineStructureEffect(def->fs_density, def->fs_min_r, def->fs_max_r, def->fs_min_coe, def->fs_max_coe);

   number_of_features = 9;
   features = new CFeature*[number_of_features];
   for (int i=0; i<number_of_features; i++) features[i] = NULL; // important if feature constructor throws

   //IF_CO cout << "Placing Rectangle" << endl << "[";

   CFeature *fes[number_of_features];
   for (int i=0; i<number_of_features; i++) fes[i] = NULL; 

#define STEP 2*def->tsize

   fes[0] = new CRectangle(def->lsize, def->tsize, 0);
   fes[0] -> move_by(CVector(STEP,0));

   fes[1] = new CRectangle(def->lsize, def->tsize, 0);
   fes[1] -> move_by(CVector(3*STEP,0));

   fes[2] = new CSnake(def->tsize, (def->lsize+def->tsize)/2, 2*STEP, (def->lsize+def->tsize)/2, 0);
   fes[2] -> move_by(CVector(5*STEP,0));

   fes[3] = new CRectangle(def->lsize/2 - def->tsize/2 - def->tsize, def->tsize, 0);
   fes[3] -> move_by(CVector(5*STEP,def->lsize/2+1.5*def->tsize));

   fes[4] = new CRectangle(def->lsize/2 - 1.5*def->tsize, def->tsize, 0);
   fes[4] -> move_by(CVector(6*STEP,0));

   fes[5] = new CRectangle(def->lsize/2 - 1.5*def->tsize, def->tsize, 0);
   fes[5] -> move_by(CVector(7*STEP,0));

   fes[6] = new CRectangle(def->lsize, def->tsize, 0);
   fes[6] -> move_by(CVector(8*STEP,0));

   fes[7] = new CRectangle(def->lsize, def->tsize, 0);
   fes[7] -> move_by(CVector(9*STEP,0));

   fes[8] = new CRectangle(def->lsize, def->tsize, 0);
   fes[8] -> move_by(CVector(10*STEP,0));



   for (int i=0; i<number_of_features; i++) if (fes[i]){
      fes[i]->set_base_gray_level(def->base_level*(1+rand()*def->base_level_variation/RAND_MAX));
      fes[i]->add_effect_chain(effects,number_of_effects);
      if (add_feature(fes[i]) == SA_ADD_OVERLAP) {
	 cout << i << " overlaps!!!" << cout;
	 throw (int) AIG_EX_FEATURE_OVERLAP;
      }
   }
   rotate(CVector(6*STEP, def->lsize/2), def->rotation); 
   //IF_CO cout << "@";
   //IF_CO cout << "]" << endl;
}/*}}}*/

/////////////////// CPeriodicCrossSample /////////////////

CPeriodicCrossSample::CPeriodicCrossSample(t_crs_definition *def):CSample(def->sizex, def->sizey){/*{{{*/
   sender_id = "CPeriodicCrossSample";
   CCross *cc;

   int countx = def->sizex / def->distance;
   int county = def->sizey / def->distance;

   number_of_features = countx * county;

   features = new CFeature*[number_of_features];
   for (int i=0; i<number_of_features; i++) features[i] = NULL; // important if feature constructor throws

   //IF_CO cout << "Placing Crosses" << endl << "[";
   for (int j=0; j<county; j++)
      for (int i=0; i<countx; i++){
	 cc = new CCross(def->lsize, def->tsize, def->rotation);
	 cc->move_by(CVector(i*def->distance,j*def->distance));
	 cc->set_base_gray_level(def->base_level*(1+rand()*def->base_level_variation/RAND_MAX));

	 effects = new CEffect*[2];
	 number_of_effects = 2;
	 for (int i=0; i<number_of_effects; i++) effects[i] = NULL; // important if effect constructor throws
	 effects[0] = new CEdgeEffect(def->ee_coefficient, def->ee_top_above_base, def->ee_thickness);
	 effects[1] = new CFineStructureEffect(def->fs_density, def->fs_min_r, def->fs_max_r, def->fs_min_coe, def->fs_max_coe);

	 cc->add_effect_chain(effects,2);
	 if (add_feature(cc) == SA_ADD_OVERLAP) throw (int) AIG_EX_FEATURE_OVERLAP;
	 //IF_CO cout << "@";
      }
   //IF_CO cout << "]" << endl;
}/*}}}*/



// vim: cindent
