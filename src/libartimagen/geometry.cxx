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
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include "artimagen.h"
#include <cassert>
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
 

namespace artimagen {
   /* cross product of two vectors */
DIST_TYPE cross(CVector a, CVector b){/*{{{*/
   return(a.x * b.y - a.y * b.x);
}/*}}}*/

   /* dot product of two vectors */
DIST_TYPE dot(CVector a, CVector b){/*{{{*/
   return(a.x * b.x + a.y * b.y);
}/*}}}*/
}

using namespace artimagen;
////////////////// CVector //////////////////

/* empty CVector constructor */
CVector::CVector(){/*{{{*/
   x = 0;
   y = 0;
}/*}}}*/

/* CVector constructor with value initialization */
CVector::CVector(DIST_TYPE newx, DIST_TYPE newy){/*{{{*/
   x = newx;
   y = newy;
}/*}}}*/

CVector CVector::operator + (CVector a){/*{{{*/
	return CVector(x+a.x, y+a.y);
}/*}}}*/

CVector CVector::operator - (CVector a){/*{{{*/
	return CVector(x-a.x, y-a.y);
}/*}}}*/

CVector CVector::operator * (DIST_TYPE a){/*{{{*/
	return CVector(a*x, a*y);
}/*}}}*/

   CVector CVector::operator = (CVector a){/*{{{*/
   x = a.x;
   y = a.y;
   return *this;
}/*}}}*/

/* Squared length of the vector method */
DIST_TYPE CVector::sq_length(){/*{{{*/
	return x*x+y*y;
}/*}}}*/

/* Relation of linear dependence of two vectors */
int CVector::is_colinear(CVector a){/*{{{*/
   if (cross(*this,a) == 0) return 1;
   return 0;
}/*}}}*/


/* String-returning function for debugging */
string CVector::p(){/*{{{*/
	ostringstream s;
   	s << "V( " << x << "," << y << " )";
	return s.str();
}/*}}}*/

/* Vector-painting function for debugging 
 * Draws a little cross */
void CVector::paint(CImage *im){/*{{{*/
   	CLine l1(CVector(x-3,y-3),CVector(x+3,y+3));
   	CLine l2(CVector(x-3,y+3),CVector(x+3,y-3));
	l1.paint(im);
	l2.paint(im);
}/*}}}*/

/* Vector rotation (around the center) function */
void CVector::rotate(CVector center, double angle){/*{{{*/
   CVector v = *this - center;
   CVector vv;
   vv.x = v.x*cos(angle) + v.y*sin(angle);
   vv.y = v.y*cos(angle) - v.x*sin(angle);
   vv = vv + center;
   x = vv.x;
   y = vv.y;
}/*}}}*/

/////////////////// CLine ///////////////////

/* Empty line constructor, the p1=(1,1) so the length is >0 */
CLine::CLine(){/*{{{*/
   p0 = CVector(0,0);
   p1 = CVector(1,1);
}/*}}}*/

/* Constructor of a line from two vectors */
CLine::CLine(CVector p0, CVector p1){/*{{{*/
   this->p0 = p0;
   this->p1 = p1;
   if ((p0.x == p1.x) && (p0.y == p1.y)) {
      throw (int) AIG_EX_ZERO_LINE; // zero-length lines are illegal
   }
}/*}}}*/

/* Line assignment */
CLine CLine::operator = (CLine a){/*{{{*/
   p0 = a.p0;
   p1 = a.p1;
   return *this;
}/*}}}*/

/* Squared-lenght function */
DIST_TYPE CLine::sq_length(){/*{{{*/
   CVector d = p1 - p0; 
   return sqrt(d.sq_length());
}/*}}}*/

/* Line parallelism relation */
bool CLine::operator || (CLine a){/*{{{*/
   CVector d = (p1 - p0); // direction vector of this 
   CVector ad = (a.p1 - a.p0); // direction vector of a 
   
   if (d.is_colinear(ad)){ // are the dir vectors linearly dependent?
      return true; // if so, the lines are parallel
   }
   return false;
}/*}}}*/

/* Line overlap relation */
bool CLine::operator ^ (CLine a){/*{{{*/
   CVector d = (p1 - p0);

   if (*this || a){ // for overlap, lines must be parallel
      if (d.is_colinear(a.p0-p0)) return true; // and there starting-point vectors must be colinear
   }
   return false;
}/*}}}*/

/* Line intersections relation*/
double CLine::operator % (CLine a){/*{{{*/
   double t, u;
   if (intersection(a, &t, &u) < 0) return -1; // error, bad dir vectors
   
   if ((u<0) || (u>1) || (t<0)) return 0; // the intersection is in negative t or u outside segment
   return 1; // all OK and the intersection is within the tested segment.
}/*}}}*/

/* Line intersection function */
int CLine::intersection(CLine a, double *tt, double *uu){/*{{{*/
   if (*this || a){ // lines must not be parallel
      return -1; 
   }
   CVector A=p1-p0; // dir vector of this
   CVector B=p0;
   CVector C=a.p1-a.p0; // dir vector of a
   CVector D=a.p0;
	
   if (((C.x == 0) && (C.y == 0)) || ((A.x == 0) && (A.y == 0))){
      return 0; // this should never happen, because zero lines are illegal
   }

   double t = (C.y * (B.x - D.x) - C.x * (B.y - D.y)) / (C.x * A.y - C.y * A.x); // X = At + B, X is the intersection
   double u;
   if (abs(C.x) > abs(C.y)) u = (A.x * t + B.x - D.x) / C.x; else u = (A.y * t + B.y - D.y) / C.y; // X = Cu + D

   // if tt, or uu are NULL, that means, that they are not needed and should not
   // be set.
   if (tt) *tt=t;  // if tt is not NULL, set the value
   if (uu) *uu=u;  // the same for uu

   return 1; // all OK
}/*}}}*/

/* debug print-string function */
string CLine::p(){/*{{{*/
	ostringstream s;
   	s << "L( " << p0.x << "," << p0.y << " ) -> " << "( " << p1.x << "," << p1.y << " )" << endl;
	return s.str();
}/*}}}*/

/* debugging paint function */
void CLine::paint(CImage *im){/*{{{*/
	double step=0.7/sq_length();
	CVector a = p1 - p0;
	for (double t=0; t<1; t+=step){
	   CVector pix = p0 + a*t;
	   im->dot(pix.x, pix.y, 1);
	}
}/*}}}*/

/* parametrize line segment. 0 <= t <= 1 */
CVector CLine::parametrize(double t){/*{{{*/
   return p0+(p1-p0)*t;
}/*}}}*/

/* Squared distance from a point (vector) to a line segment.
 * it's the perpendicular distance, or the distance to the closest 
 * end-point, whatever is valid */
DIST_TYPE CLine::sq_distance_from(CVector v){/*{{{*/
   CVector a = p1 - p0; // direction vector of this line
   CVector b = v - p0; // direction vector of P0 -> V

   DIST_TYPE t = dot(a,b)/a.sq_length(); // projection of b to a
   if (t>1) t=1; // if it's beyond of the end-points, use the nearest end-point.
   if (t<0) t=0;

   CVector x = p0 + a*t; // x = the projection point
   return (v-x).sq_length(); // return squared length of x -> v
}/*}}}*/

/////////////////// CTriangle ///////////////////

/* not yet used */
CTriangle::CTriangle(CVector va, CVector vb, CVector vc){/*{{{*/
   this->va = va;
   if (cross(vb-va, vc-va) > 0){ // A, B, C must be counter-clockwise for further processing
      this->vb = vb;
      this->vc = vc;
   } else {
      this->vc = vb;
      this->vb = vc;
   }

   ea = CLine(vb, vc);
   eb = CLine(vc, va);
   ec = CLine(va, vb);
}/*}}}*/

/////////////////// CCurve ///////////////////
/* CCurve is a general curve class; all specific curve-classes are
 * dreived from it. */

CCurve::CCurve(){/*{{{*/
   number_of_vertices = 0; 
   vertices = NULL;
}/*}}}*/

/* initialization function; same for all curves */
void CCurve::init(){
   approximate_by_vertices();
   calculate_bounding_box();
}

CCurve::~CCurve(){/*{{{*/
   if (vertices) delete [] vertices;
}/*}}}*/

CVector CCurve::parametrize(DIST_TYPE t){/*{{{*/
   return CVector(0,0);
}/*}}}*/

/* polygon-approximation function */
void CCurve::approximate_by_vertices(){/*{{{*/
   vertices = new CVector[number_of_vertices];
   for (int i=0; i<number_of_vertices; i++){
      DIST_TYPE t = (DIST_TYPE) i / (number_of_vertices-1); // vertices are placed equidistantly (parameter-t-wise)
      vertices[i] = parametrize(t);
   }
}/*}}}*/

/* minumum distance from a curve */
DIST_TYPE CCurve::sq_distance_from(CVector v){/*{{{*/
	int vert;
	DIST_TYPE min_dist=9e99; // some big number exceeding any possible distance 
	for (vert=1; vert < number_of_vertices; vert++){ //number approx lines = n.of vertices -1
		CLine l(vertices[vert-1],vertices[vert]);
		DIST_TYPE ll = l.sq_distance_from(v); // find the distance from segment
		if (min_dist > ll) min_dist = ll; // find minimum
	}
	return	min_dist; // return the minimum distance from all segments
}/*}}}*/

/* find the bounding box */
void CCurve::calculate_bounding_box(){/*{{{*/
   int i;
   CVector *tl = &bounding_box_tl; // tl = top-left corner
   CVector *br = &bounding_box_br; // br = bottom-right corner

   tl->x=vertices[0].x; // initialization
   tl->y=vertices[0].y;
   br->x=vertices[0].x;
   br->y=vertices[0].y;
   for (i=1; i<number_of_vertices; i++){
      //move tl or br if a vertice is outside the current bounding-box
      if (vertices[i].x < tl->x) tl->x = vertices[i].x; 
      if (vertices[i].y < tl->y) tl->y = vertices[i].y;
      if (vertices[i].x > br->x) br->x = vertices[i].x;
      if (vertices[i].y > br->y) br->y = vertices[i].y;
   }
}/*}}}*/

/* export bounding-box vectors */
void CCurve::give_bounding_box(CVector *tl, CVector *br){/*{{{*/
	*tl = bounding_box_tl;
	*br = bounding_box_br;
}/*}}}*/

int CCurve::give_number_of_vertices(){/*{{{*/
   return number_of_vertices;
}/*}}}*/

CVector *CCurve::give_vertices(){/*{{{*/
   return this->vertices;
}/*}}}*/

/* debugging paint function */
void CCurve::paint(CImage *im){/*{{{*/
  int i;
  for (i=0; i<number_of_vertices; i++) {
     im->dot(vertices[i].x, vertices[i].y,1);
  }
}/*}}}*/

/* detection of intersection of a curve-segment with a testline */
int CCurve::is_segment_hit(CVector p1, CVector p2, CVector x, CVector d){/*{{{*/
   CLine testline(x,d);
   CLine segment(p1,p2);

   if (testline || segment) { // are lines parallel?
      if (testline ^ segment) return CU_SEGMENT_HIT_ERR; // if parallel and overlapping, return ERR
      return CU_SEGMENT_NOT_HIT; // if parallel but no overlap, then NOT_HIT
   }

   double u = 0;
   double t = 0;
   testline.intersection(segment, &t, &u);

   // testline must be hit in positive values of t (from x thru x+d and further
   // but not from x to the other side:
   // TESTLINE:
   // ------------[X]---------->[X+D]----------[infinity]
   //             t=0            t=1 
   //   NOT HIT    |      HIT
   //  
   //	Segment must be hit within the end-points. if close to the end-points,
   //	HIT_VERTEX is returned
   //
   //  SEGMENT:
   //  -------------------------[P0]----------------------[P1]---------------
   //          NOT_HIT      HIT_VERTEX        HIT     HIT_VERTEX   NOT_HIT

   if ((t > 0) && (u<=1) && (u>=0)){ 
      // if the intersection is close to any of the endpoints, the result may not be valid.
      if ( ((u < 0.0001) && (u > -0.00001 )) || ((u > 0.99999) && (u < 1.000001)) ) {  
	 return CU_SEGMENT_HIT_VERTEX; 

      }
      return CU_SEGMENT_HIT;
   }
   return CU_SEGMENT_NOT_HIT;
}/*}}}*/


/* retunrs number of intersection of the testline with the curve or error, 
 * if the testline is not suitable */
int CCurve::is_hit(CVector x, CVector d){/*{{{*/
   int i;
   int hit, hits;
   hits = 0;

   hit = CU_SEGMENT_NOT_HIT;
   i = 1;
   try{
      while (i<number_of_vertices){ 
	 hit = is_segment_hit(vertices[i-1],vertices[i],x,d);
	 if (hit == CU_SEGMENT_HIT_ERR) throw -1; 	// error output - loop broken
	 if (hit == CU_SEGMENT_HIT_VERTEX) throw -2; 	// test line hits a vertex -
	    							// unusable for inside/outside detection.
	 if (hit == CU_SEGMENT_HIT) hits++;
	 i++;
      }
   }
   catch (int err){
      return err;
   }
return hits;
}/*}}}*/

/* detect intersections of the curve with a line segment
 * useful for construction of grid maps of features */
int CCurve::hits_line(CLine l){/*{{{*/
   int hits = 0;
   for (int i=1; i<number_of_vertices; i++){
      CLine l2(vertices[i-1], vertices[i]);

      if ((l || l2) && ( l ^ l2)) { // lines are parallel and are overlapping.
	 return -1; //  this means ERROR
      }

      double t;
      if (l.intersection(l2, &t, NULL) <=0) return -1; 
      if ( (t >=0) && (t <=1)) hits++;
   }
   return hits;
}/*}}}*/

/* move all vertices using the shift-vector mv */
void CCurve::move_by(CVector mv){/*{{{*/
   for (int i=0; i< number_of_vertices; i++){
      vertices[i] = vertices[i] + mv;
   }
   calculate_bounding_box();
}/*}}}*/

/* rotate all vertices */
void CCurve::rotate(CVector center, double angle){/*{{{*/
   for (int i=0; i< number_of_vertices; i++){
      vertices[i].rotate(center, angle);
   }
   calculate_bounding_box();
}/*}}}*/

/////////////////// CStraightLine //////////////////
/* Straight line is a simplest curve */

CStraightLine::CStraightLine(CLine l){/*{{{*/
   line = l;
   number_of_vertices = 3; // whatever, 3 is OK for subdivision of a straight line.
   init();
}/*}}}*/

CVector CStraightLine::parametrize(DIST_TYPE t){/*{{{*/
   return line.parametrize(t); // easy peasy
}/*}}}*/

/////////////////// CBezier //////////////////

CBezier::CBezier(CVector p1, CVector p2, CVector p3, CVector p4)/*{{{*/
{
   // just set the control points and run init.
   control_points[0] = p1;
   control_points[1] = p2;
   control_points[2] = p3;
   control_points[3] = p4;
   number_of_vertices = 5; // how many linear segments approximate this curve
   init();
}/*}}}*/

/* Bezier Curve parametrization - see e.g. Mathworld for the definition*/
CVector CBezier::parametrize(DIST_TYPE t){/*{{{*/
   return control_points[0]*(1-t)*(1-t)*(1-t) +
   control_points[1]*3*t*(1-t)*(1-t) +
   control_points[2]*3*t*t*(1-t) +
   control_points[3]*t*t*t;
}/*}}}*/

//////////////////// CPolygon /////////////////
/* CPolygon is derived from CCurve. All Features are in turn
 * converted to polygons. The curve is approximated by linear
 * segments, so polygon is the approximation of a closed curve*/

CPolygon::CPolygon(){/*{{{*/
   map_division = 5; // just a guess for now. :-)
   map = NULL;
}/*}}}*/

/* Detection intersection of the polygon  with a testline */
int CPolygon::is_hit(CVector x, CVector d){/*{{{*/
   int hits = 0;
   do {
      hits = CCurve::is_hit(x, d);
      if (hits < 0) {
	 d = x + CVector((float)rand()*10/RAND_MAX,(float)rand()*10/RAND_MAX); // in case of error, make a random change
      }
   } while (hits < 0);
   return hits;
}/*}}}*/

/* Decection if the point p is inside the polygon using testline and
 * intersections */
int CPolygon::is_inside_with_hits(CVector p){/*{{{*/
   CVector d = p + CVector(121,111); // this is some ugly vector that should not collide with the vertices
   if (is_hit(p, d) % 2 == 1) return 1;
   return 0;
}/*}}}*/

/* Is the piont p inside the polygon? Used map first, then hits */
int CPolygon::is_inside(CVector p){/*{{{*/
   CVector bbtl, bbbr;
   int i,j;

   give_bounding_box(&bbtl, &bbbr);

   // Test 1: is p inside the bounding box?
   if ((p.x < bbtl.x) || (p.y < bbtl.y) || (p.x >= bbbr.x) || (p.y >= bbbr.y)) return 0; // p is outside the boundign box.


   // Test 2: map test
   //
   // this is the map test (speeds up the detection significantly). 2 = whole segment is
   // inside, 0 means definitely outside.
   // Find the map indexes corresponding to p:
   i = floor ((p.x - bbtl.x) * map_division / (bbbr.x - bbtl.x));
   j = floor ((p.y - bbtl.y) * map_division / (bbbr.y - bbtl.y));
   if ((j * map_division + i) >= sq(map_division)) {
      cout << "map_division = " << map_division << ", i=" << i << "; j=" << j <<endl;
   }
   if (map[j * map_division + i] == 2) return 1; // the whole map segment is inside, thus p is inside.
   if (map[j * map_division + i] == 0) return 0; //  segment outside, thus p outside

   // Test 3: the low-level test using testlines and intersections.
   if (is_inside_with_hits(p)) return 1; 
   return 0;
}/*}}}*/

/* detection of overlap with other polygon fe */
bool CPolygon::overlaps(CPolygon *fe){/*{{{*/
   int overlaps = 0;
   for (int j=0; j< number_of_vertices; j++){
      // if any fe's vertex is inside this polygon, we've got the overlap
      if (fe->is_inside(give_vertices()[j])) overlaps = 1; 
   }
   if (overlaps) return true;
   return false;

}/*}}}*/

/* the mapping function */
void CPolygon::create_map(){/*{{{*/
// map values: 
// 	0 segment outside polygon
// 	1 uncertain
// 	2 segment inside polygon

   map = new char[sq(map_division)];
   // zero it
   for (unsigned int i=0; i < sq(map_division); i++) map[i] = 0;


   // using lines - sides of the map segments find those segments
   // that intersect with the polygon. the map segments sharing that
   // side are uncertain.
   //
   // Horizontal lines:
   for (unsigned int j=0; j<= map_division; j++)  
      // horizontal lines tested for intersections with the feature's segments.
      // n+1 lines, n columns.
      for (unsigned int i=0; i < map_division; i++){
	 DIST_TYPE x,y,xp;
	 x = bounding_box_tl.x + (bounding_box_br.x - bounding_box_tl.x)*i/(map_division);
	 y = bounding_box_tl.y + (bounding_box_br.y - bounding_box_tl.y)*j/(map_division);
	 xp = x + (bounding_box_br.x - bounding_box_tl.x)/(map_division);// xp is x(i+1)

	 CLine testline(CVector(x,y), CVector(xp,y)); // create the test line.
	 int hits = hits_line(testline);
	 if (hits != 0) { 
	    // we have a hit or overlap => make map segments neighboring with
	    // the test line uncertain
	    if (j>0) map[index(i,j-1,map_division)] = 1; 
	    if (j<map_division) map[index(i,j,map_division)] = 1;
	    // if we are within the image, make the element uncertain
	 }
      }

   // Vertical lines:
   for (unsigned int j=0; j< map_division; j++) 
	 // vertical lines tested for intersections with the feature's segments.
	 // n lines, n+1 columns.
      for (unsigned int i=0; i <= map_division; i++){
	 DIST_TYPE x,y,yp;
	 x = bounding_box_tl.x + (bounding_box_br.x - bounding_box_tl.x)*i/(map_division);
	 y = bounding_box_tl.y + (bounding_box_br.y - bounding_box_tl.y)*j/(map_division);
	 yp = y + (bounding_box_br.y - bounding_box_tl.y)/(map_division);// yp is y(i+1)

	 CLine testline(CVector(x,y), CVector(x,yp)); // create the test line.
	 if (hits_line(testline) != 0) 
	 { 
	    // we have a hit or overlap => make segments neighboring with
	    // the test line uncertain
	    if (i>0) map[index(i-1,j,map_division)] = 1; 
	    if (i<map_division) map[index(i,j,map_division)] = 1; 
	    // if we are within the image, make the element uncertain
	 }
      }
   //  mapping which 0 elements are inside and which outside.

   // The map segments set to 0 are now certainly inside or outside the polygon.
   // So inside or outside? The following decides.

   for (unsigned int j=0; j< map_division; j++) 
      for (unsigned int i=0; i < map_division; i++){
	 // x,y are in the middle of the segment, unlike before.
	 DIST_TYPE x = bounding_box_tl.x + (bounding_box_br.x - bounding_box_tl.x)*((float) i+0.5)/(map_division);
	 DIST_TYPE y = bounding_box_tl.y + (bounding_box_br.y - bounding_box_tl.y)*((float) j+0.5)/(map_division);
	 if (map[index(i,j,map_division)] == 0) { // the segment is not uncertain
	    if (is_inside_with_hits(CVector(x,y))){ // find out if inside
	       map[index(i,j,map_division)] = 2; // yep, we're inside
	    }
	}
      }

}/*}}}*/

void CPolygon::set_map_division(int md){/*{{{*/
   map_division = md;
}/*}}}*/

CPolygon::~CPolygon(){/*{{{*/
   if (map) delete [] map;
}/*}}}*/


//////////////////// CVoronoi /////////////////
// NOT YET USED

CVoronoi::CVoronoi(CPolygon *polygon, CVector *free_points, int num_of_free){/*{{{*/
   this->polygon = polygon;
   polygon->give_bounding_box(&bounding_box_tl, &bounding_box_br);
   this->free_points = free_points;
   number_of_free = num_of_free;
}/*}}}*/

CVoronoi::~CVoronoi(){/*{{{*/
   delete [] map;
}/*}}}*/

void CVoronoi::calculate_map(int hsplit, int vsplit){/*{{{*/
   this->hsplit = hsplit;
   this->vsplit = vsplit;

   if (map) delete [] map;
   map = new int[hsplit * vsplit];

   for (int j=0; j<vsplit; j++)
      for (int i=0; i<hsplit; i++){
	 CVector t = map2vec(i,j);
	 CVector v = t - free_points[0];
	 DIST_TYPE d = v.sq_length();
	 for (int q=1; q<number_of_free; q++) {
	    v = t - free_points[q];
	    DIST_TYPE dd = v.sq_length();
	    if (dd < d) { 
	       d = dd;
	       map[index(i,j,hsplit)] = q;
	    }
	 }
	 try {
	    for (int q=0; q < polygon->give_number_of_vertices(); q++) {
	       v = t - polygon->give_vertices()[q];
	       DIST_TYPE dd = v.sq_length();
	       if (dd < d) { 
		  map[index(i,j,hsplit)] = -1;
		  throw -1; // this point is closest to the bordering point - no mor porcessing needed.
	       }
	    }
	 }
	 catch (int err){
	 }

      }
}/*}}}*/

void CVoronoi::lloyd(){/*{{{*/
   CVector *new_free = new CVector[number_of_free];
   int *counts = new int[number_of_free];
   memset (new_free,0, sizeof(new_free));
   memset (counts,0, sizeof(counts));

   for (int j=0; j<vsplit; j++)
      for (int i=0; i<hsplit; i++){
	 int free_point_index = map[index(i,j,hsplit)];
	 if (free_point_index >-1){
	    new_free[free_point_index] = new_free[free_point_index] + map2vec(i,j);
	    counts[free_point_index]++;
	 }
      }
   for (int i=0; i<number_of_free; i++) {
      CVector new_point = new_free[i] * ((float) 1 / counts[i]);
      if (polygon->is_inside(new_point)) free_points[i] = new_point;
   }
   delete [] counts;
   delete [] new_free;

}/*}}}*/

CVector CVoronoi::map2vec(int x, int y){/*{{{*/
   DIST_TYPE xo= bounding_box_tl.x + x*(bounding_box_br.x - bounding_box_tl.x)/hsplit;
   DIST_TYPE yo= bounding_box_tl.y + y*(bounding_box_br.y - bounding_box_tl.y)/vsplit;

   return CVector(xo,yo);
}/*}}}*/

// vim: cindent
