/* 
    Artificial SEM Image Generator (ArtImaGen)
    2009  Petr Cizmar @ National Institute of Standards and Technology
    E-mail: petr.cizmar@nist.gov

    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to Gnu.org public domain is compatible with GPL.

 */


#include "../../config.h"
#include "artimagen_i.h"
#include <vector>
#include <iostream>


#ifdef HAVE_LUA
extern "C" {
#include <lua5.1/lualib.h>
#include <lua5.1/lauxlib.h>
}
#endif

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

const DIST_TYPE overhead = 70;

extern "C" void *generate_gc_sample(t_gc_definition *def){/*{{{*/
   def->sizex += 2*overhead;
   def->sizey += 2*overhead;
   CSample *sam = new CGoldOnCarbonSample(def);
   def->sizex -= 2*overhead;
   def->sizey -= 2*overhead;
   return (void *) sam;
}/*}}}*/

extern "C" void destroy_gc_sample(void *poi){/*{{{*/
   CSample *sam = (CGoldOnCarbonSample *) poi;
   delete sam;
}/*}}}*/

extern "C" void *generate_standard_image(void *sample, t_std_image_def_struct *def){/*{{{*/

   DIST_TYPE o_sizex = def->sizex + 2*overhead;
   DIST_TYPE o_sizey = def->sizey + 2*overhead;
   CImage *im = new CImage(o_sizex, o_sizey);

   CWavyBackgroud back(def->bg_min_gl, def->bg_max_gl, def->bg_dens_x, def->bg_dens_y);
   back.apply(im);

   //generate_snake_structure_image(&sam);
   //generate_rectangle_structure_image(&sam);
   //generate_corner_structure_image(&sam);
   //generate_cross_structure_image(&sam);
   ((CSample *)sample)->paint(im);

   im->set_ifft_blocking(IM_IFFT_BLOCK);
   CGaussianPsf psf(o_sizex, o_sizey, def->beam_sigma, def->beam_astig_ratio , def->beam_astig_angle); //sizex, sizey, sigma, astig. ratio, ast. angle.
   psf.apply(im);
   im->shift(def->shift_x, def->shift_y);
   im->calculate_ifft();

   // create the drift-distortion function and apply it
   t_vib_definition vibdef;
   vibdef.pixel_dwell_time = def->vib_pixel_dwell_time;
   vibdef.min_frequency = def->vib_min_frequency;
   vibdef.max_frequency = def->vib_max_frequency;
   vibdef.max_amplitude = def->vib_max_amplitude;
   vibdef.number_of_frequencies = def->vib_number_of_frequencies;
   vibdef.pixel_dead_time = def->vib_pixel_dead_time;
   vibdef.line_dead_time = def->vib_line_dead_time;
   CVibration vib(&vibdef);
   vib.apply(im, 0);
   im->crop(overhead, overhead, o_sizex-overhead, o_sizey-overhead);
   CGaussianNoise gn(def->noise_sigma);
   gn.apply(im);
   return (void *) im;
}/*}}}*/

extern "C" void destroy_image(void *poi){/*{{{*/
   CImage *im = (CImage *) poi;
   delete im;
}/*}}}*/

extern "C" void save_image(void *image, char* filename, char *comment){/*{{{*/
   CImage *im = (CImage *) image;
   im->tiff_write(filename, comment, BI_8B);
}/*}}}*/

class CLuaMessenger : public CObject{/*{{{*/
   public:
      CLuaMessenger(const char *comment);
      CLuaMessenger(int msg_id, const char *comment);
};

CLuaMessenger::CLuaMessenger(const char *comment){
   sender_id = "LUA";
   send_message(AIG_MSG_LUA_ERROR, comment);
}/*}}}*/


CLuaMessenger::CLuaMessenger(int msg_id, const char *comment){ /*{{{*/
   sender_id = "LUA";
   send_message(msg_id, comment);
}/*}}}*/

#ifdef HAVE_LUA

void print_lua_stack(lua_State *L){
   int i;
   
   printf("\n ********* LUA STACK *********\n");

   for (i = lua_gettop(L); i>=-100; i--){
      int typ = lua_type(L,i);
      if (typ == LUA_TFUNCTION) break;
      printf ("type of [%d]: %s = ", i, lua_typename(L, typ));
      switch (typ) {
         case LUA_TNUMBER:
            printf("%g\n", lua_tonumber(L,i));
            break;
         case LUA_TSTRING:
            printf("%s\n", lua_tostring(L,i));
            break;
         default:
            printf("(unknown)\n");
      }
   }
   printf("\n ***** LUA STACK END *********\n");
}


static void report_lua_error(lua_State *L, int exc){/*{{{*/
      const char *comment = "";
      switch (exc){
	 case AIG_LUA_ERR_NUMBER_OF_ARGUMENTS:
	    comment = "invalid number of arguments";
	    break;
	 case AIG_LUA_ERR_ARGUMENT_TYPE:
	    comment = "invalid type of arguments";
	    break;
	 case AIG_LUA_ERR_INVALID_POINTER:
	    comment = "invalid pointer";
	    break;
	 case AIG_LUA_ERR_INCOMPATIBLE_OBJECT:
	    comment = "incompatible object";
	    break;
	 case AIG_LUA_ERR_CURVE_INSERTION_ERROR:
	    comment = "curve insertion error";
	    break;
	 case AIG_LUA_ERR_ILLEGAL_OVERLAP:
	    comment = "features overlap";
	    break;
	 case AIG_LUA_ERR_INVALID_QUALIFIER:
	    comment = "invalid qualifier";
	    break;
	 case AIG_LUA_ERR_FATAL:
	    comment = "unspecified fatal error";
	    break;
	 default:
	    comment = "undescribed error - This is a bug, please report it.";
      }
      luaL_error(L,"%s",comment);
}/*}}}*/

// ------------- Sample generation ------------------

static int l_new_curve(lua_State *L){/*{{{*/
   try{
      if (lua_gettop(L) != 2) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_isstring(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // curve type ("bezier" or "straight")

      //////////////////////// STRAIGHT SEGMENT ////////////////////////
      if (!strcmp(lua_tolstring(L,1,NULL),"segment")) {  // straight line will be created
	 // the next parameter must be table of two tables of two numbers
	 // e.g.: aig_create_curve("segment",{{1,2},{3,4}})

	 if (!lua_istable(L,2)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // bad parameters
	 DIST_TYPE coords[2][2];
	 for (int j=1; j<=2; j++)
	    for (int i=1; i<=2; i++){
	       lua_pushnumber(L,j); // get the j-th item
	       lua_gettable(L,2); // the coordinate table is the second parameter
	       lua_pushnumber(L,i); // get the i-th item of the table which is
	       // the j-th item of the outer table
	       if (!lua_istable(L,-2)) throw AIG_LUA_ERR_ARGUMENT_TYPE;
	       lua_gettable(L,-2);
	       if (!lua_isnumber(L, -1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // the item must be a number or else throw
	       coords[j-1][i-1] = (DIST_TYPE) lua_tonumber(L, -1);
	       lua_pop(L,2);
	    }
	 CCurve *curve = new CStraightLine(CLine(
	       CVector(coords[0][0],coords[0][1]),
	       CVector(coords[1][0],coords[1][1])
	       ));

	 lua_pushlightuserdata(L, (void *) curve);
	 return 1;
      }

      //////////////////////// BEZIER CURVE ////////////////////////
      if (!strcmp(lua_tolstring(L,1,NULL),"bezier")) {  // bezier curve will be created
	 // the next parameter must be table of four tables of two numbers
	 // e.g.: aig_create_curve("bezier",{{1,2},{3,4},{5,6},{7,8}})

	 if (!lua_istable(L,2)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // bad parameters
	 DIST_TYPE coords[4][2];
	 for (int j=1; j<=4; j++)
	    for (int i=1; i<=2; i++){
	       lua_pushnumber(L,j); // get the j-th item
	       lua_gettable(L,2); // the coordinate table is the second parameter
	       lua_pushnumber(L,i); // get the i-th item of the table which is
	       // the j-th item of the outer table
	       if (!lua_istable(L,-2)) throw AIG_LUA_ERR_ARGUMENT_TYPE;
	       lua_gettable(L,-2); 
	       if (!lua_isnumber(L, -1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // the item must be a number or else throw
	       coords[j-1][i-1] = (DIST_TYPE) lua_tonumber(L, -1);
	       lua_pop(L,2);
	    }
	 CCurve *curve = new CBezier(
	       CVector(coords[0][0],coords[0][1]),
	       CVector(coords[1][0],coords[1][1]),
	       CVector(coords[2][0],coords[2][1]),
	       CVector(coords[3][0],coords[3][1])
	       );

	 lua_pushlightuserdata(L, (void *) curve);
	 return 1;
      }
      throw AIG_LUA_ERR_INVALID_QUALIFIER;

   }
   
   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/

static int l_new_feature(lua_State *L){/*{{{*/
   // pars: array of pointers to curve, array of pointers to effect
   try{
      if (lua_gettop(L) != 3) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_istable(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // table of curves
      if (!lua_istable(L, 2)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // table of effects
      if (!lua_isnumber(L, 3)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // base graylevel

      IM_STORE_TYPE base_gl = lua_tonumber(L, 3);
      
      lua_pushnil(L);  /* first key */
      vector <CCurve *> curves;
      while (lua_next(L, 1) != 0) { // table is at index 1
	 if (!lua_islightuserdata(L,-1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // bad agrument type
	 CObject *ob = (CObject *) lua_topointer(L, -1);
	 if (!ob) throw AIG_LUA_ERR_INVALID_POINTER;

	 if (!ob->check_id(AIG_ID_CURVE)) throw AIG_LUA_ERR_INCOMPATIBLE_OBJECT; // non-curve object
	 curves.push_back((CCurve *)ob);
	 lua_pop(L, 1); // clean up
      }

      lua_pushnil(L);  /* first key */
      vector <CEffect *> effects;
      while (lua_next(L, 2) != 0) { // table is at index 2
	 if (!lua_islightuserdata(L,-1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // bad agrument type
	 CObject *ob = (CObject *) lua_topointer(L, -1);
	 if (!ob) throw AIG_LUA_ERR_INVALID_POINTER;

	 if (!ob->check_id(AIG_ID_EFFECT)) throw AIG_LUA_ERR_INCOMPATIBLE_OBJECT; // non-curve object
	 effects.push_back((CEffect *)ob);
	 lua_pop(L, 1); // clean up
      }

      CFeature *fe = new CFeature(curves, effects);
      fe->set_base_gray_level(base_gl);
      lua_pushlightuserdata(L, (void *) fe);
      return 1;
   }

   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/

static int l_paint_feature(lua_State *L){/*{{{*/
   // this function is mostly for debugging.
   // parameters: image pointer, feature pinter

   try{
      if (lua_gettop(L) != 2) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_islightuserdata(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // pointer to image
      if (!lua_islightuserdata(L, 2)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // pointer to feature

	 CObject *ob = (CObject *) lua_topointer(L, 1);
	 if (!ob) throw AIG_LUA_ERR_INVALID_POINTER;
	 if (!ob->check_id(AIG_ID_IMAGE)) throw AIG_LUA_ERR_INCOMPATIBLE_OBJECT; // bad object
	 CImage *im = (CImage *) ob;

	 ob = (CObject *) lua_topointer(L, 2);
	 if (!ob) throw AIG_LUA_ERR_INVALID_POINTER;
	 if (!ob->check_id(AIG_ID_FEATURE)) throw AIG_LUA_ERR_INCOMPATIBLE_OBJECT; // bad object
	 CFeature *fe = (CFeature *) ob;

	 fe->paint(im);
	 return 0;
   }


   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/

static int l_move_feature(lua_State *L){/*{{{*/
   // this function is mostly for debugging.
   // parameters: feature pointer, {dx, dy}

   try{
      if (lua_gettop(L) != 2) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_islightuserdata(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // pointer to feature
      if (!lua_istable(L, 2)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // shift vector (table)

      CObject *ob = (CObject *) lua_topointer(L, 1);
      if (!ob) throw AIG_LUA_ERR_INVALID_POINTER;
      if (!ob->check_id(AIG_ID_FEATURE)) throw AIG_LUA_ERR_INCOMPATIBLE_OBJECT; // bad object
      CFeature *fe = (CFeature *) ob;

      DIST_TYPE coords[2];
      for (int i=1; i<=2; i++){
	 lua_pushnumber(L,i); // get the i-th item of the table 
	 lua_gettable(L,2); // the coordinate table is the second parameter
	 if (!lua_isnumber(L, -1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // the item must be a number or else throw
	 coords[i-1] = (DIST_TYPE) lua_tonumber(L, -1);
	 lua_pop(L,1);
      }

      CVector shv(coords[0],coords[1]);
      fe->move_by(shv);
      return 0;
   }


   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/

static int l_new_sample(lua_State *L){/*{{{*/
   // args: sizex, sizey, array of pointers to feature
   try{
      if (lua_gettop(L) != 3) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_isnumber(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // parameter is sizex
      if (!lua_isnumber(L, 2)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // parameter is sizex
      if (!lua_istable(L, 3)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // parameter is the table of features


      DIST_TYPE sizex = lua_tonumber(L, 1);
      DIST_TYPE sizey = lua_tonumber(L, 2);

      lua_pushnil(L);  /* first key */
      vector <CFeature *> features;
      while (lua_next(L, 3) != 0) { // table is at index 3
	 if (!lua_islightuserdata(L,-1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // bad agrument type
	 CObject *ob = (CObject *) lua_topointer(L, -1);
	 if (!ob) throw AIG_LUA_ERR_INVALID_POINTER;

	 if (!ob->check_id(AIG_ID_FEATURE)) throw AIG_LUA_ERR_INCOMPATIBLE_OBJECT; // non-feature object
	 features.push_back((CFeature *)ob);
	 lua_pop(L, 1); // clean up
      }

      CSample *sa = NULL;
      try {
	 sa = new CSample(sizex, sizey, features);
      }
      catch (int sample_exception) {
	 if (sample_exception == AIG_EX_FEATURE_OVERLAP) throw AIG_LUA_ERR_ILLEGAL_OVERLAP;
      }
      lua_pushlightuserdata(L, (void *) sa);
      return 1;
   }

   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/
 
static int l_delete_sample(lua_State *L){/*{{{*/
   // arg: pointer to sample
   try{
      if (lua_gettop(L) != 1) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_islightuserdata(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // pointer to sample

      CSample *sa = (CSample *) lua_topointer(L, 1);
      if (!sa) throw AIG_LUA_ERR_INVALID_POINTER;
      if (((CObject *)sa)->check_id(AIG_ID_SAMPLE)) throw AIG_LUA_ERR_INCOMPATIBLE_OBJECT;

      delete sa;
      return 0;
   }

   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/

static int l_paint_sample(lua_State *L){/*{{{*/
   try{
      if (lua_gettop(L) != 2) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_islightuserdata(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // pointer to image
      if (!lua_islightuserdata(L, 2)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // pointer to sample

      CObject *ob = (CObject *) lua_topointer(L, 1);
      if (!ob) throw AIG_LUA_ERR_INVALID_POINTER;
      if (!ob->check_id(AIG_ID_IMAGE)) throw AIG_LUA_ERR_INCOMPATIBLE_OBJECT; // bad object
      CImage *im = (CImage *) ob;

      ob = (CObject *) lua_topointer(L, 2);
      if (!ob) throw AIG_LUA_ERR_INVALID_POINTER;
      if (!ob->check_id(AIG_ID_SAMPLE)) throw AIG_LUA_ERR_INCOMPATIBLE_OBJECT; // bad object
      CSample *sa = (CSample *) ob;

      sa->paint(im);
      return 0;
   }

   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/

static int l_new_effect(lua_State *L){/*{{{*/
   try{
      CEffect *ee = NULL;
      if (lua_gettop(L) < 2) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_isstring(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // type qualifier 
      // ( at this moment only "edge", "finestructure")

      if (!strcmp(lua_tolstring(L, 1, NULL), "edge")){
	 // generating edge effect
	 // args: "edge", ee_coe, ee_top_gl
	 
	 if (lua_gettop(L) != 3) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
	 if (!lua_isnumber(L, 2)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // ee coefficient
	 if (!lua_isnumber(L, 3)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // edge top gray level - base

	 float coef = (float) lua_tonumber(L, 2);
	 IM_STORE_TYPE edge_top_gl = (IM_STORE_TYPE) lua_tonumber(L, 3);

	 ee = new CEdgeEffect(coef, edge_top_gl);
	 if (!ee) throw AIG_LUA_ERR_FATAL;
	 lua_pushlightuserdata(L, (void *) ee);
	 return 1;
      }

      if (!strcmp(lua_tolstring(L, 1, NULL), "finestructure")){
	 // generating fine stricture effect
	 // args: "finestructure", density, min_r, max_r, min_coe, max_coe

	 if (lua_gettop(L) != 6) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
	 for (int i=2; i<=6; i++) if (!lua_isnumber(L, i)) throw AIG_LUA_ERR_ARGUMENT_TYPE;
	 
	 float density = lua_tonumber(L, 2);
	 DIST_TYPE min_r = lua_tonumber(L, 3);
	 DIST_TYPE max_r = lua_tonumber(L, 4);
	 float min_coe = lua_tonumber(L, 5);
	 float max_coe = lua_tonumber(L, 6);

	 ee = new CFineStructureEffect(density, min_r, max_r, min_coe, max_coe);

	 lua_pushlightuserdata(L, (void *) ee);
	 return 1;
      }

      throw AIG_LUA_ERR_INVALID_QUALIFIER;

   }

   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/

// ------------- Image processing ------------------

static int l_new_image(lua_State *L){/*{{{*/
   // pars: sizex, sizey
   try{
      if (lua_gettop(L) != 2) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_isnumber(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE;
      if (!lua_isnumber(L, 2)) throw AIG_LUA_ERR_ARGUMENT_TYPE;

      DIST_TYPE sizex = lua_tonumber(L, 1);
      DIST_TYPE sizey = lua_tonumber(L, 2);

      CImage *im = new CImage(sizex, sizey);

      lua_pushlightuserdata(L, (void *) im);
      return 1;
   }
   
   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/

static int l_delete_image(lua_State *L){/*{{{*/
   // pars: image pointer
   try{
      if (lua_gettop(L) != 1) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_islightuserdata(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // pointer to image

      CImage *im = (CImage *) lua_topointer(L, 1);
      if (!im) throw AIG_LUA_ERR_INVALID_POINTER;

      delete im;

      return 0;
   }
   
   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/

static int l_save_image(lua_State *L){/*{{{*/
   // pars: image pointer, filename, comment
   try{
      if (lua_gettop(L) != 3) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_islightuserdata(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // pointer to image
      if (!lua_isstring(L, 2)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // filename
      if (!lua_isstring(L, 3)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // comment

      CObject *ob = (CObject *) lua_topointer(L, 1);
      if (!ob) throw AIG_LUA_ERR_INVALID_POINTER;
      if (!ob->check_id(AIG_ID_IMAGE)) throw AIG_LUA_ERR_INCOMPATIBLE_OBJECT; // non-curve object
      CImage *im = (CImage *) ob;
      const char *fn = lua_tolstring(L, 2, NULL);
      const char *comment = lua_tolstring(L, 3, NULL);

      im->tiff_write(fn, comment, BI_8B);

      return 0;
   }
   
   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/

static int l_add_background_image(lua_State *L){/*{{{*/
   // pars: image pointer, {xdnesity,ydensity}, {x*y values}
   try{
      if (lua_gettop(L) != 3) throw AIG_LUA_ERR_NUMBER_OF_ARGUMENTS;
      if (!lua_islightuserdata(L, 1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // pointer to image
      if (!lua_istable(L, 2)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // filename
      if (!lua_istable(L, 3)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // comment

      CObject *ob = (CObject *) lua_topointer(L, 1);
      if (!ob) throw AIG_LUA_ERR_INVALID_POINTER;
      if (!ob->check_id(AIG_ID_IMAGE)) throw AIG_LUA_ERR_INCOMPATIBLE_OBJECT; // non-curve object
      CImage *im = (CImage *) ob;

      // getting densities {x,y}

      vector <double> densities;
      lua_pushnil(L);  /* first key */
      while (lua_next(L, 2) != 0) { // table is at index 2
	 if (!lua_isnumber(L,-1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // bad agrument type
	 densities.push_back(lua_tonumber(L,-1));
	 lua_pop(L, 1); // clean up
      }
      if (densities.size() != 2 ) throw AIG_LUA_ERR_ARGUMENT_TYPE;

      int densx = densities[0];
      int densy = densities[1];

      // getting values {v_1, ......., v_x*y}
      
      vector <IM_STORE_TYPE> values;
      lua_pushnil(L);  /* first key */
      while (lua_next(L, 3) != 0) { // table is at index 3
	 if (!lua_isnumber(L,-1)) throw AIG_LUA_ERR_ARGUMENT_TYPE; // bad agrument type
	 values.push_back(lua_tonumber(L,-1));
	 lua_pop(L, 1); // clean up
      }
      if (values.size() != densx*densy) throw AIG_LUA_ERR_ARGUMENT_TYPE;
      
      // construction the background
      
      CWavyBackgroud bg(densx, densy, values);
      bg.apply(im);
      return 0;
   }
   
   catch (t_aig_lua_err ex){
      report_lua_error(L, ex);
   }
   return 0;
}/*}}}*/


// lua function registration and lua code execution
int exec_lua_file(const char *fn){/*{{{*/
   lua_State *L = lua_open();
   luaL_openlibs(L);
   lua_register(L, "aig_new_curve", l_new_curve);
   lua_register(L, "aig_new_feature", l_new_feature);
   lua_register(L, "aig_paint_feature", l_paint_feature);
   lua_register(L, "aig_move_feature", l_move_feature);
   lua_register(L, "aig_new_sample", l_new_sample);
   lua_register(L, "aig_delete_sample", l_delete_sample);
   lua_register(L, "aig_paint_sample", l_paint_sample);
   lua_register(L, "aig_new_effect", l_new_effect);

   lua_register(L, "aig_new_image", l_new_image);
   lua_register(L, "aig_delete_image", l_delete_image);
   lua_register(L, "aig_save_image", l_save_image);
   lua_register(L, "aig_add_background_image", l_add_background_image);

   int err = luaL_dofile(L, fn);
   if (err) CLuaMessenger(AIG_MSG_FATAL_ERROR,lua_tostring(L,-1));
   lua_close(L);
   return err;
}/*}}}*/
#endif

// vim: cindent
