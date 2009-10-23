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

CLuaMessenger::CLuaMessenger(int msg_id, const char *comment){
   sender_id = "LUA";
   send_message(msg_id, comment);
}/*}}}*/

#ifdef HAVE_LUA
static int l_new_image(lua_State *L){/*{{{*/
   try{
      if (lua_gettop(L) != 2) throw -1;
      if (!lua_isnumber(L, 1)) throw -2;
      if (!lua_isnumber(L, 2)) throw -2;

      DIST_TYPE sizex = lua_tonumber(L, 1);
      DIST_TYPE sizey = lua_tonumber(L, 2);

      CImage *im = new CImage(sizex, sizey);

      lua_pushlightuserdata(L, (void *) im);
      return 1;
   }
   
   catch (int ex){

      const char *comment = "";
      switch (ex){
	 case -1:
	    comment = "aig_new_image - invalid number of arguments";
	    break;
	 case -2:
	    comment = "aig_new_image - invalid type of arguments";
	    break;
      }
      CLuaMessenger m((const char *)comment);
      lua_error(L);
   }
   return 0;
}/*}}}*/

static int l_delete_image(lua_State *L){/*{{{*/
   try{
      if (lua_gettop(L) != 1) throw -1;
      if (!lua_islightuserdata(L, 1)) throw -2; // pointer to image

      CImage *im = (CImage *) lua_topointer(L, 1);
      if (!im) throw -99;

      delete im;

      return 0;
   }
   
   catch (int ex){

      const char *comment = "";
      switch (ex){
	 case -1:
	    comment = "aig_save_image - invalid number of arguments";
	    break;
	 case -2:
	    comment = "aig_save_image - invalid type of arguments";
	    break;
	 case -99:
	    comment = "aig_save_image - invalid image pointer";
	    break;
      }
      CLuaMessenger m((const char *)comment);
      lua_error(L);
   }
   return 0;
}/*}}}*/

static int l_save_image(lua_State *L){/*{{{*/
   try{
      if (lua_gettop(L) != 3) throw -1;
      if (!lua_islightuserdata(L, 1)) throw -2; // pointer to image
      if (!lua_isstring(L, 2)) throw -2; // filename
      if (!lua_isstring(L, 3)) throw -2; // comment

      CImage *im = (CImage *) lua_topointer(L, 1);
      if (!im) throw -99;
      const char *fn = lua_tolstring(L, 2, NULL);
      const char *comment = lua_tolstring(L, 3, NULL);

      im->tiff_write(fn, comment, BI_8B);

      return 0;
   }
   
   catch (int ex){

      const char *comment = "";
      switch (ex){
	 case -1:
	    comment = "aig_save_image - invalid number of arguments";
	    break;
	 case -2:
	    comment = "aig_save_image - invalid type of arguments";
	    break;
	 case -99:
	    comment = "aig_save_image - invalid image pointer";
	    break;
      }
      CLuaMessenger m((const char *)comment);
      lua_error(L);
   }
   return 0;
}/*}}}*/

static int l_new_curve(lua_State *L){/*{{{*/
   try{
      if (lua_gettop(L) != 2) throw -1;
      if (!lua_isstring(L, 1)) throw -2; // curve type ("bezier" or "straight")

      //////////////////////// STRAIGHT SEGMENT ////////////////////////
      if (!strcmp(lua_tolstring(L,1,NULL),"segment")) {  // straight line will be created
	 // the next parameter must be table of two tables of two numbers
	 // e.g.: aig_create_curve("segment",{{1,2},{3,4}})

	 if (!lua_istable(L,2)) throw -2; // bad parameters
	 DIST_TYPE coords[2][2];
	 for (int j=1; j<=2; j++)
	    for (int i=1; i<=2; i++){
	       lua_pushnumber(L,j); // get the j-th item
	       lua_gettable(L,2); // the coordinate table is the second parameter
	       lua_pushnumber(L,i); // get the i-th item of the table which is
	       // the j-th item of the outer table
	       if (!lua_istable(L,-2)) throw -2;
	       lua_gettable(L,-2); // the coordinate table is the second parameter
	       if (!lua_isnumber(L, -1)) throw -2; // the item must be a number or else throw
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

	 CLuaMessenger("aig-new-curve Creating Bezier");
	 if (!lua_istable(L,2)) throw -2; // bad parameters
	 DIST_TYPE coords[4][2];
	 for (int j=1; j<=4; j++)
	    for (int i=1; i<=2; i++){
	       lua_pushnumber(L,j); // get the j-th item
	       lua_gettable(L,2); // the coordinate table is the second parameter
	       lua_pushnumber(L,i); // get the i-th item of the table which is
	       // the j-th item of the outer table
	       if (!lua_istable(L,-2)) throw -2;
	       lua_gettable(L,-2); // the coordinate table is the second parameter
	       if (!lua_isnumber(L, -1)) throw -2; // the item must be a number or else throw
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

   }
   
   catch (int ex){

      const char *comment = "";
      switch (ex){
	 case -1:
	    comment = "aig_new_curve - invalid number of arguments";
	    break;
	 case -2:
	    comment = "aig_new_curve - invalid type of arguments";
	    break;
      }
      CLuaMessenger m(comment);
      lua_error(L);
   }
   return 0;
}/*}}}*/

static int l_new_feature(lua_State *L){/*{{{*/
   try{
      if (lua_gettop(L) != 1) throw -1;
      if (!lua_istable(L, 1)) throw -2; // parameter is the table of curves

      lua_pushnil(L);  /* first key */
      vector <CCurve *> curves;
      while (lua_next(L, 1) != 0) { // table is at index 1
	 if (!lua_islightuserdata(L,-1)) throw -2; // bad agrument type
	 CObject *ob = (CObject *) lua_topointer(L, -1);

	 if (!ob->check_id(AIG_ID_CURVE)) throw -3; // non-curve object
	 curves.push_back((CCurve *)ob);
	 lua_pop(L, 1); // clean up
      }

      return 0;
   }

   catch (int ex){

      const char *comment = "";
      switch (ex){
	 case -1:
	    comment = "aig_new_feature - invalid number of arguments";
	    break;
	 case -2:
	    comment = "aig_new_feature - invalid type of arguments";
	    break;
	 case -3:
	    comment = "aig_new_feature - non-curve object in table";
	    break;
      }
      CLuaMessenger m(comment);
      lua_error(L);
   }
   return 0;
}/*}}}*/

int exec_lua_file(const char *fn){/*{{{*/
   lua_State *L = lua_open();
   luaL_openlibs(L);
   lua_register(L, "aig_new_image", l_new_image);
   lua_register(L, "aig_delete_image", l_delete_image);
   lua_register(L, "aig_save_image", l_save_image);
   lua_register(L, "aig_new_curve", l_new_curve);
   lua_register(L, "aig_new_feature", l_new_feature);
   int err = luaL_dofile(L, fn);
   lua_close(L);
   return err;
}/*}}}*/
#endif

// vim: cindent
