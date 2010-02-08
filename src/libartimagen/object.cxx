/* 
    Artificial Charged-Particle-Microscope Image Generator (ARTIMAGEN)
    2009  Petr Cizmar @ National Institute of Standards and Technology
    E-mail: petr.cizmar@nist.gov

    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to GNU.org public domain is compatible with GPL.

 */
#define AIGAPP_DECLARATION 1  // this file declares the global AIGApp variable

#include "artimagen_i.h"
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <vector>
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
 
artimagen::CApp *AIGApp = NULL; // global variable.

using namespace std;
using namespace artimagen;

CApp::CApp(){/*{{{*/
   assert(!AIGApp); // Attempt to create another CApp class. Only once CApp instance is allowed per application.
   AIGApp = this;
   call_back = NULL; // the call_back function is not set (later used as indicator of being set)

   srandom(time(0)); // initialize the random generator.
#ifdef HAVE_FFTW3_THREADS
   nthreads = 1;
   fftw_init_threads(); // initialize the FFTW3 threads
#endif
}/*}}}*/

void CApp::set_message_call_back(void (*f)(t_message *)){/*{{{*/
   call_back = f;
}/*}}}*/

void CApp::broadcast_message(t_message message){/*{{{*/
   if (call_back) call_back(&message); // executes the call_back function
   // then the message is automatically destroyed.
}/*}}}*/

void CApp::set_number_of_threads(int nthreads){/*{{{*/
#ifdef HAVE_FFTW3_THREADS
   this->nthreads = nthreads;
#endif
}/*}}}*/

int CApp::get_number_of_threads(){/*{{{*/
#ifdef HAVE_FFTW3_THREADS
   return this->nthreads;
#else
   return -1; // no FFTW3 threads support
#endif
}/*}}}*/

CObject::CObject(){/*{{{*/
//   sender_id = "ILLEGAL"; // each class must rewrite this.
//   object_id = AIG_ID_OBJECT;
   messaging_enabled = 1;
}/*}}}*/

void CObject::send_message(int message, const char *comment){/*{{{*/
   if (messaging_enabled){
      t_message mes;
      mes.sender_id = sender_id; // build the message
      mes.message = message;
      mes.comment = comment;
      if (AIGApp) AIGApp->broadcast_message(mes); // broadcast it
   }
}/*}}}*/

void CObject::enable_messages(){/*{{{*/
   messaging_enabled = 1;
}/*}}}*/

void CObject::disable_messages(){/*{{{*/
   messaging_enabled = 0;
}/*}}}*/

int CObject::check_id(int id){/*{{{*/
   if (id == object_id) return AIG_CHECK_ID_EXACT;
   for (int i=0; i < inherited_ids.size(); i++){
      if(inherited_ids[i] == id) return AIG_CHECK_ID_INHERITED;
   }
   return AIG_CHECK_ID_NO_MATCH;
}/*}}}*/

void CObject::ident(int id){/*{{{*/
   object_id = id;
   inherited_ids.push_back(id);

}/*}}}*/

//
// vim: cindent
