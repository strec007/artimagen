/* 
    Artificial SEM Image Generator (ArtImaGen)
    2009  Petr Cizmar @ National Institute of Standards and Technology
    E-mail: petr.cizmar@nist.gov

    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to Gnu.org public domain is compatible with GPL.

 */
#include "artimagen_i.h"
#include <cassert>
#include <ctime>
#include <cstdlib>
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
 
using namespace std;
using namespace artimagen;

CApp *AIGApp = NULL; // Global variable.

CApp::CApp(){/*{{{*/
   assert(!AIGApp); // Attempt to create another CApp class. Only once CApp instance is allowed per application.
   AIGApp = this;
   call_back = NULL; // the call_back function is not set (later used as indicator of being set)

   srandom(time(0)); // initialize the random generator.
}/*}}}*/

void CApp::set_message_call_back(void (*f)(t_message *)){/*{{{*/
   call_back = f;
}/*}}}*/

void CApp::broadcast_message(t_message message){/*{{{*/
   if (call_back) call_back(&message); // executes the call_back function
   // then the message is automatically destroyed.
}/*}}}*/

CObject::CObject(){/*{{{*/
   sender_id = "ILLEGAL"; // each class must rewrite this.
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

//
// vim: cindent