/* 
    Artificial SEM Image Generator (ArtImaGen)
    2009  Petr Cizmar @ National Institute of Standards and Technology
    E-mail: petr.cizmar@nist.gov

    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to Gnu.org public domain is compatible with GPL.

 */


#include <config.h>
#include <iostream>
#include <cstdlib>
#include "libartimagen/artimagen_i.h"

using namespace std;
using namespace artimagen;

void message_call_back(t_message *msg){ // this is the call-back function that will print messages
   cout << msg->sender_id << ": " << msg->message << " and comments: " << msg->comment << endl;
}


int main(int argc, char **argv){
   cout << "\n **** ArtImaGen v." << PACKAGE_VERSION << " ****\n * by Petr Cizmar @ NIST\n * 2008-2009\n\n";
   if (argc != 2){
      cerr << "Usage: " << argv[0] << " defs.aig" << endl << endl;
      exit(-1);
   }

   CApp app;
   app.set_message_call_back(message_call_back); // set the call-back function

   int err = exec_lua_file(argv[1]);
   if (err) {
      cerr << "Error #" << err << " in the .aig file." << endl << endl;
   };

return 0;
}

// vim: cindent
