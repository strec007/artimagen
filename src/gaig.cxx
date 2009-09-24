/*
 * =====================================================================================
 *
 *       Filename:  gaig.cxx
 *
 *    Description:  GUI wrapper for the libartimagen library.
 *
 *        Version:  0.1
 *        Created:  08/02/2009 06:57:25 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Dr. Petr Cizmar , petr@g.cizmar.org
 *     
 *
 * =====================================================================================
 */


#define IM_BORDER 100

#include <wx/wx.h>
#include <wx/notebook.h>
#include "libartimagen/artimagen.h"
#include "../config.h"
#include <string>
#include <vector>


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

enum{
   GAIG_TYPE_INT,
   GAIG_TYPE_STORE_TYPE,
   GAIG_TYPE_DIST_TYPE,
   GAIG_TYPE_FLOAT,
   GAIG_TYPE_DOUBLE
};

typedef struct{
   wxString ID;
   wxString descr;
   double initial_value;
   wxString format;
   wxTextCtrl *control;
} t_gaig_def;

void fill_def_struct(t_gaig_def *def, std::string ID, std::string descr, double initial_value, std::string format){
   def->ID = wxString::From8BitData(ID.c_str());
   def->descr = wxString::From8BitData(descr.data());
   def->initial_value = initial_value;
   def->format = wxString::From8BitData(format.data());
   def->control = NULL;
}

double get_def_value(std::vector <t_gaig_def> def, std::string ID){
   wxString IDwx = wxString::From8BitData(ID.c_str());
   for (unsigned int i=0; i<def.size(); i++){
      if (def[i].ID == IDwx){
	 return atof((const char *) def[i].control->GetValue().mb_str(wxConvUTF8));
      }
   }
   throw("Bad ID!");
}


using namespace artimagen;

// declarations of the artimagen variables
t_gc_definition gcdef;
IM_COORD_TYPE xsize;
IM_COORD_TYPE ysize;
IM_STORE_TYPE bg_min;
IM_STORE_TYPE bg_max;
IM_COORD_TYPE psf_sigma;
float astig_coef;
float astig_dir;
float dwell_time;
float pix_dead;
float line_dead;
float min_freq;
float max_freq;
float max_ampl;
int num_freq;
IM_STORE_TYPE noise_sigma;

enum{
   AIG_PREVIEW_BUT=9900,
   AIG_GENERATE_BUT
};

class GaigApp: public wxApp
{
   virtual bool OnInit();
};

class GaigFrame: public wxFrame/*{{{*/
{
   public:

      GaigFrame(const wxString& title, const wxPoint& pos, const wxSize& size);
      void Preview_Button_Pressed(wxCommandEvent &event);
      void Generate_Button_Pressed(wxCommandEvent &event);
      void generate(int preview);
      void setup_artimagen_parameters(int prev);


      void OnQuit(wxCommandEvent& event);
      void OnAbout(wxCommandEvent& event);

      std::vector <t_gaig_def> gaig_image_def;
      std::vector <t_gaig_def> gaig_gc_def;

      wxStaticBitmap *bm;

      DECLARE_EVENT_TABLE()
   private:
      void draw_form(std::vector <t_gaig_def> *pdef, wxScrolledWindow *pane, wxFlexGridSizer *sizer);
};/*}}}*/

enum
{
   ID_Quit = 1,
   ID_About,
};

   BEGIN_EVENT_TABLE(GaigFrame, wxFrame)
   EVT_MENU(ID_Quit, GaigFrame::OnQuit)
   EVT_MENU(ID_About, GaigFrame::OnAbout)
   EVT_BUTTON(AIG_PREVIEW_BUT, GaigFrame::Preview_Button_Pressed)
   EVT_BUTTON(AIG_GENERATE_BUT, GaigFrame::Generate_Button_Pressed)
   END_EVENT_TABLE()

IMPLEMENT_APP(GaigApp)

bool GaigApp::OnInit()/*{{{*/
{
   srandom(time(0));

   GaigFrame *frame = new GaigFrame( _T("gAIG - Artificial SEM Image Gerenrator"), wxPoint(50,50), wxSize(1024,700) );

   frame->Show(TRUE);
   SetTopWindow(frame);
   frame->generate(1);
   return TRUE;
} /*}}}*/

void GaigFrame::draw_form(vector <t_gaig_def> *pdef, wxScrolledWindow *pane, wxFlexGridSizer *sizer){/*{{{*/
   vector <t_gaig_def> def;
   def = *pdef;
   for (unsigned int i=0; i<def.size(); i++){
      if (def[i].ID != wxT("label")) {
	 sizer->Add(new wxStaticText(pane, -1, def[i].descr),0, wxALIGN_CENTRE_VERTICAL, 0);
	 def[i].control = new wxTextCtrl(pane, -1,  wxT(""));
	 sizer->Add(def[i].control, 0, wxLEFT | wxRIGHT, 4);
	 def[i].control->SetValue(wxString::Format(def[i].format, def[i].initial_value));
      } else {
	 wxStaticText *stxt = new wxStaticText(pane, -1, def[i].descr);
	 sizer->Add(stxt ,0, wxALIGN_CENTRE_VERTICAL, 0);
	 stxt->SetFont(wxFont(10, wxDECORATIVE, wxNORMAL, wxBOLD, 0, wxT("")));
	 sizer->Add(new wxStaticText(pane, -1, wxT("")),0, wxALIGN_CENTRE_VERTICAL, 0);

      }
   }
   *pdef = def;
}/*}}}*/

   GaigFrame::GaigFrame(const wxString& title, const wxPoint& pos, const wxSize& size)/*{{{*/
: wxFrame((wxFrame *)NULL, -1, title, pos, size)
{
   

   t_gaig_def def;

   fill_def_struct(&def,"label","Size:", 0, ""); gaig_image_def.push_back(def);
   fill_def_struct(&def,"sizex","Image size (one side) [pix]:", 1024, "%0.0f"); gaig_image_def.push_back(def);
   //fill_def_struct(&def,"sizey","Image Y-size [pix]:", 1024, "%0.0f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"label","Background:", 0, ""); gaig_image_def.push_back(def);
   fill_def_struct(&def,"bg_min","Background minimum [GL]:", 0.1, "%0.2f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"bg_max","Background maximum [GL]:", 0.3, "%0.2f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"label","Blur:", 0, ""); gaig_image_def.push_back(def);
   fill_def_struct(&def,"psf_sigma","Blur-PSF sigma [pix]:", 1, "%0.2f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"astig_coef","Astigmatism ratio [1]:", 1, "%0.2f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"astig_dir","Astigmatism direction [deg]:", 60, "%0.2f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"label","Drift:", 0, ""); gaig_image_def.push_back(def);
   fill_def_struct(&def,"dwell_time","Pixel dwell time [ns]:", 50, "%0.1f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"pix_dead","Pixel dead time [ns]:", 50, "%0.1f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"line_dead","Line dead time [ns]:", 250, "%0.1f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"min_freq","Minimun drift frequency [Hz]:", 3, "%0.1f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"max_freq","Maximun drift frequency [Hz]:", 180, "%0.1f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"max_ampl","Maximum drift amplitude [pix]:", 0.2, "%0.2f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"num_freq","Number of drift frequencies:", 10, "%0.0f"); gaig_image_def.push_back(def);
   fill_def_struct(&def,"label","Noise:", 0, ""); gaig_image_def.push_back(def);
   fill_def_struct(&def,"noise_sigma","Gaussian noise sigma [GL]:", 0.07, "%0.3f"); gaig_image_def.push_back(def);
   
   fill_def_struct(&def,"label","Edge Effect:", 0, ""); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"ee_coe","Edge-effect coefficient", 0.2, "%0.2f"); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"ee_gl","Edge-effect level [GL]:", 0.3, "%0.2f"); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"base_gl","Base gray-level [GL]:", 0.5, "%0.2f"); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"base_var","Base-level variation [GL]:", 0.1, "%0.2f"); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"label","Grain Sizes and Density:", 0, ""); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"gc_min","Gold-grain min size [pix]:", 5, "%0.1f"); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"gc_max","Gold-grain max size [pix]:", 15, "%0.1f"); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"gc_dens","Grain density:", 0.5, "%0.2f"); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"label","Fine Structure:", 0, ""); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"fs_dens","Fine-structure density x10^-4 [spots/pix]:", 7, "%0.2f"); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"fs_minrad","Fine-structure min radius [pix]:", 6, "%0.1f"); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"fs_maxrad","Fine-structure max radius [pix]:", 10, "%0.1f"); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"fs_minamp","Fine-structure min amplif. :", 0.95, "%0.2f"); gaig_gc_def.push_back(def);
   fill_def_struct(&def,"fs_maxamp","Fine-structure max amplif. :", 1.05, "%0.2f"); gaig_gc_def.push_back(def);

   // MENU 
   wxMenu *menuFile = new wxMenu;

   menuFile->Append( ID_About, _T("&About...") );
   menuFile->AppendSeparator();
   menuFile->Append( ID_Quit, _T("E&xit") );

   wxMenuBar *menuBar = new wxMenuBar;
   menuBar->Append( menuFile, _T("&File") );

   SetMenuBar( menuBar );

   // STATUS BAR
 //    CreateStatusBar();
 //    SetStatusText( _T("Welcome to wxWindows!") );
	
   // ALL CONTROLS
   wxBoxSizer *appsizer = new wxBoxSizer(wxHORIZONTAL);
   
        //LEFT SIDE 
   wxBoxSizer *leftsizer = new wxBoxSizer(wxVERTICAL);
   appsizer -> Add(leftsizer, 0, wxALL, 10);

   bm = new wxStaticBitmap(this, -1, 0, wxDefaultPosition, wxSize(512,512), 0, wxT(""));
   leftsizer -> Add(bm, 0, wxALL, 10);

   wxBoxSizer *butsizer = new wxBoxSizer(wxHORIZONTAL);
   leftsizer -> Add(butsizer, 0, wxALL, 10);
   
   butsizer -> Add(new wxButton(this, AIG_PREVIEW_BUT, wxT("&Preiew")), 0, wxALL, 5);
   butsizer -> Add(new wxButton(this, AIG_GENERATE_BUT, wxT("&Generate")), 0, wxALL, 5);


   	//RIGHT SIDE
   wxNotebook *def_notebook = new wxNotebook(this, wxID_ANY, wxDefaultPosition, wxSize(480,-1), wxNB_LEFT);
   appsizer -> Add(def_notebook, 0, wxALL, 10);

		// IMAGE PANE
   wxScrolledWindow *image_pane = new wxScrolledWindow(def_notebook, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL);
   def_notebook->AddPage(image_pane, wxT("Image defs."));
   wxFlexGridSizer *image_sizer = new wxFlexGridSizer(13, 2, 2, 2); // rows, cols, vgap, hgap 

   draw_form(&gaig_image_def, image_pane, image_sizer);
   image_pane->SetSizer(image_sizer);

		// SAMPLE PANE
   wxScrolledWindow *sample_pane = new wxScrolledWindow(def_notebook, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL);
   def_notebook->AddPage(sample_pane, wxT("GoC defs."));
   wxFlexGridSizer *sample_sizer = new wxFlexGridSizer(12, 2, 2, 2); // rows, cols, vgap, hgap 

   draw_form(&gaig_gc_def, sample_pane, sample_sizer);
   sample_pane->SetSizer(sample_sizer);

// FINISH dialog definition
   SetSizer(appsizer);
}/*}}}*/

void GaigFrame::OnQuit(wxCommandEvent& WXUNUSED(event))/*{{{*/
{
   Close(TRUE);
}/*}}}*/

void GaigFrame::OnAbout(wxCommandEvent& WXUNUSED(event))/*{{{*/
{
   wxMessageBox(_T("This is the gaig. Graphical interface to the\n libartimagen SEM artificial image generator.\n By Petr Cizmar <pcizmar@nist.gov>, 2009\n\nThis program is public domain."),
	 _T("About gaig"), wxOK | wxICON_INFORMATION, this);
}/*}}}*/

void GaigFrame::Preview_Button_Pressed(wxCommandEvent &event)/*{{{*/
{
   event.Skip();
   generate(1); // 1 - preview will be generated 
}/*}}}*/

void GaigFrame::Generate_Button_Pressed(wxCommandEvent &event)/*{{{*/
{
   event.Skip();
   generate(0); //0 - generate and save image
}/*}}}*/

void GaigFrame::generate(int preview){/*{{{*/
   wxBusyCursor wait;

   setup_artimagen_parameters(preview); 

   gcdef.sizex +=IM_BORDER;
   gcdef.sizey +=IM_BORDER;

   //create the image
   CImage *im = new CImage(gcdef.sizex, gcdef.sizey);
   //create background
   CWavyBackgroud *bg = new CWavyBackgroud(bg_min, bg_max, 5,5);
   //apply backgroud
   bg->apply(im);
   delete bg;
 
   //create the sample
   try {
      CGoldOnCarbonSample *sam = new CGoldOnCarbonSample(&gcdef);
      //paint it to the image
      sam->paint(im);
      delete sam;
   }
   catch (int ex){
      if (ex == AIG_EX_GOLDONCARBON_TOO_MANY_FAILS) wxMessageBox(wxT("Too many failed feature placements.\n Try to decrease the sizes of \nthe features, or decrease the density."));
      return;
   }

   // create the PSF
   if (psf_sigma >= 1) {
      CGaussianPsf *psf = new CGaussianPsf(gcdef.sizex, gcdef.sizey, psf_sigma, astig_coef, astig_dir); 
      psf->apply(im);
      delete psf;
   } else {
      if (psf_sigma > 0) wxMessageBox(wxT("The blur-sigma (beam size) is too small. Minimum is 1 pixel.\nContinuing with no blur."));
   }

   // create the drift-distortion function and apply it
   t_vib_definition vibdef;
   vibdef.pixel_dwell_time = dwell_time;
   vibdef.min_frequency = min_freq;
   vibdef.max_frequency = max_freq;
   vibdef.max_amplitude = max_ampl;
   vibdef.number_of_frequencies = num_freq;
   vibdef.pixel_dead_time = pix_dead;
   vibdef.line_dead_time = line_dead;

   CVibration *vib = new CVibration(&vibdef);
   vib->apply(im, 0);
   delete vib;

   gcdef.sizex -=IM_BORDER;
   gcdef.sizey -=IM_BORDER;

   //crop the image
   im->crop(IM_BORDER/2, IM_BORDER/2, gcdef.sizex+IM_BORDER/2, gcdef.sizey+IM_BORDER/2);

   // apply Gaussian noise
   

   CGaussianNoise *gn = new CGaussianNoise(noise_sigma);
   gn->apply(im);
   delete gn;

   if (preview){
      // convert image to 8-bit (originally stored in doubles)
      CSaveBuffer<TYPE_8B> *sb = new CSaveBuffer<TYPE_8B>(im);

      TYPE_8B *wx_image_data = (TYPE_8B *) malloc(3*gcdef.sizex*gcdef.sizey); // must use malloc, wxImage frees it itself
      TYPE_8B *buffer_data = (TYPE_8B *) sb->give_save_buffer();
      for (int i=0; i<gcdef.sizex*gcdef.sizey; i++) 
	 for (int j=0; j<3; j++) wx_image_data[3*i+j] = buffer_data[i]; 
      delete sb;

      wxImage *wxim = new wxImage(gcdef.sizex, gcdef.sizey, wx_image_data, false);

      wxBitmap *bitmap = new wxBitmap(*wxim);
      delete wxim;

      bm->SetBitmap(*bitmap);
      delete bitmap;
   } else {
      wxFileDialog dlg(this, wxT("Save Image..."), wxT(""), wxT(""), wxT("TIFF files (*.tiff;*.tif)|*.tiff;*.tif"), wxFD_SAVE);
      if (dlg.ShowModal() == wxID_OK){
	 im->tiff_write((const char *) dlg.GetPath().mb_str(wxConvUTF8),"Saved by gAIG; (c) Petr Cizmar", BI_8B);
      }
   }

   delete im;

}/*}}}*/

void GaigFrame::setup_artimagen_parameters(int prev){/*{{{*/

   double sizex, sizey;
   if (prev) { // preview image will be generated -> size = 512x512
      sizex = 512;
      sizey = 512;
   } else {
      sizex = get_def_value(gaig_image_def, "sizex");
      sizey = get_def_value(gaig_image_def, "sizex"); // FIX this after fixinf the size bug
 //     sizex = FLVAL(xsize_ctl);
 //     sizey = FLVAL(ysize_ctl);
   }

   gcdef.sizex = sizex;
   gcdef.sizey = sizey;
   gcdef.ee_coefficient = get_def_value(gaig_gc_def, "ee_coe");
   gcdef.ee_top_above_basic = get_def_value(gaig_gc_def, "ee_gl");
   gcdef.basic_level = get_def_value(gaig_gc_def, "base_gl");
   gcdef.basic_level_variation = get_def_value(gaig_gc_def, "base_var");
   gcdef.grain_min_size = get_def_value(gaig_gc_def, "gc_min");
   gcdef.grain_max_size = get_def_value(gaig_gc_def, "gc_max");
   gcdef.number_of_grains = sizex*sizey/(50*50)*get_def_value(gaig_gc_def, "gc_dens");
   gcdef.fs_density = 1e-4 * get_def_value(gaig_gc_def, "fs_dens");
   gcdef.fs_min_r = get_def_value(gaig_gc_def, "fs_minrad");
   gcdef.fs_max_r = get_def_value(gaig_gc_def, "fs_maxrad");
   gcdef.fs_min_coe = get_def_value(gaig_gc_def, "fs_minamp");
   gcdef.fs_max_coe = get_def_value(gaig_gc_def, "fs_maxamp");

   bg_min = get_def_value(gaig_image_def, "bg_min");
   bg_max =  get_def_value(gaig_image_def, "bg_max");
   psf_sigma =  get_def_value(gaig_image_def, "psf_sigma");
   astig_coef =  get_def_value(gaig_image_def, "astig_coef");
   astig_dir =  get_def_value(gaig_image_def, "astig_dir");
   dwell_time =  get_def_value(gaig_image_def, "dwell_time");
   pix_dead =  get_def_value(gaig_image_def, "pix_dead");
   line_dead =  get_def_value(gaig_image_def, "line_dead");
   min_freq =  get_def_value(gaig_image_def, "min_freq");
   max_freq =  get_def_value(gaig_image_def, "max_freq");
   max_ampl =  get_def_value(gaig_image_def, "max_ampl");
   num_freq =  get_def_value(gaig_image_def, "num_freq");
   noise_sigma =  get_def_value(gaig_image_def, "noise_sigma");

   }/*}}}*/
