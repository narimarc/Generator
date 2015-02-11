//____________________________________________________________________________
/* !

\program gvld_e_hadro_atten

\brief   Compares GENIE with hadron attenuation measurements from semi-inclusive
         electron-nucleus scattering.

         Syntax:
           gvld_e_hadro_atten 
                 [-g genie_input_file_list]  
                 [-d data_archive_location]

         Options:

           [] Denotes an optional argument.

           -d Full path to the electron scattering archive.
              By default, will pick the one distributed with GENIE.

           -g An XML file with GENIE inputs (cross sections and event samples).
              If not set, only data -no GENIE predictions- will be displayed.
              Multiple models can be included in the input file, each identified 
              by a "name" (all model predictions will be overlayed).
              For info on the XML file format see the GSimFiles class documentation.
              Notes:
              - The input event files are `gst' summary ntuples generated by 
                GENIE gntpc utility.
              - The inputs for this benchmark test are prepared by the script
                submit_eA_hadron_attenuation_validation_mc_jobs.pl
             
\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created June 06, 2008 

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <sstream>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TGraphAsymmErrors.h>
#include <TPostScript.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TText.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TChain.h>

#include "Conventions/GBuild.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/StringUtils.h"
#include "Utils/Style.h"
#include "Utils/SystemUtils.h"
#include "Utils/GSimFiles.h"

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::pdg;
using namespace genie::utils;

//____________________________________________________________________________
// Utility class to hold info on plotted datasets
class eHadroMultiplicityRatioDataSetDescription
{ 
public:
  eHadroMultiplicityRatioDataSetDescription(
     double E, int hpdg, int tgtpdg, 
     string citation, 
     double v_min  =  0,  double v_max  = 99999999,  /* define cuts if necessary */
     double z_min  =  0,  double z_max  =  1,        /* define cuts if necessary */
     double xF_min = -1,  double xF_max = +1,        /* define cuts if necessary */
     bool show = true) :
    fE         (E),
    fHadronPdg (hpdg),
    fTgtPdg    (tgtpdg),
    fvMin      (v_min),
    fvMax      (v_max),
    fzMin      (z_min),
    fzMax      (z_max),
    fxFMin     (xF_min),
    fxFMax     (xF_max),
    fCitation  (citation),
    fShow      (show)
  {
  }
  eHadroMultiplicityRatioDataSetDescription()
  {
  }
  int    TgtPdg   (void) const { return fTgtPdg; }
  int    TgtZ     (void) const { return pdg::IonPdgCodeToZ(fTgtPdg); }
  int    TgtA     (void) const { return pdg::IonPdgCodeToA(fTgtPdg); }
  int    HadronPdg(void) { return fHadronPdg; }
  string Citation (void) const { return fCitation; }
  double E        (void) const { return fE;     }
  double vMin     (void) const { return fvMin;  }
  double vMax     (void) const { return fvMax;  }
  double zMin     (void) const { return fzMin;  }
  double zMax     (void) const { return fzMax;  }
  double xFMin    (void) const { return fxFMin; }
  double xFMax    (void) const { return fxFMax; }
  string LabelTeX (void) const {  
    ostringstream label;
    label << "e + "<< this->TgtNameTeX() << " #rightarrow e + X, ";
    label << "E = " << fE << " GeV, R_{" << this->HadronNameTeX() << "}";
    label << " [" << fCitation << "]";
    return label.str();
  }
  string YAxisLabelTeX (void) const {  
    ostringstream label;
    label << "R_{" << this->HadronNameTeX() << "}";
    return label.str();
  }
  string TgtNameTeX  (void) const {
    if(fTgtPdg == 1000020040) return "^{4}He";
    if(fTgtPdg == 1000100200) return "^{20}Ne";
    if(fTgtPdg == 1000360840) return "^{84}Kr";
    if(fTgtPdg == 1000541320) return "^{132}Xe";
    return "Other";
  }                 
  string HadronNameTeX  (void) const {
    if(fHadronPdg == kPdgPiP   ) return "#pi^{+}";
    if(fHadronPdg == kPdgKP    ) return "K^{+}";
    if(fHadronPdg == kPdgProton) return "p";
    return "Other";
  }                 
private:
  double fE;          //   
  int    fHadronPdg;  //
  int    fTgtPdg;     //
  double fvMin;       //   
  double fvMax;       //   
  double fzMin;       //   
  double fzMax;       //   
  double fxFMin;      //   
  double fxFMax;      //   
  string fCitation;   // 
  bool   fShow;       //
};

/* 
..............................................................................
DATA
..............................................................................
*/
const int kNumOfDataSets = 12;

eHadroMultiplicityRatioDataSetDescription * kDataSet[kNumOfDataSets] = 
{
  // HERMES, 4He
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgPiP,    1000020040, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true),
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgKP,     1000020040, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true),
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgProton, 1000020040, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true),

  // HERMES, 20Ne
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgPiP,    1000100200, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true),
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgKP,     1000100200, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true),
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgProton, 1000100200, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true),

  // HERMES, 84Kr
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgPiP,    1000360840, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true),
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgKP,     1000360840, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true),
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgProton, 1000360840, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true),

  // HERMES, 132Xe
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgPiP,    1000541320, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true),
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgKP,     1000541320, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true),
  new eHadroMultiplicityRatioDataSetDescription(27.6, kPdgProton, 1000541320, "Airapetian:2007vu", 6.0, 999999., 0.2, 1.0, 0.0, 1.0, true)
};


// function prototypes
void               Init               (void);
void               End                (void);
TGraphAsymmErrors* Data               (int iset);
TH1D *             Model              (int iset, int imodel);
void               Draw               (int iset);
void               GetCommandLineArgs (int argc, char ** argv);
void               PrintSyntax        (void);

// command-line arguments
GSimFiles gOptGenieInputs;
string    gOptDataFilename = "";

// dbase information
const char * kDefDataFile = "data/validation/eA/hadronics/multiplicity_ratios/Rh.root";

// globals
TFile *       gRhDataFile = 0;
TTree *       gRhDataTree = 0;
TPostScript * gPS = 0;
TCanvas *     gC  = 0;
bool          gShowModel = false;

// plotting consts

const int kNCx = 2; // number of columns in TCanvas::Divide()
const int kNCy = 2; // number of rows    in TCanvas::Divide()

// model line styles
const int kNMaxNumModels = 5;
const int kLStyle    [kNMaxNumModels] = { 
   1,       2,        3,        5,            6 
}; 
string    kLStyleTxt [kNMaxNumModels] = { 
  "solid", "dashed", "dotted", "dot-dashed", "dot-dot-dashed" 
};

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc,argv);

  Init();

  for(int iset = 0; iset < kNumOfDataSets; iset++) 
  {
    Draw(iset);
  }

  End ();
  
  LOG("gvldtest", pINFO)  << "Done!";
  return 0;
}
//_________________________________________________________________________________
void Init(void)
{
  LOG("gvldtest", pNOTICE) << "Initializing...";;

  // Set GENIE style
  utils::style::SetDefaultStyle();

  // Get TTree with hadron multiplicity ratio data
  if( ! utils::system::FileExists(gOptDataFilename) ) {
      LOG("gvldtest", pFATAL)
         << "Can not find file: " << gOptDataFilename;
      gAbortingInErr = true;
      exit(1);
  }
  gRhDataFile = new TFile(gOptDataFilename.c_str(),"read");
  gRhDataTree = (TTree *) gRhDataFile->Get("hmr");
  if(!gRhDataTree) {
      LOG("gvldtest", pFATAL)
         << "Can not find TTree `hmr' in file: " << gOptDataFilename;
      gAbortingInErr = true;
      exit(1);
  }

  // Create plot canvas
  gC = new TCanvas("c","",20,20,500,650);
  gC->SetBorderMode(0);
  gC->SetFillColor(0);
  gC->SetGridx();
  gC->SetGridy();

  // Get local time to tag outputs
  string lt_for_filename   = utils::system::LocalTimeAsString("%02d.%02d.%02d_%02d.%02d.%02d");
  string lt_for_cover_page = utils::system::LocalTimeAsString("%02d/%02d/%02d %02d:%02d:%02d");

  // Create output ps file
  string filename  = Form("genie-hadron_multiplicity_ratio_data_comp-%s.ps",lt_for_filename.c_str());
  gPS = new TPostScript(filename.c_str(), 111);

  // Add cover page
  gPS->NewPage();
  gC->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText("GENIE comparison with hadron multiplicity ratio data from eA->eX");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  for(int imodel=0; imodel< gOptGenieInputs.NModels(); imodel++) {
    ostringstream stream;
    stream << "model tag: " << gOptGenieInputs.ModelTag(imodel)
           << " (" << kLStyleTxt[imodel] << " line)";
    hdr.AddText(stream.str().c_str());
  }
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(lt_for_cover_page.c_str());
  hdr.Draw();
  gC->Update();

  gC->SetLogx();
  gC->SetLogy();

}
//_________________________________________________________________________________
void End(void)
{
  LOG("gvldtest", pNOTICE) << "Cleaning up...";

  gPS->Close();

  delete gC;
  delete gPS;

  gRhDataFile->Close();
}
//_________________________________________________________________________________
TH1D* Model(int iset, int imodel)
{
  if(!gShowModel) return 0;

  LOG("gvldtest", pNOTICE) 
    << "Getting GENIE prediction (model ID = " 
    << imodel << ", data set ID = " << iset << ")";

  TChain * event_chain = gOptGenieInputs.EvtChain(imodel);
  if(!event_chain) {
     LOG("gvldtest", pNOTICE) 
        << "No corresponding event chain.";
     return 0;
  }

  const int    nbins = 20;
  const double Q2min =  0;
  const double Q2max = 10;

  TH1D * RhA = new TH1D("RhA","",nbins,Q2min,Q2max);
  TH1D * yD  = new TH1D("yD", "",nbins,Q2min,Q2max);

  double dE    = 1E-3;
  double E     = kDataSet[iset]->E();
  int    hpdg  = kDataSet[iset]->HadronPdg();
  int    A     = kDataSet[iset]->TgtA();
  int    Z     = kDataSet[iset]->TgtZ();
  double v_min = kDataSet[iset]->vMin();
  double z_min = kDataSet[iset]->zMin();

  ostringstream cut2He;
  ostringstream cutA;

  if      (hpdg == kPdgPiP) { cut2He << "nfpip*"; }
  else if (hpdg == kPdgPiP) { cut2He << "nfkp*";  }
  else if (hpdg == kPdgPiP) { cut2He << "nfp*";   }
  else { exit(1); }

  if      (hpdg == kPdgPiP) { cutA   << "nfpip*"; }
  else if (hpdg == kPdgPiP) { cutA   << "nfkp*";  }
  else if (hpdg == kPdgPiP) { cutA   << "nfp*";   }
  else { exit(1); }

  cut2He << "((abs(Ev-" << E << ")<" << dE 
         << ")&&(Z==1" 
         << ")&&(A==2" 
         << ")&&(Ev-El>" << v_min 
         << ")&&(Ef/(y*Ev)>" << z_min 
         << ")&&(pdgf==" << hpdg << "))";

  cutA   << "((abs(Ev-" << E << ")<" << dE 
         << ")&&(Z==" << Z
         << ")&&(A==" << A
         << ")&&(Ev-El>" << v_min 
         << ")&&(Ef/(y*Ev)>" << z_min 
         << ")&&(pdgf==" << hpdg << "))";

  event_chain->Draw("Q2>>yD", cut2He.str().c_str(),"goff");
  LOG("gvldtest", pNOTICE) 
     << "Used " << event_chain->GetSelectedRows() << " entries";
  event_chain->Draw("Q2>>RhA", cutA.str().c_str(),"goff");
  LOG("gvldtest", pNOTICE) 
        << "Used " << event_chain->GetSelectedRows() << " entries";
  RhA->Divide(yD);

  return RhA;
}
//_________________________________________________________________________________
TGraphAsymmErrors * Data(int iset)
{
  LOG("gvldtest", pNOTICE) 
    << "Loading experimental data set ID = " << iset;

  int hpdg  = kDataSet[iset]->HadronPdg();
  int A     = kDataSet[iset]->TgtA();
  int Z     = kDataSet[iset]->TgtZ();

  ostringstream cut;
  cut << "(Z==" << Z << ")&&(A==" << A << ")&&(hadron==" << hpdg << ")";

  gRhDataTree->Draw("Q2:R:dRp:dRm",cut.str().c_str(),"GOFF");
  LOG("gvldtest", pNOTICE)
    << "Applying selection: " << cut.str();

  int n = gRhDataTree->GetSelectedRows();
  LOG("gvldtest", pNOTICE)
    << "Found " << n << " data points in the data archive";

  if(n == 0) return 0; // return null graph

  // Data returned by TTree::Draw() are not necessarily ordered in Q2
  // Do the ordering here before building the graph
  int    *  idx  = new int   [n];
  double *  xv   = new double[n];
  double *  yv   = new double[n];
  double *  dypv = new double[n];
  double *  dymv = new double[n];

  TMath::Sort(n,gRhDataTree->GetV1(),idx,false);

  for(int i=0; i<n; i++) {
    int ii = idx[i];
    xv  [i] = gRhDataTree->GetV1()[ii];
    yv  [i] = gRhDataTree->GetV2()[ii];
    dypv[i] = gRhDataTree->GetV3()[ii];
    dymv[i] = gRhDataTree->GetV4()[ii];
  }

  TGraphAsymmErrors * gr = new TGraphAsymmErrors(n,xv,yv,0,0,dypv,dymv);

  utils::style::Format(gr, kBlack, kSolid, 1, kBlack, 20, 0.75);

  delete [] xv;
  delete [] yv;
  delete [] dypv;
  delete [] dymv;

  return gr;
}
//_________________________________________________________________________________
void Draw(int iset)
{
  // get data for the current comparison
  TGraph * data = Data(iset);

  // get the corresponding GENIE model prediction
  vector<TH1D*> models;
  for(int imodel=0; imodel< gOptGenieInputs.NModels(); imodel++) {
    TH1D * h = Model(iset,imodel);
    models.push_back(h);
  }

  if(models.size()==0 && !data) return;

  int plots_per_page = kNCx * kNCy;
  int iplot = 1 + iset % plots_per_page;
    
  if(iplot == 1) {   
     gPS->NewPage();
     gC -> Clear();
     gC -> Divide(kNCx,kNCy);
  }
  
  gC -> GetPad(iplot) -> Range(0,0,100,100);
  gC -> GetPad(iplot) -> SetFillColor(0);
  gC -> GetPad(iplot) -> SetBorderMode(0);
  gC -> GetPad(iplot) -> cd();

  double xmin = 0.0, scale_xmin = 0.5;
  double xmax = 0.0, scale_xmax = 1.2;
  double ymin = 0.0, scale_ymin = 0.4;
  double ymax = 0.0, scale_ymax = 1.2;
       
  TH1F * hframe = 0;
  bool have_frame = false;

  // have data points to plot?
  if(data) {
    // create frame from the data point range
    xmin  = ( data->GetX() )[TMath::LocMin(data->GetN(),data->GetX())];
    xmax  = ( data->GetX() )[TMath::LocMax(data->GetN(),data->GetX())];
    ymin  = ( data->GetY() )[TMath::LocMin(data->GetN(),data->GetY())];
    ymax  = ( data->GetY() )[TMath::LocMax(data->GetN(),data->GetY())];
    hframe = (TH1F*) gC->GetPad(iplot)->DrawFrame(
        scale_xmin*xmin, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
    hframe->Draw();
    have_frame = true;

    // draw current data set
    data->Draw("P");
  }//dbtable?

  // have model prediction to plot?
  if(models.size()>0) {
     for(int imodel=0; imodel<gOptGenieInputs.NModels(); imodel++) {

       TH1D * plot = models[imodel];
       if(plot) {
         int lsty = kLStyle[imodel];     
         utils::style::Format(plot,1,lsty,2,1,1,1);
         plot->Draw("CSAME");
       }
     }
  }//model?

  // axes label
  hframe->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  hframe->GetYaxis()->SetTitle(kDataSet[iset]->YAxisLabelTeX().c_str());

  // title
  TLatex * title = new TLatex(
     scale_xmin*xmin + 0.2*(scale_xmax*xmax-scale_xmin*xmin),
    1.01*scale_ymax*ymax,kDataSet[iset]->LabelTeX().c_str());
  title->SetTextSize(0.022);
  title->Draw();
 
  // update    
  gC->GetPad(iplot)->Update();
  gC->Update();
}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get data archive
  if(parser.OptionExists('d')){
     string filename = parser.ArgAsString('d');
     gOptDataFilename = filename;
  } else {
     if(gSystem->Getenv("GENIE")) {
        string base_dir = string( gSystem->Getenv("GENIE") );
        string filename = base_dir + kDefDataFile;
        gOptDataFilename = filename;
     } else {
        LOG("gvldtest", pFATAL)
          << "\n Please make sure that $GENIE is defined, or use the -d option"
          << "\n You didn't specify a data file and I can not pick the default one either";
        gAbortingInErr = true;
        exit(1);
     }
  }

  // get GENIE inputs
  gShowModel = true;
  if(parser.OptionExists('g')) {
     string inputs = parser.ArgAsString('g');
     bool ok = gOptGenieInputs.LoadFromFile(inputs);
     if(!ok) {
        LOG("gvldtest", pFATAL) 
           << "Could not read your validation program inputs from: " << inputs;
        gAbortingInErr = true;
        exit(1);
     }
  } else {
     gShowModel = false;
     LOG("gvldtest", pNOTICE) << " *** You didn't specify any GENIE MC outputs!";
     LOG("gvldtest", pNOTICE) << " *** Will plot only digitized expt data";
  }
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "  gvld_e_hadro_atten [-g genie_inputs.xml] [-d data_archive]\n";
}
//_________________________________________________________________________________

