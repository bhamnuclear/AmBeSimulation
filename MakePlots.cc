/*
 * Root macro?
 *
 */

#include <iostream>
#include <vector>
//ROOT libraries
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TSystem.h>



//constants
Double_t RunTime = 100;//s

int MakePlots(TFile* _file0);
int MakePlots(TString inname="");


//code
int MakePlots(TFile* _file0)
{
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame tree_FullBath("FullBath", _file0);

  std::vector<TString> BranchList = {"EmergingNeutrons", "initN", "secondaryNeut","secondaryNeutEmerging","FissIon", "FissNeut", "FissNeutEmerging"};

  TString filename = _file0->GetName();
  TString date = filename(0,4)+filename(5,2)+filename(8,2);
  TString filename_base = filename(11,filename.Length()-5-11);
  filename_base.ReplaceAll("-","_");
  TString filename_new = filename_base+"_"+date+"_";

  TCanvas c("c","x hist");
  ROOT::RDF::RResultPtr<TH1D> myHisto;

  for (unsigned int i=0; i<BranchList.size(); i++)
  {
    myHisto = tree_FullBath.Histo1D({"histName", "histTitle", 120, 0., 24.}, BranchList[i]);
    myHisto->Draw();
    c.SaveAs(filename_new+BranchList[i]+".C");
  }

  return 0;
}

int MakePlots(TString inname)
{
 ///---OPEN AND READ FILE---///
  TString filename;
  if (inname=="")
  {
    std::cout << "Enter filename to analyse" << "\n";
    std::cin >> filename;
  }
  else
  {
    filename=inname;
  }
  std::cout <<"Now analysing: " << filename << "\n";//MSG
  //add root extension if missing
  if (filename(filename.Length()-5,filename.Length()) != ".root") filename=filename+".root";

  //initialise TTree read
  if(gSystem->AccessPathName(filename))
  {
    std::cout << "File does not exist" << std::endl;//MSG
    return 1;
  }
  TFile* file = new TFile(filename, "read");

  return MakePlots(file);
}

int main()
{
  MakePlots();
  return 0;
}
