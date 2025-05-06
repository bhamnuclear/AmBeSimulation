#include "RunAction.hh"
#include "Analysis.hh"

#include <G4Gamma.hh>
#include <G4Electron.hh>
#include <G4Neutron.hh>
#include <G4AccumulableManager.hh>

#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>

#include <G4Threading.hh>
#include <TROOT.h>
#include <TString.h>
#include <TGraph2DErrors.h>


#include "G4AutoLock.hh"

namespace
{
  G4Mutex aMutex = G4MUTEX_INITIALIZER;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(int rank, int NumberOfThreads, bool FissFragments, bool InitialNeutrons, G4String RadioIsotope):
  G4UserRunAction(),
  fRank(rank),
  fNumberOfThreads(NumberOfThreads),
  fFissFragments(FissFragments),
  fInitialNeutrons(InitialNeutrons),
  fRadioIsotope(RadioIsotope)
{
  //Messenger
  RunActionMessenger();


  fThreadid = G4Threading::G4GetThreadId();//omp_get_thread_num();
  if(fThreadid<0) fThreadid=0;

  //Multithreaded protected Code
  if (true)
  {
    ROOT::EnableThreadSafety();//see https://root-forum.cern.ch/t/mutexes-when-running-inside-geant4-threads/42060
    //mutex lock for multithreaded unsafe code
    G4AutoLock l(&aMutex);
    //Hard-code diff Xsection for 9Be(a,n)12C - Part of PrimaryGenerator
    gXS_0 = new TGraph("AmBeData/XS_0.txt");
    gXS_0->SetName("gXS_0"+fThreadid);
    gXS_0->SetTitle("gXS_0"+fThreadid);
    gXS_1 = new TGraph("AmBeData/XS_1.txt");
    gXS_1->SetName("gXS_1"+fThreadid);
    gXS_1->SetTitle("gXS_1"+fThreadid);
    gXS_2 = new TGraph("AmBeData/XS_2.txt");
    gXS_2->SetName("gXS_2"+fThreadid);
    gXS_2->SetTitle("gXS_2"+fThreadid);
    gstopping = new TGraph("AmBeData/stopping_power.txt");
    gstopping->SetName("gstopping"+fThreadid);
    gstopping->SetTitle("gstopping"+fThreadid);
    gXS_t = new TGraph("AmBeData/XS_t.txt"); //Atomic Data and Nuclear Data Tables, Vol. 15, No. 1, January 1975
    gXS_t->SetName("gXS_t"+fThreadid);
    gXS_t->SetTitle("gXS_t"+fThreadid);
    //G4cout << "First TGraph2D thread" <<fThreadid <<G4endl;//DBG

    g_gs = loadAng("AmBeData/n0_ang.dat");
    g_gs->SetName("g_gs"+fThreadid);
    g_gs->SetTitle("g_gs"+fThreadid);
    //G4cout << "second TGraph2D thread" <<fThreadid <<G4endl;//DBG
    gp = loadAng("AmBeData/n1_ang.dat");
    gp->SetName("gp"+fThreadid);
    gp->SetName("gp"+fThreadid);
    //G4cout << "third TGraph2D thread " <<fThreadid <<G4endl;//DBG
    gpp = loadAng("AmBeData/n2_ang.dat");
    gpp->SetName("gpp"+fThreadid);
    gpp->SetTitle("gpp"+fThreadid);
  }


  //Analysis via Root TTrees
  FileOutNameRank = fRadioIsotope.substr(3,5) + FileOutName+"-N"+fRank;
  fout.resize(fNumberOfThreads);
  tree.resize(fNumberOfThreads);
  if(fout[fThreadid]==0)
  {
    G4cout<<"THREAD: "<<FileOutNameRank+"-T"+fThreadid+".root"<<G4endl;
    fout[fThreadid]=new TFile(FileOutNameRank+"-T"+fThreadid+".root","RECREATE");

    tree[fThreadid] = new TTree(TString("FullBath"), "FullBath");//TTree name


    tree[fThreadid]->Branch("EmergingNeutrons",&fEmergingNeutrons);
    if (false)
    {
      tree[fThreadid]->Branch("secondaryGamma",&fbranchSecondaryGamma);
      tree[fThreadid]->Branch("secondaryElec",&fbranchSecondaryElectron);
    }
    tree[fThreadid]->Branch("secondaryNeut",&fbranchSecondaryNeutron);
    tree[fThreadid]->Branch("secondaryNeutEmerging",&fbranchSecondaryNeutronEmerging);
    if (fInitialNeutrons)
    {
      tree[fThreadid]->Branch("initN",&fbranchEInitN);
    }

    if (fFissFragments)
    {
      tree[fThreadid]->Branch("FissIon",&fbranchFissIon);
      tree[fThreadid]->Branch("FissNeut",&fbranchFissNeut);
      tree[fThreadid]->Branch("FissNeutEmerging",&fbranchFissNeutEmerging);
    }
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  G4cout<<"Trying to close thread "<<fThreadid<<" for "<<fout[fThreadid]<<G4endl;
  if (fout[fThreadid]!=0) fout[fThreadid]->Close();

  delete gXS_0, gXS_1, gXS_2, gstopping, gXS_t, g_gs, gp, gpp;

  if (IsMaster())//merge into one file per rank
  {
    //Merge root files
    G4String hadd_cmd = "hadd -f " + FileOutNameRank + ".root " + FileOutNameRank + "-T*.root";
    G4cout << "Running hadd command" << hadd_cmd.c_str() << G4endl;
    int haddFiles = system((hadd_cmd).c_str());

    if (haddFiles==0)
    {
      int removeFiles = system(("rm -vf "+FileOutNameRank+"-T*.root").c_str());
      if (removeFiles==0) G4cout << "Files merged and removed successfully" << G4endl;
      else G4cout << "Files merged but NOT removed" << G4endl;
    }
    else G4cout << "Files NOT merged - check manually" << G4endl;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  // Reset all accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  ClearBranches();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  //retrieve the number of events produced in the run
  G4int nofEvents = run->GetNumberOfEvent();

  //do nothing, if no events were processed
  if (nofEvents == 0) return;

  //Merge accumulables
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  G4cout<<"Finished thread: "<<fThreadid<<G4endl;
  //G4cout<<"Run ID: "<<thr<<G4endl;

  ///*
  if (IsMaster())
  {
    G4cout
     << "\n--------------------End of Global Run-----------------------"
     << " \n The run was " << nofEvents << " events " << G4endl;
  }

  //DO NOT replace with TreeWriteMidRun
  //Close File
  if(fout[fThreadid]!=0 && fThreadid!=-1)
  {
    //G4cout<<"FILL "<<fThreadid<<"\t"<<tree[fThreadid]<<"\t"<<fout[fThreadid]<<G4endl;
    //tree[fThreadid]->Fill();//Don't fill here as it would add an extra event
    tree[fThreadid]->Print();
    G4cout<<"Write"<<G4endl;
    tree[fThreadid]->Write();
    G4cout<<"CLOSE"<<G4endl;
  }
  ClearBranches();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TreeFill(G4int eventNumber)//end of event
{
  fEventNo = eventNumber;//this increases also for empty events, which are not recorded
  if (fout[fThreadid]!=0)
  {
    tree[fThreadid]->Fill();
    //G4cout<<"FILL "<<fThreadid<<"\t"<<tree[fThreadid]<<"\t"<<fout[fThreadid]<<G4endl;
  }
  //reset variables
  ClearBranches();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ClearBranches()
{
  //non-vecotr data can be simply assigned and overwritten
  fbranchEInitN = -1;
  fEventNo = -1;
  fbranchValid = false;

  fbranchSecondaryNeutron.clear();
  fbranchSecondaryNeutronEmerging.clear();
  fbranchSecondaryElectron.clear();
  fbranchSecondaryGamma.clear();

  fbranchFissionProduct.clear();
  fbranchFissionProductMass.clear();
  fbranchFissIon.clear();
  fbranchFissNeut.clear();
  fbranchFissNeutEmerging.clear();

  fEmergingNeutrons.clear();

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::RecordSecondaries(std::vector<Double_t> neutron, std::vector<Double_t> electron, std::vector<Double_t> gamma)
{
  fbranchSecondaryNeutron = neutron;
  fbranchSecondaryElectron = electron;
  fbranchSecondaryGamma = gamma;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::RecordSecondariesEmerging(std::vector<Double_t> neutron, std::vector<Double_t> electron, std::vector<Double_t> gamma)
{
  fbranchSecondaryNeutronEmerging = neutron;
  fbranchSecondaryElectronEmerging = electron;
  fbranchSecondaryGammaEmerging = gamma;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::loadFile(TString filename, std::vector<double>& Z_array, std::vector<double>& X_array, std::vector<double>& Y_array, std::vector<double>& ex_array, std::vector<double>& ey_array)
{
  std::ifstream file(filename);

  if (!file.is_open())
  {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  std::string line;
  while (std::getline(file, line))
  {
    std::stringstream ss(line);
    double Z, X, Y, ex, ey;

    if (ss >> Z >> X >> Y >> ex >> ey)
    {
      Z_array.push_back(Z);
      X_array.push_back(X);
      Y_array.push_back(Y);
      ex_array.push_back(ex);
      ey_array.push_back(ey);
    }
    else
    {
      std::cerr << "Warning: Could not parse line: " << line << std::endl;
    }
  }

  file.close();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::RunActionMessenger()
{
  fMessenger = new G4GenericMessenger(this, "/AmBe/run/", "Run control");

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TGraph2DErrors* RunAction::loadAng(TString filename)
{
    std::vector<double> Z_array, X_array, Y_array, ex_array, ey_array, ez_array;
    loadFile(filename, Z_array, X_array, Y_array, ex_array, ey_array);
	for (int i=0; i<Z_array.size(); i++) ez_array.push_back(0);

    // Display loaded data
	TGraph2DErrors *tg = new TGraph2DErrors(Z_array.size(),&Z_array[0],&X_array[0],&Y_array[0],&ez_array[0],&ex_array[0],&ey_array[0]);
	return tg;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
