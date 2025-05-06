#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>

#include <G4RunManager.hh>
#include <G4Region.hh>
#include <G4RegionStore.hh>
#include <G4ProcessManager.hh>
#include <G4AutoDelete.hh>
#include <G4ParallelWorldPhysics.hh>

#include <G4PhysListFactory.hh>
#include <G4VPhysicsConstructor.hh>

#include <G4EmStandardPhysics_option4.hh>
#include <G4EmStandardPhysics.hh>
#include <G4EmExtraPhysics.hh>
#include <G4EmLivermorePhysics.hh>

#include <G4StoppingPhysics.hh>
#include <G4DecayPhysics.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4NeutronTrackingCut.hh>
#include <G4LossTableManager.hh>

#include <G4IonPhysics.hh>
#include <G4IonElasticPhysics.hh>

#include <G4HadronPhysicsQGSP_BIC_HP.hh>
#include <G4HadronPhysicsQGSP_BIC.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4HadronElasticPhysicsHP.hh>
#include <G4HadronPhysicsQGSP_BIC_AllHP.hh>
#include <G4HadronElasticPhysicsVI.hh>
#include <G4HadronPhysicsQGSP_BERT_HP.hh>
#include <G4HadronicInteraction.hh>
#include <G4HadronElasticPhysicsHPT.hh>

#include <G4Neutron.hh>



PhysicsList::PhysicsList(G4int ver):
    G4VModularPhysicsList()
{
  pMessenger = new PhysicsListMessenger(this);


  if ( ver > 0 ) G4cout << "<<< Geant4 Physics List simulation engine: QGSP_BERT_HPT with BinaryReaction" << G4endl << G4endl;

  defaultCutValue = 0.7*CLHEP::mm;
  SetVerboseLevel(ver);


  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics_option4(ver) );
  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver) );
  // Decays
  RegisterPhysics( new G4DecayPhysics(ver) );
  RegisterPhysics( new G4RadioactiveDecayPhysics(ver) );
  // Hadron Elastic scattering
  RegisterPhysics( new G4HadronElasticPhysicsHPT(ver) ); //Better for thermal neutrons w.r.t G4HadronElasticPhysicsVI
  // Hadron Physics
  RegisterPhysics( new G4HadronPhysicsQGSP_BERT_HP(ver) );//BERT
  //RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP(ver) );
  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver) );
  // Ion Physics
  RegisterPhysics( new G4IonElasticPhysics(ver) );
  RegisterPhysics( new G4IonPhysics(ver) );
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete pMessenger;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
  if (verboseLevel >1) G4cout << "QGSP_BERT_HP::SetCuts:";
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets  the default cut value for all particle types

  SetCutsWithDefault();

  //Set proton cut value to 0 for producing low energy recoil nucleus
  SetCutValue(0, "proton");
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void PhysicsList::ConstructProcess()
{
  // Transportation
  //
  AddTransportation();


  // Electromagnetic physics
  emPhysicsList->ConstructProcess();
  em_config.AddModels();
  decay_List->ConstructProcess();
  radioactiveDecay_List->ConstructProcess();

  // Hadronic physics
  for(size_t i=0; i < hadronPhys.size(); i++)
  {
    hadronPhys[i]->ConstructProcess();
  }
  // step limitation (as a full process)
  //AddStepMax();
  //Parallel world sensitivity
  //G4ParallelWorldPhysics* pWorld = new G4ParallelWorldPhysics("DetectorROGeometry");
  //pWorld->ConstructProcess();

  return;

  /*
  // example of GetHadronicModel (due to bug in QGSP_BIC_AllHP)
  //
  G4ProcessManager* pManager = G4Neutron::Neutron()->GetProcessManager();
  G4HadronicProcess* process = dynamic_cast<G4HadronicProcess*>(pManager->GetProcess("nCapture"));
  G4HadronicInteraction* model = process->GetHadronicModel("nRadCapture");
  if(model) model->SetMinEnergy(19.9*MeV);
  */
//}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>1)
  {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }
  G4cout << "NOT YET IMPLEMENTED"<<G4endl;
}
/*
  if (name == emName) return;

  ///////////////////////////////////
  //   ELECTROMAGNETIC MODELS
  ///////////////////////////////////
  if (name == "standard_opt4")
  {
    emName = name;
    delete emPhysicsList;
    hadronPhys.clear();
    emPhysicsList = new G4EmStandardPhysics_option4();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option4" << G4endl;
  }
  else if (name == "standard_opt3")
  {
    emName = name;
    delete emPhysicsList;
    hadronPhys.clear();
    emPhysicsList = new G4EmStandardPhysics_option3();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option3" << G4endl;
  }
  else if (name == "local_ion_ion_inelastic")
  {
    hadronPhys.push_back(new LocalIonIonInelasticPhysic());
    locIonIonInelasticIsRegistered = true;
  }
  else if (name == "standard")
  {
    emName = name;
    delete emPhysicsList;
    hadronPhys.clear();
    emPhysicsList = new G4EmStandardPhysics();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option4" << G4endl;
  }
  else if (name == "HADRONTHERAPY_1")
  {
    //The HADRONTHERAPY_1 physics list corresponds to the actual QGSP_BIC_HP list but with the following differences:
    AddPhysicsList("standard_opt4");
    hadronPhys.push_back( new G4DecayPhysics());
    hadronPhys.push_back( new G4RadioactiveDecayPhysics());
    hadronPhys.push_back( new G4IonBinaryCascadePhysics());
    hadronPhys.push_back( new G4EmExtraPhysics());
    hadronPhys.push_back( new G4HadronElasticPhysicsHP());
    hadronPhys.push_back( new G4StoppingPhysics());
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP());
    hadronPhys.push_back( new G4NeutronTrackingCut());


/*
 66   // EM Physics
 67   RegisterPhysics( new G4EmStandardPhysics_option4(ver) );
 68
 69   // Synchroton Radiation & GN Physics
 70   RegisterPhysics( new G4EmExtraPhysics(ver) );
 71
 72   // Decays
 73   RegisterPhysics( new G4DecayPhysics(ver) );
 74   RegisterPhysics( new G4RadioactiveDecayPhysics(ver) );
 75
 76   // Hadron Elastic scattering
 77   RegisterPhysics( new G4HadronElasticPhysicsHP(ver) );
 78
 79   // Hadron Physics
 80   RegisterPhysics(  new G4HadronPhysicsQGSP_BIC_HP(ver));
 81
 82   // Stopping Physics
 83   RegisterPhysics( new G4StoppingPhysics(ver) );
 84
 85   // Ion Physics
 86   RegisterPhysics( new G4IonElasticPhysics(ver) );
 87   RegisterPhysics( new G4IonPhysics(ver));

*/
/*
    G4cout << "HADRONTHERAPY_1 PHYSICS LIST has been activated" << G4endl;
  }
  else if (name == "HADRONTHERAPY_2")
  {
  */
    /*
         // The HADRONTHERAPY_2 physics list corresponds to the actual QGSP_BIC_HP list
235         // but with the following differences:
236         // --> G4EmStandardPhysics_option4 for the electromagnetic processes
237         //     is used in place of the less accurate G4EmStandardPhysics
238         // --> The G4RadioactiveDecayPhysics is added
239         // --> The 'local_ion_ion_inelastic' physics is used in place of the
240         //     G4IonBinaryCascadePhysics(): it used the QMD model to treat
241         //     the ion-ion inelastic interactions
242
243         AddPhysicsList("standard_opt4");
244         AddPhysicsList("local_ion_ion_inelastic");
245         hadronPhys.push_back( new G4DecayPhysics());
246         hadronPhys.push_back( new G4RadioactiveDecayPhysics());
247         hadronPhys.push_back( new G4EmExtraPhysics());
248         hadronPhys.push_back( new G4HadronElasticPhysicsHP());
249         hadronPhys.push_back( new G4StoppingPhysics());
250         hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP());
251
*/
    /*
    // HP models are switched off
    AddPhysicsList("standard_opt4");
    hadronPhys.push_back( new G4DecayPhysics());
    hadronPhys.push_back( new G4RadioactiveDecayPhysics());
    hadronPhys.push_back( new G4IonBinaryCascadePhysics());
    hadronPhys.push_back( new G4EmExtraPhysics());
    hadronPhys.push_back( new G4HadronElasticPhysics());
    hadronPhys.push_back( new G4StoppingPhysics());
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC());
    hadronPhys.push_back( new G4NeutronTrackingCut());

    G4cout << "HADRONTHERAPY_2 PHYSICS LIST has been activated" << G4endl;
  }
  else if (name == "QGSP_BERT_HP")
  {
    //Quite similar to Stock QGSP_BERT_HP
    defaultCutValue = 0.7*CLHEP::mm;
    cutForGamma = defaultCutValue;
    cutForElectron = defaultCutValue;
    cutForPositron = defaultCutValue;

    if(verboseLevel > 0)
    {
      G4cout << "<<< Geant4 Physics List simulation engine: QGSP_BERT_HP"<<G4endl;
      G4cout <<G4endl;
    }
    AddPhysicsList("standard");
    // Decays
    hadronPhys.push_back( new G4DecayPhysics(verboseLevel) );
    hadronPhys.push_back( new G4RadioactiveDecayPhysics(verboseLevel) );
    // Synchroton Radiation & GN Physics
    hadronPhys.push_back( new G4EmExtraPhysics(verboseLevel) );
    // Hadron Elastic scattering
    hadronPhys.push_back( new G4HadronElasticPhysicsVI(verboseLevel) );
    // Hadron Physics
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BERT_HP(verboseLevel));
    // Stopping Physics
    hadronPhys.push_back( new G4StoppingPhysics(verboseLevel));
    // Ion Physics
    hadronPhys.push_back( new G4IonPhysics(verboseLevel));
    hadronPhys.push_back( new G4IonElasticPhysics(verboseLevel) );

    //for (int i=0; i<hadronPhys.size(); i++)
    //  RegisterPhysics(hadronPhys[i]);
    G4cout << "QGSP_BERT_HP PHYSICS LIST has been activated" << G4endl;
  }
  else
  {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
    << " is not defined"
    << G4endl;
  }
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void PhysicsList::AddStepMax()
{
    // Step limitation seen as a process
    // This process must exist in all threads.
    //
    StepMax* stepMaxProcess  = new StepMax();
    
    auto particleIterator = GetParticleIterator();
    particleIterator->reset();
    while ((*particleIterator)()){
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        
        if (stepMaxProcess->IsApplicable(*particle) && pmanager)
        {
            pmanager ->AddDiscreteProcess(stepMaxProcess);
        }
    }
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
