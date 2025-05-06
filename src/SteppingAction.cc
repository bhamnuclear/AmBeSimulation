#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"

#include <string.h>
#include <G4Step.hh>

#include <G4VPhysicalVolume.hh>
#include <G4VTouchable.hh>
#include <G4LogicalVolume.hh>

#include <G4Event.hh>
#include <G4RunManager.hh>

#include <G4UnitsTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4IonTable.hh>

SteppingAction::SteppingAction(EventAction* eventAction, bool FissFragments)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
  fFissFragments(FissFragments)
{
  fDBG = false;
  fDecayLimit = 30*year;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  if (aStep->GetTrack()->GetGlobalTime() > fDecayLimit)//kill step is time above fTimeLimit (30*years)
  {
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);//Disables scoring from any neutron
    return;
  };
  const G4StepPoint* postPoint = aStep->GetPostStepPoint();//EndPoint of step
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();//EndPoint of step
  const G4Track* aTrack = aStep->GetTrack();


  const DetectorConstruction* detectorConstruction = static_cast<const DetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  ///--- Volumes handling
  //get volume of the current step
  G4LogicalVolume* startVolume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  //if (startVolume->GetName() == "WorldLogical" and prePoint->GetStepStatus() == fGeomBoundary) return; //assume out of boundary
  if (!aTrack->GetNextVolume())
  {
    //G4cout << "particle out of World"<<G4endl;//DEBUG
    return; //next volume does not exist, i.e. OutOfWorld
  }

  G4LogicalVolume* endVolume = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4String startVolumeName = startVolume->GetName();
  G4String endVolumeName = endVolume->GetName();
  G4double eDepTotal = aStep->GetTotalEnergyDeposit();
  //G4double eDep = aStep->GetEnergyDeposit();
  G4double eDep = eDepTotal;


  G4String particleName = aTrack->GetParticleDefinition()->GetParticleName();

  if (particleName=="neutron" and startVolumeName=="ContainerLogical" and endVolumeName=="WorldLogical")
  {
    fEventAction->ScoreEmergingNeutron(aStep->GetPostStepPoint()->GetKineticEnergy());
  }

  //Tracking emerging secondaries if not Fission
  if (aTrack->GetParentID()>0)
  {
    if (startVolumeName=="ContainerLogical" and endVolumeName=="WorldLogical")
    {
      if(aTrack->GetParticleDefinition()->GetParticleName()=="neutron" and aTrack->GetCreatorProcess()->GetProcessName()!="nFissionHP")
      {
        fEventAction->AddSecondaryEmerging(aTrack->GetParticleDefinition(), aTrack->GetKineticEnergy());
      }
      if (fFissFragments)
      {
        if(aTrack->GetParticleDefinition()->GetParticleName()=="neutron" and aTrack->GetCreatorProcess()->GetProcessName()=="nFissionHP")
        {
          fEventAction->AddFissNeutEmerging(aTrack->GetKineticEnergy()/MeV);
        }
      }
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
