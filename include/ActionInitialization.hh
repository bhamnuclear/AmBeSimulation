#ifndef ACTION_INITIALIZATION_HH
#define ACTION_INITIALIZATION_HH

#include <G4VUserActionInitialization.hh>
#include <globals.hh>
#include "DetectorConstruction.hh"

class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization(int rank, int NumberOfThreads, DetectorConstruction* detectorConstruction, bool FissFragments, G4String RadioIsotope, bool InitialNeutrons);
  ~ActionInitialization();
  //virtual ~ActionInitialization();

  void Build() const override;

  void BuildForMaster() const override;

private:
  G4int fNumberOfThreads { -1 };
  G4bool fFissFragments=false;
  G4int fRank = 0;
  G4String fRadioIsotope;
  G4bool fInitialNeutrons;

  DetectorConstruction *fDetectorConstruction = nullptr;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
