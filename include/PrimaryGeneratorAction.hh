#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4ParticleGun.hh>
#include <G4GeneralParticleSource.hh>
#include <globals.hh>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <G4SystemOfUnits.hh>
#include <G4GenericMessenger.hh>

#include "RunAction.hh"
#include "NNDCLoader.hh"
#include "DetectorConstruction.hh"

class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class G4Box;
class G4GenericMessenger;
class RunAction;
class NNDCLoader;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(RunAction* runAction, DetectorConstruction* detectorConstruction, G4String IsotopeString);
  ~PrimaryGeneratorAction();
  virtual void GeneratePrimaries(G4Event* anEvent) override;
  double GetNeutron(double alphaEnergy);

  const G4ParticleGun *GetParticleGun() const { return fParticleGun; }
  //const G4GeneralParticleSource* GetParticleGun() const { return fGPS; }

  double InteractionE(double);
  double InteractionTheta(double);
  void makesumgraph();
  double Eloss(double,double);
  double ConvertCMtoLab(double Eint,double thetaCM);
  double NeutronEnergy(double Eint, double thetalab, double thetaCM);

  void PrimaryGeneratorMessenger();

private:
  G4GenericMessenger* fMessenger = nullptr;
  G4ParticleGun *fParticleGun;
  G4GeneralParticleSource *fGPS;
  G4Box *fEnvelopeBox;
  G4SingleParticleSource *fCarbonGPS;
  G4SingleParticleSource *fNeutronGPS;
  G4SingleParticleSource *fGamma_59_54, *fGamma_26_34, *fGamma_13_9 = nullptr;
  G4bool f241AmGammaEnabled;
  G4String fRadioIsotope = "241Am";

  DetectorConstruction *fDetectorConstruction = nullptr;
  RunAction *fRunAction = nullptr;
  NNDCLoader *fNNDCAlpha = nullptr;

  G4int fDBG {0};

  //NtupleIDs
  G4int fNtupleDetLayer = 0;

  TGraph *gXS_0;
  TGraph *gXS_1;
  TGraph *gXS_2;
  TGraph *gstopping;
  TGraph *gXS_t;
  TGraph2DErrors *g_gs;
  TGraph2DErrors *gp;
  TGraph2DErrors *gpp;

  double d2r=atan(1.)/45.;
  double pi=4.*atan(1.);
  double mb=4.002603;//Mass 4He
  double mt=9.012182;//Mass 9Be
  double ml=1.0008665;//Mass neutron
  double mh=12.;//Mass 12C
  double amu=931.4941;//MeV/u
  double Q=+5.702;//MeV
  double max_depth = 0.03639234;//mm according to LISE++ at E=Einitial for 241Am:110=15.999O:220=9.012Be:890 and density 2.9265704*g/cm3
  double Einitial=5.48556;//MeV
  double Eth=0. ;//Threshold energy
  double thetaCM,thetalab;
  bool secondks=false;
  bool excited=false;
  bool excited2=false;

  G4double r_max, z_max;

  G4double alphaPerNeutron = 35890/2.27;

  //12C data
  G4int Z = 6;
  G4int A = 12;
  G4double ionCharge = 0.*eplus;
  G4double C_Excit_GND = 0.*MeV;
  G4double C_Excit_1ST = 4.43982*MeV;
  G4double C_Excit_2ND = 7.65407*MeV;
  G4double C_Excit_GND_MeV = 0.;
  G4double C_Excit_1ST_MeV = 4.43982;
  G4double C_Excit_2ND_MeV = 7.65407;

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
