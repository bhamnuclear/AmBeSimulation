#include "PrimaryGeneratorAction.hh"
#include "Analysis.hh"
#include "NNDCLoader.hh"

#include <G4ParticleTable.hh>
#include <G4IonTable.hh>
#include <G4Event.hh>
#include <G4SystemOfUnits.hh>
#include <G4Geantino.hh>
#include <G4ParticleGun.hh>
#include <Randomize.hh>
#include <G4GeneralParticleSource.hh>
#include <G4Threading.hh>

#include <TMath.h>

//using namespace std;

PrimaryGeneratorAction::PrimaryGeneratorAction(RunAction* runAction, DetectorConstruction*  detectorConstruction, G4String IsotopeString)
: G4VUserPrimaryGeneratorAction(),
  fGPS(0),
  fEnvelopeBox(0),
  fRunAction(runAction),
  fDetectorConstruction(detectorConstruction),
  fRadioIsotope(IsotopeString)
{
  //messenger
  //f241AmGammaEnabled = true;
  f241AmGammaEnabled = false;
  PrimaryGeneratorMessenger();

  //GPS
  fGPS = new G4GeneralParticleSource();

  //G4ParticleDefinition* proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  //G4ParticleDefinition* electron = G4ParticleTable::GetParticleTable()->FindParticle("electron");
  G4ParticleDefinition* geantino = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
  G4ParticleDefinition* gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  G4ParticleDefinition* neutron = G4ParticleTable::GetParticleTable()->FindParticle("neutron");

  fGPS->GetCurrentSource()->SetNumberOfParticles(1);
  fGPS->GetCurrentSource()->SetParticleDefinition(neutron);
  fNeutronGPS = fGPS->GetCurrentSource();
  fNeutronGPS ->GetPosDist()->SetPosDisType("Point");
  fNeutronGPS ->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.));//mm
  //fNeutronGPS ->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0.,-1.,0.));//DEBUG
  fNeutronGPS ->GetAngDist()->SetAngDistType("iso");
  fNeutronGPS ->GetEneDist()->SetEnergyDisType("Mono");
  //fGPS->SetMultipleVertex(false);
  fGPS->SetMultipleVertex(true);

  fGPS->AddaSource(1);
  fGPS->GetCurrentSource()->SetNumberOfParticles(1);
  fCarbonGPS = fGPS->GetCurrentSource();
  fCarbonGPS->SetParticleDefinition(geantino);
  fCarbonGPS->GetPosDist()->SetPosDisType("Point");
  fCarbonGPS->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.));//mm
  //fCarbonGPS->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0.,0.,0.));//DEBUG
  fCarbonGPS->GetAngDist()->SetAngDistType("iso");
  fCarbonGPS->GetEneDist()->SetEnergyDisType("Mono");

  if (f241AmGammaEnabled)
  {
    fGPS->AddaSource(1);
    fGamma_59_54 = fGPS->GetCurrentSource();
    fGPS->GetCurrentSource()->SetNumberOfParticles(alphaPerNeutron*0.37);
    fGamma_59_54->SetParticleDefinition(gamma);
    fGamma_59_54->GetPosDist()->SetPosDisType("Point");
    fGamma_59_54->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.));//mm
    fGamma_59_54->GetAngDist()->SetAngDistType("iso");
    fGamma_59_54->GetEneDist()->SetEnergyDisType("Mono");
    fGamma_59_54->GetEneDist()->SetMonoEnergy(59.5409*keV);

    fGPS->AddaSource(1);
    fGamma_26_34 = fGPS->GetCurrentSource();
    fGPS->GetCurrentSource()->SetNumberOfParticles(alphaPerNeutron*0.0227);
    fGamma_26_34->SetParticleDefinition(gamma);
    fGamma_26_34->GetPosDist()->SetPosDisType("Point");
    fGamma_26_34->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.));//mm
    fGamma_26_34->GetAngDist()->SetAngDistType("iso");
    fGamma_26_34->GetEneDist()->SetEnergyDisType("Mono");
    fGamma_26_34->GetEneDist()->SetMonoEnergy(26.3446*keV);

    fGPS->AddaSource(1);
    fGamma_13_9 = fGPS->GetCurrentSource();
    fGPS->GetCurrentSource()->SetNumberOfParticles(alphaPerNeutron*0.359);
    fGamma_13_9->SetParticleDefinition(gamma);
    fGamma_13_9->GetPosDist()->SetPosDisType("Point");
    fGamma_13_9->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., 0.));//mm
    fGamma_13_9->GetAngDist()->SetAngDistType("iso");
    fGamma_13_9->GetEneDist()->SetEnergyDisType("Mono");
    fGamma_13_9->GetEneDist()->SetMonoEnergy(13.9*keV);
  }

  //Hard-code diff Xsection for 9Be(a,n)12C - Part of PrimaryGenerator
  //renaming these is useless as this simply gets the pointer to something else
  if (true)
  {
    gXS_0 = fRunAction->GetgXS_0();
    gXS_1 = fRunAction->GetgXS_1();
    gXS_2 = fRunAction->GetgXS_2();
    gstopping = fRunAction->Getgstopping();
    gXS_t = fRunAction->GetgXS_t();
    g_gs = fRunAction->Getg_gs();
    gp = fRunAction->Getgp();
    gpp = fRunAction->Getgpp();
  }

  if (fRadioIsotope=="241Am")
    fNNDCAlpha = new NNDCLoader("AmBeData/241AmAlpha.dat");
  else if (fRadioIsotope=="239Pu")
    fNNDCAlpha = new NNDCLoader("AmBeData/239PuAlpha.dat");

  r_max = fDetectorConstruction->GetAmBeSolid()->GetOuterRadius() - max_depth*mm;
  z_max = fDetectorConstruction->GetAmBeSolid()->GetZHalfLength() - max_depth*mm;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGPS, fNNDCAlpha;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //Randomise Centre //TODO
  G4double z = G4UniformRand()*(2*z_max)-z_max;
  G4double r = G4UniformRand()*r_max;
  G4double theta = G4UniformRand()*2*pi;
  G4double x, y;
  x = cos(theta)*r;
  y = sin(theta)*r;
  G4ThreeVector PositionVector = G4ThreeVector(x,y,z);

  fNeutronGPS ->GetPosDist()->SetCentreCoords(PositionVector);
  fCarbonGPS->GetPosDist()->SetCentreCoords(PositionVector);

  if (f241AmGammaEnabled)
  {
    G4ThreeVector GammaVector = PositionVector;//G4ThreeVector(cos(theta)*r_max,sin(theta)*r_max,z);
    fGamma_59_54->GetPosDist()->SetCentreCoords(GammaVector);
    fGamma_26_34->GetPosDist()->SetCentreCoords(GammaVector);
    fGamma_13_9->GetPosDist()->SetCentreCoords(GammaVector);
  }

  /* 241Am
   * 84.8% 5.48556 MeV
   * 13.1% 5.44280 MeV
   * Tot: 97.9%. Re-averaged: 87.62% and 13.38%
   */
  double QValue = 5.702049;//MeV - 9Be(a,n) Q value
  double Einitial = 5.48556;//Main peak
  //if(G4UniformRand()<0.1338) Einitial = 5.44280;//Satellite peak
  Einitial = fNNDCAlpha->GetEnergy(G4UniformRand());
  if (fDBG) G4cout << "PrimaryGenerator: Einitial="<<Einitial<<G4endl;//Debug
  /* //DBG - Check ratios
  if (true)
  {
    std::ofstream NNDCFile("NNDCOutput.dat", std::ios::out | std::ios::app);
    if (NNDCFile)
    {
      NNDCFile << Einitial<<"\n";
    }
    NNDCFile.close();
  }
  //*/
  if (Einitial==-1)
  {
    G4cout << "PrimaryGenerator: Invalid energy provided from NNDCLoader"<<G4endl;
    return;
  }
  double alphaEnergy = InteractionE(Einitial);//alpha of energy when it interacts - after Eloss

  G4double EnerNeutron = (GetNeutron(alphaEnergy))*MeV;
  if (EnerNeutron<0.) return;
  fNeutronGPS->GetEneDist()->SetMonoEnergy(EnerNeutron);
  //Save initital Neutron energy in the Ntuple
  fRunAction->FillInitialNeutron(EnerNeutron);
  //do not close row toa void data loss (otherwise try with closing every row, i.e. not conditioned on activation of detectors in EventAction)

  double _CarbonEnergy = alphaEnergy+QValue-EnerNeutron/MeV;//-ExcitationEnergy
  double CarbonEnergy = _CarbonEnergy*MeV;
  if (CarbonEnergy<0.*keV)
  {
    if (CarbonEnergy>-1.0*keV)//0>CarbonEnergy>-1.0keV - acceptable range, set to 0
    {CarbonEnergy=0.0*keV;}
    else
    {
      G4cout << "PrimaryGenerator: ";
      if (excited) G4cout << "Carbon in 1st excited state; ";
      if (excited2) G4cout << "Carbon in 2nd excited state; ";
      else G4cout << "Carbon in ground state; ";
      G4cerr << "Carbon Energy negative: "<< CarbonEnergy << "; neutron energy:"<< EnerNeutron <<G4endl;
    }
  }
  //Gamma from Carbon:
  //First test with just the gamma, then create a proper carbon (because of doppler broadening)
  //Also, for 2nd excited, check if hoyle state decays by alpha or just gamma. If just gamma, apply the branching ratio to gamma by creating this state using G4UniformRand
  if (excited or excited2)//1st carbon excited
  {
    //fGPS->SetMultipleVertex(false);
    if (excited)
    {
      if (CarbonEnergy-C_Excit_1ST<0) G4cerr << "PrimaryGenerator: Carbon energy below 1st excited"<<G4endl;
      //fCarbonGPS->GetEneDist()->SetMonoEnergy(4.43982*MeV);
      G4ParticleDefinition* C_Ion_1ST = G4IonTable::GetIonTable()->GetIon(Z,A,C_Excit_1ST);
      fCarbonGPS->SetParticleDefinition(C_Ion_1ST);
      fCarbonGPS->SetParticleCharge(ionCharge);
      fCarbonGPS->GetEneDist()->SetMonoEnergy(CarbonEnergy-C_Excit_1ST);
      if (fDBG) G4cout << "PrimaryGenerator: CarbonEnergy="<<CarbonEnergy/MeV<<" MeV, KE: "<<(CarbonEnergy-C_Excit_1ST)/MeV<<" MeV, 1st excited"<< G4endl;//Debug
    }
    else if (excited2)
    {
      //check if is within hoyle state
      if(G4UniformRand()<0.00416)//gamma from Hoyle
      {
        if (CarbonEnergy-C_Excit_2ND<0) G4cerr << "PrimaryGenerator: Carbon energy below 2nd excited"<<G4endl;
        //fCarbonGPS->GetEneDist()->SetMonoEnergy(7.65407*MeV);
        G4ParticleDefinition* C_Ion_2ND = G4IonTable::GetIonTable()->GetIon(Z,A,C_Excit_2ND);
        fCarbonGPS->SetParticleDefinition(C_Ion_2ND);
        fCarbonGPS->SetParticleCharge(ionCharge);
        fCarbonGPS->GetEneDist()->SetMonoEnergy(CarbonEnergy-C_Excit_2ND);
        if (fDBG) G4cout << "PrimaryGenerator: CarbonEnergy="<<CarbonEnergy/MeV<<" MeV, KE: "<<(CarbonEnergy-C_Excit_2ND)/MeV<<" MeV, 2nd excited"<< G4endl;//Debug
      }
      else//alpha decay
      {
        //set as null carbon as it should alpha decay
        G4ParticleDefinition* C_Ion_GND = G4IonTable::GetIonTable()->GetIon(Z,A,C_Excit_GND);
        fCarbonGPS->SetParticleDefinition(C_Ion_GND);
        fCarbonGPS->SetParticleCharge(ionCharge);
        //fCarbonGPS->GetEneDist()->SetMonoEnergy(CarbonEnergy-C_Excit_2ND);
        //G4cout << "PrimaryGenerator: CarbonEnergy="<<CarbonEnergy/MeV<<" MeV, KE: "<<(CarbonEnergy)/MeV<<" MeV, GND state"<< G4endl;//Debug
        fCarbonGPS->GetEneDist()->SetMonoEnergy(0);
        if (fDBG) G4cout << "PrimaryGenerator: Null Carbon"<< G4endl;//Debug
      }
    }
    else
    {
      G4cerr << "PrimaryGenerator: ERROR! - flags have an issue in PrimaryGeneratorAction::GeneratePrimaries()"<<G4endl;
      return;
    }
    //if (excited) G4cout << "Excited 1st , Primaries : " << anEvent->GetNumberOfPrimaryVertex() << G4endl;//DEBUG
    //else if (excited2) G4cout << "Excited 2nd , Primaries : " << anEvent->GetNumberOfPrimaryVertex() << G4endl;//DEBUG
  }
  else //GND state
  {
    //fCarbonGPS->GetEneDist()->SetMonoEnergy(0.0*MeV);
    G4ParticleDefinition* C_Ion_GND = G4IonTable::GetIonTable()->GetIon(Z,A,C_Excit_GND);
    fCarbonGPS->SetParticleDefinition(C_Ion_GND);
    fCarbonGPS->SetParticleCharge(ionCharge);
    fCarbonGPS->GetEneDist()->SetMonoEnergy(CarbonEnergy-C_Excit_GND);
    if (fDBG) G4cout << "PrimaryGenerator: CarbonEnergy="<<CarbonEnergy/MeV<<" MeV, KE: "<<(CarbonEnergy)/MeV<<" MeV, GND state"<< G4endl;//Debug
  }
  fGPS->GeneratePrimaryVertex(anEvent);
  //G4cout << "Primaries: " << anEvent->GetNumberOfPrimaryVertex() << G4endl;//DEBUG
  //if (anEvent->GetNumberOfPrimaryVertex()!=2) G4cout << "ERROR!" << G4endl; //-- All good //DBG
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::GetNeutron(double alphaEnergy)
{
  double Eint = alphaEnergy;//energy of interacting alpha
  //Sample theta CM
  thetaCM = InteractionTheta(Eint);
  if(thetaCM==0) return 0;//Skip this guy?
  else if (thetaCM == -1)
  {
    G4cout << "PrimaryGenerator::GetNeutron\tError when determining the reaction angle in the Centre of Mass frame. Sampling failed in InteractionTheta"<<G4endl;
    return -1;
  }
  thetalab = ConvertCMtoLab(Eint,thetaCM);
  double En = NeutronEnergy(Eint,thetalab,thetaCM);
  if(En==0) return 0;//Skip this guy
  return En;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PrimaryGeneratorAction::InteractionE(double Ep)
{
  //Rejection sampling for z_interaction based on XS(Ep) and stopping power
  //Then return Ep for z
  bool trials=true;
  //double XSmax = 1000;
  double XSmax = 700; //1.2*(TMath::MaxElement(gXS_0->GetN(),gXS_0->GetY())+TMath::MaxElement(gXS_1->GetN(),gXS_1->GetY())+TMath::MaxElement(gXS_2->GetN(),gXS_2->GetY()));
  /*
   * XSmax  x1    x1.2
   * 0+1+2  573   687
   * ...+3  638   766
   * ...+b  1072  1286
   */
  int loopcounter=0;
  while(trials) //Randomly choose a depth into the target (rather than energy to keep same target thickness)
  {
    double trial_z = max_depth*G4UniformRand();
    //Calc XS at Z
    double energy = Eloss(Ep,trial_z);
    if(energy<=Eth) continue;
    //double XS = gXS_t->Eval(energy);
    double XS = (gXS_0->Eval(energy)+gXS_1->Eval(energy)+gXS_2->Eval(energy));
    //cout<<gXS_1->Eval(energy)<<"\t"<<gXS_0->Eval(energy)<<endl;
    double choose_chan = G4UniformRand();
    double BRn0=gXS_0->Eval(energy)/XS;
    double BRn1=gXS_1->Eval(energy)/XS;
    double BRn2=gXS_2->Eval(energy)/XS;
    excited=false;
    excited2=false;
    if(choose_chan<BRn0)
    {
      excited=false;
      excited2=false;
    }
    else if(choose_chan<BRn0+BRn1)
    {
      excited=true;
      excited2=false;
    }
    else
    {
      excited=false;
      excited2=true;
    }
    if(XS>XSmax) {G4cout<<"PrimaryGenerator::IneractionE\tMax XS is exceeded!!"<<G4endl;}
    if(XS>XSmax*G4UniformRand())
    {
      return energy;
    }
    loopcounter++;
    if(loopcounter>1000)
    {
      G4cout<<"PrimaryGenerator::InteractionE\tRejection sampling not working properly - 1000 samples taken"<<G4endl;
      return 0;
    }
  }
  return 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PrimaryGeneratorAction::Eloss(double Ea, double z) //Energy of alpha of initial energy Ea after z mm in Be
{
  /*
   * Manually extracted from Lise++
   * Using 241Am at stechio 110
   * 16O at stechio 220
   * 9Be at stechio 890
   *
   * and alpha energy at max-> 5.48556
   *
   * -> max range 36.82486*um;
  */
  //if (fRadioIsotope=="241Am")
  //{
    //double range=0.000467*Ea*Ea+0.00238*Ea+0.000481;//mm
    double range=0.000528638*Ea*Ea+0.00355149*Ea+0.000974915;//mm
    if(z>range) return 0;//Fully stopped before Z
    range-=z;
    //double Eout=-3232.85*range*range+282.754*range+0.0932;//MeV
    double Eout=-1775.23*range*range+216.585*range-0.0837652;//MeV
  //}
  if (fRadioIsotope=="239Pu")
  {
    //double range=0.000467*Ea*Ea+0.00238*Ea+0.000481;//mm
    double range=0.000530474*Ea*Ea+0.00355288*Ea+0.00105571;//mm
    if(z>range) return 0;//Fully stopped before Z
    range-=z;
    //double Eout=-3232.85*range*range+282.754*range+0.0932;//MeV
    double Eout=-1609.33*range*range+211.354*range-0.0815178;//MeV
  }
  return Eout;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PrimaryGeneratorAction::InteractionTheta(double Ep) //COM theta
{
  //return (100+10*G4UniformRand())*d2r;
  double random_theta = acos(-1+2*G4UniformRand());
  //random_theta = 2.*pi*G4UniformRand();//NOTE:check this
  random_theta = pi*G4UniformRand();//TEST -> It is effetively as in two lines before acos(-1+2*G4UniformRand());
  double A[3]={1,0,0};
  bool sample=true;
  double XSmax=1.2;//Max is 1, XS data are in realtive yield
  int counter=0;
  TGraph2DErrors *select;
  if(!excited && !excited2) {select=g_gs;}
  if(excited && !excited2) {select=gp;}
  if(!excited && excited2) {select=gpp;}
  while(sample)
  {
    double XS_sample = XSmax*G4UniformRand();
    random_theta = acos(-1+2.*G4UniformRand());
    double XS=select->Interpolate(Ep,random_theta/d2r);
    if(XS>XSmax) G4cout<<"PrimaryGenerator::InteractionTheta\tSampling error for angle"<<G4endl;
    if(XS==0 || excited2) XS=1;//Outside of the energy range for this data set or Hoyle (poor data)
    if(XS_sample<XS)
    {
      return random_theta;//Isotropic for now
    }
    counter++;
    if(counter>100)
    {
      G4cout<<"PrimaryGenerator::InteractionTheta\tOver 100 samples taken - no solution for Ealpha = "<<Ep<<"\t"<<excited<<"\t"<<excited2<<G4endl;
      return 0;
    }
  }
  return -1;//should not reach //TODO:check if correct
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PrimaryGeneratorAction::ConvertCMtoLab(double Ep, double thetaCM) //Convert thetaCM to lab
{
  double Ex=C_Excit_GND_MeV;
  if(excited) Ex=C_Excit_1ST_MeV;
  if(excited2) Ex=C_Excit_2ND_MeV;

  double gamma = sqrt(((mb*ml)/(mt*mh))*(Ep/(Ep+(Q-Ex)*(1.+mb/mt))));//Calculate gamma for the conversion to lab COM
  gamma = gamma*((mb+mt)/(ml+mh));//Account for change in the CM velocity
  double thetalab = atan2(sin(thetaCM),(cos(thetaCM)+gamma));//
  secondks = false;
  if(gamma*cos(thetaCM)<=-1.) secondks=true;//Neutron is going bac
  //thetalab=(12*G4UniformRand())*d2r;//TEST
  //if(G4UniformRand()>0.5) secondks=true;//TEST
  return thetalab;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PrimaryGeneratorAction::NeutronEnergy(double Eb, double theta, double thetaCM) //Get neutron energy for a given Ep and theta
{
  double En;
  double Ex=C_Excit_GND_MeV;
  if(excited) Ex=C_Excit_1ST_MeV;//MeV
  if(excited2) Ex=C_Excit_2ND_MeV;

  double ECM = mt*Eb/(mt+mb)+Q-Ex;
  if(ECM<0)
  {
    G4cout<<"PrimaryGenerator::NeutronEnergy\tBelow threshold! Something is wrong!"<<G4endl;
    return 0;
  }
  double p_nT = sqrt(2.*ECM*ml*(mh/(ml+mh)));
  double phi=2*pi*G4UniformRand();
  double pn[3]={p_nT*sin(thetaCM)*cos(phi),p_nT*sin(thetaCM)*sin(phi),p_nT*cos(thetaCM)};//COM neutron momentum
  pn[2]+=(ml/(mh+ml))*sqrt(2.*mb*Eb);//Add the momentum boost for the lab frame
  double En_new = 0;
  for(int i=0;i<3;i++) En_new+=pn[i]*pn[i];
  En_new /= (2.*ml);
  En = En_new;
  theta = acos(pn[2]/p_nT);
  return En;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::PrimaryGeneratorMessenger()
{
  fMessenger = new G4GenericMessenger(this, "/AmBe/PrimaryGenerator/", "PrimaryGenerator control");
  fMessenger->DeclareProperty("241AmGamma",f241AmGammaEnabled);
  fMessenger->DeclareProperty("RadioIsotope",fRadioIsotope);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
