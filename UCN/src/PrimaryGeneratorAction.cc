#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC)
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fMyDetector(myDC),
  fPosOffset(),
  fSourceRadius(0),
  fRelToSourceHolder(false)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="gamma"/*"e-"*/);

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleTime(0.0*ns);	// Michael's has this line. Idk why.
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,+1.));
  fParticleGun->SetParticleEnergy(6.*MeV);
}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
// first line should be to set the random seed engine for Geant.

  G4double x0 = 0;	// In Michael's code he doesn't set the particle position
  G4double y0 = 0;	// unless he is throwing events from nuclear decay.
  G4double z0 = -3*m;	// Idk what it should be. Ask Brad.
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  if(anEvent -> GetEventID() == 0)
  {
    displayGunStatus();
  }
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::diskRandom(G4double radius, G4double& x, G4double& y)
{
  while(true)
  {
    x = (2.0*G4UniformRand()-1.)*radius;
    y = (2.0*G4UniformRand()-1.)*radius;
    if(x*x+y*y<=radius*radius) break;
  }
}

void PrimaryGeneratorAction::displayGunStatus()
{
  G4cout
  << fParticleGun->GetParticleDefinition()->GetParticleName()
  << " gun from " << fParticleGun->GetParticlePosition()/m
  << "m towards " << fParticleGun->GetParticleMomentumDirection()
  << " at " << fParticleGun->GetParticleTime()/ns
  << "ns : " << fParticleGun->GetParticleEnergy()/keV
  << "keV" <<
  G4endl;
}
