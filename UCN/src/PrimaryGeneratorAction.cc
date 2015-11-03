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
}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
// first line should be to set the random seed engine for Geant.

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle;

  int r1;
  r1 = rand() % 1007246 + 1;	// This bound is # of digits I want to produce
				// NOTE	not set to 100 because on nndc we get 100.7% for 391 keV gammas.

  double percentage = r1/10000.;	// This gives us 0.001 precision.
  G4cout << "random number: " << r1 << G4endl;
  G4cout << "percentage: " << percentage << G4endl;

  if((percentage >= 0) && (percentage <= 64.97))
  {
    fParticleGun -> SetParticleEnergy(391.698*keV);
    particle = particleTable->FindParticle(particleName="gamma");
  }
  else if((percentage > 64.97) && (percentage <= 93.77))
  {
    fParticleGun -> SetParticleEnergy(363.758*keV);
    particle = particleTable -> FindParticle(particleName="e-");
  }
  else if((percentage > 93.77) && (percentage <= 99.37))
  {
    fParticleGun -> SetParticleEnergy(387.461*keV);
    particle = particleTable -> FindParticle(particleName="e-");
  }
  else if((percentage > 99.37) && (percentage <= 100.507))
  {
    fParticleGun -> SetParticleEnergy(390.872*keV);
    particle = particleTable -> FindParticle(particleName="e-");
  }
  else if((percentage > 100.507) && (percentage <= 100.712))
  {
    fParticleGun -> SetParticleEnergy(391.576*keV);
    particle = particleTable -> FindParticle(particleName="e-");
  }
  else if((percentage > 100.712) && (percentage >= 100.7246))
  {
    fParticleGun -> SetParticleEnergy(391.697*keV);
    particle = particleTable -> FindParticle(particleName="e-");
  }
  else
  {
    G4cout << "Random number sampled beyond the scope of the decay." << G4endl;
  }

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleTime(0.0*ns);        // Michael's has this line. Idk why.
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,+1.));


  G4double x0 = 0;	// In Michael's code he doesn't set the particle position
  G4double y0 = 0;	// unless he is throwing events from nuclear decay.
  G4double z0 = -1*m;	// Idk what it should be. Ask Brad.
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  displayGunStatus();
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
