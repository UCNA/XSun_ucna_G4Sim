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

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,+1.));
  fParticleGun->SetParticleEnergy(6.*MeV);
}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

//  G4double envSizeX = 0;
//  G4double envSizeY = 0;
//  G4double envSizeZ = 0;

//  G4double size = 0.8;
//  G4double x0 = size * envSizeX * (G4UniformRand()-0.5);
//  G4double y0 = size * envSizeY * (G4UniformRand()-0.5);
//  G4double z0 = -0.5 * envSizeZ;

  G4double x0 = 0;
  G4double y0 = 0;
  G4double z0 = -3*m;

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

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
