#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "ElectronBindingEnergy.hh"
#include "NuclEvtGen.hh"
#include "DetectorConstruction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"	// original example used these 3
#include "G4ParticleGun.hh"
#include "globals.hh"

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4ParticleGun.hh>
#include <G4Event.hh>
#include <G4VUserEventInformation.hh>

class G4ParticleGun;
class G4Event;

/// The primary generator action class with particle gun.

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(DetectorConstruction*);
    virtual ~PrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);

    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

  private:
    G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
    DetectorConstruction* fMyDetector;	// pointer to the detector geometry class

    double fSourceRadius;

    void DiskRandom(G4double radius, G4double& x, G4double& y);
    void DisplayGunStatus();
    void Set_113SnSource();

};

#endif


