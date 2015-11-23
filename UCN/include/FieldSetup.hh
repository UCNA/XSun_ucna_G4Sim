#ifndef FieldSetup_h
#define FieldSetup_h 1

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;

class FieldSetup
{
public:
  FieldSetup();
  ~FieldSetup();

  void UpdateField();
  void SetFieldValue(G4double fieldStrength);

  G4FieldManager* GetLocalFieldManager() { return fLocalFieldManager;}
protected:
  G4FieldManager* GetGlobalFieldManager();	// a method set to return the global field manager

  G4FieldManager* fFieldManager;
  G4FieldManager* fLocalFieldManager;
  G4ChordFinder* fChordFinder;
  G4ChordFinder* fLocalChordFinder;
  G4Mag_UsualEqRhs* fEquation;
  G4Mag_UsualEqRhs* fLocalEquation;
  G4MagneticField* fMagneticField;
  G4MagneticField* fLocalMagneticField;
  G4MagIntegratorStepper* fStepper;
  G4MagIntegratorStepper* fLocalStepper;

  G4double fMinStep;


};

#endif
