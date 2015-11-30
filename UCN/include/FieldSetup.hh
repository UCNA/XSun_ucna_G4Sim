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

  G4FieldManager* GetEastLocalFieldManager() { return fEastMWPCFieldManager;}
  G4FieldManager* GetWestLocalFieldManager() { return fWestMWPCFieldManager;}
protected:
  G4FieldManager* GetGlobalFieldManager();	// a method set to return the global field manager

  G4FieldManager* fFieldManager;
  G4FieldManager* fEastMWPCFieldManager;
  G4FieldManager* fWestMWPCFieldManager;	// Idk if I need the remaining objects to be East/West
						// So for safety's sake, I'll have duplicates of all.
  G4ChordFinder* fChordFinder;
  G4ChordFinder* fEastMWPCChordFinder;
  G4ChordFinder* fWestMWPCChordFinder;
  G4Mag_UsualEqRhs* fEquation;
  G4Mag_UsualEqRhs* fEastMWPCEquation;
  G4Mag_UsualEqRhs* fWestMWPCEquation;
  G4MagneticField* fMagneticField;
  G4MagneticField* fEastMWPCMagneticField;
  G4MagneticField* fWestMWPCMagneticField;
  G4MagIntegratorStepper* fStepper;
  G4MagIntegratorStepper* fEastMWPCStepper;
  G4MagIntegratorStepper* fWestMWPCStepper;

  G4double fMinStep;


};

#endif
