#include "Run.hh"

Run::Run()
: G4Run(),
  fEdep(0.),
  fEdep2(0.)
{}


Run::~Run()
{}


void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  fEdep  += localRun->fEdep;
  fEdep2 += localRun->fEdep2;

  G4Run::Merge(run);
}


void Run::AddEdep (G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}


