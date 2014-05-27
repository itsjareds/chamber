//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
  //initialize cumulative quantities
  //
  sumEAbs = sum2EAbs =sumEGap = sum2EGap = 0.;
  sumLAbs = sum2LAbs =sumLGap = sum2LGap = 0.;
  sumEAbs2 = sum2EAbs2 = sumEGap2 = sum2EGap2 = 0.;
  sumLAbs2 = sum2LAbs2 = sumLGap2 = sum2LGap2 = 0.;
  sumEAbs3 = sum2EAbs3 = sumEGap3 = sum2EGap3 = 0.;
  sumLAbs3 = sum2LAbs3 = sumLGap3 = sum2LGap3 = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent(G4double EAbs, G4double EGap,
                                  G4double LAbs, G4double LGap,
			     G4double EAbs2, G4double EGap2,
				  G4double LAbs2, G4double LGap2,
			     G4double EAbs3, G4double EGap3,
				  G4double LAbs3, G4double LGap3)
{
  //accumulate statistic
  //
  sumEAbs += EAbs;  sum2EAbs += EAbs*EAbs;
  sumEGap += EGap;  sum2EGap += EGap*EGap;
  
  sumLAbs += LAbs;  sum2LAbs += LAbs*LAbs;
  sumLGap += LGap;  sum2LGap += LGap*LGap; 

  sumEAbs2 += EAbs2; sum2EAbs2 += EAbs2*EAbs2;
  sumEGap2 += EGap2; sum2EGap2 += EGap2*EGap2;

  sumLAbs2 += LAbs2; sum2LAbs2 += LAbs2*LAbs2;
  sumLGap2 += LGap2; sum2LGap2 += LGap2*LGap2;

  sumEAbs3 += EAbs3; sum2EAbs3 += EAbs3*EAbs3;
  sumEGap3 += EGap3; sum2EGap3 += EGap3*EGap3;

  sumLAbs3 += LAbs3; sum2LAbs3 += LAbs3*LAbs3;
  sumLGap3 += LGap3; sum2LGap3 += LGap3*LGap3; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  //compute statistics: mean and rms
  //
  sumEAbs /= NbOfEvents; sum2EAbs /= NbOfEvents;
  G4double rmsEAbs = sum2EAbs - sumEAbs*sumEAbs;
  if (rmsEAbs >0.) rmsEAbs = std::sqrt(rmsEAbs); else rmsEAbs = 0.;

  sumEAbs2 /= NbOfEvents; sum2EAbs2 /= NbOfEvents;
  G4double rmsEAbs2 = sum2EAbs2 - sumEAbs2*sumEAbs2;
  if (rmsEAbs2 >0.) rmsEAbs2 = std::sqrt(rmsEAbs2); else rmsEAbs2 =0.;

  sumEAbs3 /= NbOfEvents; sum2EAbs3 /= NbOfEvents;
  G4double rmsEAbs3 = sum2EAbs3 - sumEAbs3*sumEAbs3;
  if (rmsEAbs3 >0.) rmsEAbs3 = std::sqrt(rmsEAbs3); else rmsEAbs3 =0.;
  
  sumEGap /= NbOfEvents; sum2EGap /= NbOfEvents;
  G4double rmsEGap = sum2EGap - sumEGap*sumEGap;
  if (rmsEGap >0.) rmsEGap = std::sqrt(rmsEGap); else rmsEGap = 0.;

  sumEGap2 /= NbOfEvents; sum2EGap2 /= NbOfEvents;
  G4double rmsEGap2 = sum2EGap2 - sumEGap2*sumEGap2;
  if (rmsEGap2 >0.) rmsEGap2 = std::sqrt(rmsEGap2); else rmsEGap2 = 0.;
  
  sumEGap3 /= NbOfEvents; sum2EGap3 /= NbOfEvents;
  G4double rmsEGap3 = sum2EGap3 - sumEGap3*sumEGap3;
  if (rmsEGap3 >0.) rmsEGap3 = std::sqrt(rmsEGap3); else rmsEGap3 = 0.;

  sumLAbs /= NbOfEvents; sum2LAbs /= NbOfEvents;
  G4double rmsLAbs = sum2LAbs - sumLAbs*sumLAbs;
  if (rmsLAbs >0.) rmsLAbs = std::sqrt(rmsLAbs); else rmsLAbs = 0.;

  sumLAbs2 /= NbOfEvents; sum2LAbs2 /= NbOfEvents;
  G4double rmsLAbs2 = sum2LAbs2 - sumLAbs2*sumLAbs2;
  if (rmsLAbs2 >0.) rmsLAbs2 = std::sqrt(rmsLAbs2); else rmsLAbs2 = 0.;

  sumLAbs3 /= NbOfEvents; sum2LAbs3 /= NbOfEvents;
  G4double rmsLAbs3 = sum2LAbs3 - sumLAbs3*sumLAbs3;
  if (rmsLAbs3 >0.) rmsLAbs3 = std::sqrt(rmsLAbs3); else rmsLAbs3 = 0.;
  
  sumLGap /= NbOfEvents; sum2LGap /= NbOfEvents;
  G4double rmsLGap = sum2LGap - sumLGap*sumLGap;
  if (rmsLGap >0.) rmsLGap = std::sqrt(rmsLGap); else rmsLGap = 0.;

  sumLGap2 /= NbOfEvents; sum2LGap2 /= NbOfEvents;
  G4double rmsLGap2 = sum2LGap2 - sumLGap2*sumLGap2;
  if (rmsLGap2 >0.) rmsLGap2 = std::sqrt(rmsLGap2); else rmsLGap2 = 0.;

  sumLGap3 /= NbOfEvents; sum2LGap3 /= NbOfEvents;
  G4double rmsLGap3 = sum2LGap3 - sumLGap3*sumLGap3;
  if (rmsLGap3 >0.) rmsLGap3 = std::sqrt(rmsLGap3); else rmsLGap3 = 0.;
  
  //print
  //
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << "\n mean Energy in Absorber : " << G4BestUnit(sumEAbs,"Energy")
     << " +- "                          << G4BestUnit(rmsEAbs,"Energy")  
     << "\n mean Energy in Gap      : " << G4BestUnit(sumEGap,"Energy")
     << " +- "                          << G4BestUnit(rmsEGap,"Energy")
     << "\n mean Energy in Absorber2: " << G4BestUnit(sumEAbs2,"Energy")
     << " +- "				<< G4BestUnit(rmsEAbs2,"Energy")
     << "\n mean Energy in Gap2     : " << G4BestUnit(sumEGap2,"Energy")
     << " +- "                          << G4BestUnit(rmsEGap2,"Energy")
     << "\n mean Energy in Absorber3: " << G4BestUnit(sumEAbs3,"Energy")
     << " +- "				<< G4BestUnit(rmsEAbs3,"Energy")
     << "\n mean Energy in Gap3     : " << G4BestUnit(sumEGap3,"Energy")
     << " +- "                      	<< G4BestUnit(rmsEGap3,"Energy")
     << G4endl;
     
  G4cout
     << "\n mean trackLength in Absorber : " << G4BestUnit(sumLAbs,"Length")
     << " +- "                               << G4BestUnit(rmsLAbs,"Length")  
     << "\n mean trackLength in Gap      : " << G4BestUnit(sumLGap,"Length")
     << " +- "                               << G4BestUnit(rmsLGap,"Length")
     << "\n mean trackLength in Absorber2: " << G4BestUnit(sumLAbs2,"Length")
     << " +- "				     << G4BestUnit(rmsLAbs2,"Length")
     << "\n mean trackLength in Gap2     : " << G4BestUnit(sumLGap2,"Length")
     << " +- "				     << G4BestUnit(rmsLGap2,"Length")
     << "\n mean trackLength in Absorber3: " << G4BestUnit(sumLAbs3,"Length")
     << " +- "				     << G4BestUnit(rmsLAbs3,"Length")
     << "\n mean trackLength in Gap3     : " << G4BestUnit(sumLGap3,"Length")
     << " +- "				     << G4BestUnit(rmsLGap3,"Length")
     << "\n------------------------------------------------------------\n"
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
