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
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
     /*
     void SetAbsorberMaterial (G4String);     
     void SetAbsorberThickness(G4double);
     void SetAbsorber2Material (G4String);
     void SetAbsorber2Thickness(G4double);
     void SetAbsorber3Material (G4String);
     void SetAbsorber3Thickness(G4double);    

     void SetGapMaterial (G4String);     
     void SetGapThickness(G4double);
     void SetGap2Material (G4String);
     void SetGap2Thickness(G4double);
     void SetGap3Material (G4String);
     void SetGap3Thickness(G4double);

     void SetCalorSizeYZ(G4double);          
     void SetNbOfLayers (G4int);   
      
     void SetMagField(G4double);
     */
     
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
     void PrintCalorParameters(); 
     
     
     G4double GetWorldSizeX()           {return WorldSizeX;}; 
     G4double GetWorldSizeYZ()          {return WorldSizeYZ;};
     
     
     //G4double GetCalorThickness()       {return CalorThickness;}; 
     G4double GetCalorSizeYZ()          {return CalorSizeYZ;};
     
     /*
     G4int GetNbOfLayers()              {return NbOfLayers;}; 
     
     G4Material* GetAbsorberMaterial()  {return AbsorberMaterial;};
     G4double    GetAbsorberThickness() {return AbsorberThickness;}; 

     G4Material* GetAbsorber2Material() {return Absorber2Material;};
     G4double    GetAbsorber2Thickness() {return Absorber2Thickness;};     
     
     G4Material* GetAbsorber3Material()  {return Absorber3Material;};
     G4double    GetAbsorber3Thickness() {return Absorber3Thickness;};

     G4Material* GetGapMaterial()       {return GapMaterial;};
     G4double    GetGapThickness()      {return GapThickness;};

     G4Material* GetGap2Material()      {return Gap2Material;};
     G4double    GetGap2Thickness()     {return Gap2Thickness;};

     G4Material* GetGap3Material()      {return Gap3Material;};
     G4double    GetGap3Thickness()     {return Gap3Thickness;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
     const G4VPhysicalVolume* GetAbsorber()   {return physiAbsorber;};
     const G4VPhysicalVolume* GetAbsorber2()  {return physiAbsorber2;};
     const G4VPhysicalVolume* GetAbsorber3()  {return physiAbsorber3;};
     const G4VPhysicalVolume* GetGap()        {return physiGap;};
     const G4VPhysicalVolume* GetGap2()       {return physiGap2;};
     const G4VPhysicalVolume* GetGap3()       {return physiGap3;};
     */
                 
  private:
     /*
     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;

     G4Material*        Absorber2Material;
     G4double           Absorber2Thickness;

     G4Material*	Absorber3Material;
     G4double		Absorber3Thickness;
     
     G4Material*        GapMaterial;
     G4double           GapThickness;
     
     G4Material*        Gap2Material;
     G4double           Gap2Thickness;

     G4Material*	Gap3Material;
     G4double		Gap3Thickness;

     G4int              NbOfLayers;
     G4double           LayerThickness;
     G4double           Layer2Thickness;
     G4double		Layer3Thickness;
     */
          
     G4double           CalorSizeYZ;
     //G4double           CalorThickness;
     
     
     G4Material*        defaultMaterial;
     G4Material*        Iron;
     G4Material*        Tungsten;
     G4Material*        Plexiglas;
     G4Material*        Hydrogen;

     G4double           WorldSizeYZ;
     G4double           WorldSizeX;
            
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

     /*
     G4Box*             solidCalor;    //pointer to the solid Calor 
     G4LogicalVolume*   logicCalor;    //pointer to the logical Calor
     G4VPhysicalVolume* physiCalor;    //pointer to the physical Calor
     
     G4Box*             solidLayer;    //pointer to the solid Layer 
     G4LogicalVolume*   logicLayer;    //pointer to the logical Layer
     G4VPhysicalVolume* physiLayer;    //pointer to the physical Layer

     G4Box*             solidLayer2;
     G4LogicalVolume*   logicLayer2;
     G4VPhysicalVolume* physiLayer2;

     G4Box*		solidLayer3;
     G4LogicalVolume*	logicLayer3;
     G4VPhysicalVolume*	physiLayer3;
         
     G4Box*             solidAbsorber; //pointer to the solid Absorber
     G4LogicalVolume*   logicAbsorber; //pointer to the logical Absorber
     G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber

     G4Box*             solidAbsorber2;
     G4LogicalVolume*   logicAbsorber2;
     G4VPhysicalVolume* physiAbsorber2;

     G4Box*		solidAbsorber3;
     G4LogicalVolume*	logicAbsorber3;
     G4VPhysicalVolume*	physiAbsorber3;
     
     G4Box*             solidGap;      //pointer to the solid Gap
     G4LogicalVolume*   logicGap;      //pointer to the logical Gap
     G4VPhysicalVolume* physiGap;      //pointer to the physical Gap

     G4Box*             solidGap2;
     G4LogicalVolume*   logicGap2;
     G4VPhysicalVolume* physiGap2;

     G4Box*		solidGap3;
     G4LogicalVolume*	logicGap3;
     G4VPhysicalVolume*	physiGap3;
     */
     
     //G4UniformMagField* magField;      //pointer to the magnetic field
     
     DetectorMessenger* detectorMessenger;  //pointer to the Messenger
      
  private:
    
     void DefineMaterials();
     //void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*inline void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
     LayerThickness = AbsorberThickness + GapThickness;
     Layer2Thickness = Absorber2Thickness + Gap2Thickness;
     Layer3Thickness = Absorber3Thickness + Gap3Thickness;
  if(Layer2Thickness > 0 and Layer3Thickness > 0)
     CalorThickness = NbOfLayers*LayerThickness + Layer2Thickness + Layer3Thickness;
  else
     CalorThickness = NbOfLayers*LayerThickness;
     
     WorldSizeX = 1.2*CalorThickness; WorldSizeYZ = 1.2*CalorSizeYZ;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

