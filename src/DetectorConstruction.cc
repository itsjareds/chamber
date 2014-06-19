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

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4SubtractionSolid.hh"
#include "G4EllipticalCone.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:defaultMaterial(0),
 solidWorld(0),logicWorld(0),physiWorld(0),
 WorldSizeX(3.0*m), WorldSizeYZ(3.0*m),
 CalorSizeYZ(107.6352) //probably wrong--fix later
{
  // default parameter values of the calorimeter
  /*AbsorberThickness = 250.*um;
  GapThickness      =  25.*um;
  Absorber2Thickness=  10.*um;
  Gap2Thickness     =  10.*um;
  Absorber3Thickness=   2.*mm;
  Gap3Thickness     =  31.*mm;
  NbOfLayers        =   1;
  CalorSizeYZ       =  10.*cm;*/
  //ComputeCalorParameters();
  
  // materials
  DefineMaterials();
/*SetAbsorberMaterial("Beryllium");
  SetGapMaterial("Air");
  SetAbsorber2Material("Mylar");
  SetGap2Material("Water");
  SetAbsorber3Material("Water");
  SetGap3Material("Air");*/
 
  // create commands for interactive definition of the calorimeter
  detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String symbol;             //a=mass of a mole;
G4double a, z, density;      //z=mean number of protons;  
G4int iz, n;                 //iz=number of protons  in an isotope; 
                             // n=number of nucleons in an isotope;

G4int ncomponents, natoms;
G4double abundance, fractionmass;

//
// define Elements
//

G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
G4Element* Si = new G4Element("Silicon",symbol="Si" , z= 14., a= 28.09*g/mole);
G4Element* Be = new G4Element("Beryllium", symbol="Be" , z=04., a=9.012*g/mole);
/* weights taken from Wikipedia from here on down */
G4Element* S  = new G4Element("Sulfur", symbol="S", z=16., a=32.07*g/mole);
G4Element* Na = new G4Element("Sodium", symbol="Na", z=11., a=22.99*g/mole);
G4Element* Cl = new G4Element("Chlorine", symbol="Cl", z=17., a=35.45*g/mole);
G4Element* P  = new G4Element("Phosphorus", symbol="P", z=15., a=30.97*g/mole);
G4Element* K  = new G4Element("Potassium", symbol="K", z=19., a=39.10*g/mole);
G4Element* Fe = new G4Element("Iron", symbol="Fe", z=26., a=55.845*g/mole);
G4Element* Mg = new G4Element("Magnesium", symbol="Mg", z=12., a=24.305*g/mole);
G4Element* Ca = new G4Element("Calcium", symbol="Ca", z=20., a=40.078*g/mole);
G4Element* W  = new G4Element("Tungsten", symbol="W", z=74., a=183.84*g/mole);

//
// define an Element from isotopes, by relative abundance 
//

G4Isotope* U5 = new G4Isotope("U235", iz=92, n=235, a=235.01*g/mole);
G4Isotope* U8 = new G4Isotope("U238", iz=92, n=238, a=238.03*g/mole);

G4Element* U  = new G4Element("enriched Uranium",symbol="U",ncomponents=2);
U->AddIsotope(U5, abundance= 90.*perCent);
U->AddIsotope(U8, abundance= 10.*perCent);

//
// define simple materials
//
/*
new G4Material("Aluminium", z=13., a=26.98*g/mole, density=2.700*g/cm3);
new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
new G4Material("Lead"     , z=82., a= 207.19*g/mole, density= 11.35*g/cm3);
*/
//
// define a material from elements.   case 1: chemical molecule
//

G4Material* H2O = 
new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
H2O->AddElement(H, natoms=2);
H2O->AddElement(O, natoms=1);
// overwrite computed meanExcitationEnergy with ICRU recommended value 
H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

G4Material* Sci = 
new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
Sci->AddElement(C, natoms=9);
Sci->AddElement(H, natoms=10);

G4Material* Myl = 
new G4Material("Mylar", density= 1.397*g/cm3, ncomponents=3);
Myl->AddElement(C, natoms=10);
Myl->AddElement(H, natoms= 8);
Myl->AddElement(O, natoms= 4);

G4Material* SiO2 = 
new G4Material("quartz",density= 2.200*g/cm3, ncomponents=2);
SiO2->AddElement(Si, natoms=1);
SiO2->AddElement(O , natoms=2);

G4Material* Plastic = 
new G4Material("Plastic", density= 1*g/cm3, ncomponents=2);
Plastic->AddElement(H, natoms=8);
Plastic->AddElement(C, natoms=8);

G4Material* Beryllium =
new G4Material("Beryllium", density=1.85*g/cm3, ncomponents=1);
Beryllium->AddElement(Be, fractionmass=1);

Iron =
new G4Material("Iron", density=7.874*g/cm3, ncomponents=1);
Iron->AddElement(Fe, fractionmass=1);

Tungsten =
new G4Material("Tungsten", density=19.25*g/cm3, ncomponents=1);
Tungsten->AddElement(W, fractionmass=1);

Hydrogen =
new G4Material("Hydrogen", density=0.00008988*g/cm3, ncomponents=1);
Hydrogen->AddElement(H, fractionmass=1);

//
// define a material from elements.   case 2: mixture by fractional mass
//

G4Material* Air = 
new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
Air->AddElement(N, fractionmass=0.7);
Air->AddElement(O, fractionmass=0.3);

G4Material* Air2 = 
new G4Material("Air2"  , density= 1.290*mg/cm3, ncomponents=2);
Air2->AddElement(N, fractionmass=0.78);
Air2->AddElement(O, fractionmass=0.22);

Plexiglas =
new G4Material("Plexiglas", density= 1.19*g/cm3, ncomponents=3);
Plexiglas->AddElement(H, fractionmass=0.0805);
Plexiglas->AddElement(C, fractionmass=0.5999);
Plexiglas->AddElement(O, fractionmass=0.3196);

G4Material* Fat = 
new G4Material("icru-44 fat", density= 0.95*gram/cm3, ncomponents=7);
Fat->AddElement(H, fractionmass=0.114);
Fat->AddElement(C, fractionmass=0.598);
Fat->AddElement(N, fractionmass=0.007);
Fat->AddElement(O, fractionmass=0.278);
Fat->AddElement(Na, fractionmass=0.001);
Fat->AddElement(S, fractionmass=0.001);
Fat->AddElement(Cl, fractionmass=0.001);
/*		("H", 0.114), ("C", 0.598), ("N", 0.007), ("O", 0.278),  
		("Na", 0.001), ("S", 0.001), ("Cl", 0.001) */
	
G4Material* Blood =
new G4Material("icru-44 blood", density= 1.06*gram/cm3, ncomponents=10);
Blood->AddElement(H, fractionmass=0.102);
Blood->AddElement(C, fractionmass=0.110);
Blood->AddElement(N, fractionmass=0.033);
Blood->AddElement(O, fractionmass=0.745);
Blood->AddElement(Na, fractionmass=0.001);
Blood->AddElement(P, fractionmass=0.001);
Blood->AddElement(S, fractionmass=0.002);
Blood->AddElement(Cl, fractionmass=0.003);
Blood->AddElement(K, fractionmass=0.002);
Blood->AddElement(Fe, fractionmass=0.001);
/*		("H", 0.102), ("C", 0.110), ("N", 0.033), ("O", 0.745),  
		("Na", 0.001), ("P", 0.001), ("S", 0.002), ("Cl", 0.003),
		("K", 0.002), ("Fe", 0.001) */

G4Material* Muscle =
new G4Material("icru-44 skeletal muscle", density= 1.05*gram/cm3, ncomponents=9);
Muscle->AddElement(H, fractionmass=0.102);
Muscle->AddElement(C, fractionmass=0.143);
Muscle->AddElement(N, fractionmass=0.034);
Muscle->AddElement(O, fractionmass=0.710);
Muscle->AddElement(Na, fractionmass=0.001);
Muscle->AddElement(P, fractionmass=0.002);
Muscle->AddElement(S, fractionmass=0.003);
Muscle->AddElement(Cl, fractionmass=0.001);
Muscle->AddElement(K, fractionmass=0.004);
/*		("H", 0.102), ("C", 0.143), ("N", 0.034), ("O", 0.710),  
		("Na", 0.001), ("P", 0.002), ("S", 0.003), ("Cl", 0.001),
		("K", 0.004) */

G4Material* Tissue =
new G4Material("icru-44 soft tissue", density= 1.06*gram/cm3, ncomponents=9);
Tissue->AddElement(H, fractionmass=0.102);
Tissue->AddElement(C, fractionmass=0.143);
Tissue->AddElement(N, fractionmass=0.034);
Tissue->AddElement(O, fractionmass=0.708);
Tissue->AddElement(Na, fractionmass=0.002);
Tissue->AddElement(P, fractionmass=0.003);
Tissue->AddElement(S, fractionmass=0.003);
Tissue->AddElement(Cl, fractionmass=0.002);
Tissue->AddElement(K, fractionmass=0.003);
/*		("H", 0.102), ("C", 0.143), ("N", 0.034), ("O", 0.708),  
		("Na", 0.002), ("P", 0.003), ("S", 0.003), ("Cl", 0.002),
		("K", 0.003) */

G4Material* Lung =
new G4Material("icru-44 lung tissue", density= 0.26*gram/cm3, ncomponents=9);
Lung->AddElement(H, fractionmass=0.103);
Lung->AddElement(C, fractionmass=0.105);
Lung->AddElement(N, fractionmass=0.031);
Lung->AddElement(O, fractionmass=0.749);
Lung->AddElement(Na, fractionmass=0.002);
Lung->AddElement(P, fractionmass=0.002);
Lung->AddElement(S, fractionmass=0.003);
Lung->AddElement(Cl, fractionmass=0.003);
Lung->AddElement(K, fractionmass=0.002);
/*		("H", 0.103), ("C", 0.105), ("N", 0.031), ("O", 0.749),  
		("Na", 0.002), ("P", 0.002), ("S", 0.003), ("Cl", 0.003),
		("K", 0.002) */

G4Material* Skin =
new G4Material("icru-44 skin", density= 1.09*gram/cm3, ncomponents=9);
Skin->AddElement(H, fractionmass=0.100);
Skin->AddElement(C, fractionmass=0.204);
Skin->AddElement(N, fractionmass=0.042);
Skin->AddElement(O, fractionmass=0.645);
Skin->AddElement(Na, fractionmass=0.002);
Skin->AddElement(P, fractionmass=0.001);
Skin->AddElement(S, fractionmass=0.002);
Skin->AddElement(Cl, fractionmass=0.003);
Skin->AddElement(K, fractionmass=0.001);
/*		("H", 0.100), ("C", 0.204), ("N", 0.042), ("O", 0.645),  
		("Na", 0.002), ("P", 0.001), ("S", 0.002), ("Cl", 0.003),
		("K", 0.001) */

G4Material* Cartilage =
new G4Material("icru-44 cartilage", density= 1.10*gram/cm3, ncomponents=8);
Cartilage->AddElement(H, fractionmass=0.096);
Cartilage->AddElement(C, fractionmass=0.099);
Cartilage->AddElement(N, fractionmass=0.022);
Cartilage->AddElement(O, fractionmass=0.744);
Cartilage->AddElement(Na, fractionmass=0.005);
Cartilage->AddElement(P, fractionmass=0.022);
Cartilage->AddElement(S, fractionmass=0.009);
Cartilage->AddElement(Cl, fractionmass=0.003);
/*		("H", 0.096), ("C", 0.099), ("N", 0.022), ("O", 0.744),  
		("Na", 0.005), ("P", 0.022), ("S", 0.009), ("Cl", 0.003), */

G4Material* Nerve =
new G4Material("icru-44 brain", density= 1.04*gram/cm3, ncomponents=9);
Nerve->AddElement(H, fractionmass=0.107);
Nerve->AddElement(C, fractionmass=0.145);
Nerve->AddElement(N, fractionmass=0.022);
Nerve->AddElement(O, fractionmass=0.712);
Nerve->AddElement(Na, fractionmass=0.002);
Nerve->AddElement(P, fractionmass=0.004);
Nerve->AddElement(S, fractionmass=0.002);
Nerve->AddElement(Cl, fractionmass=0.003);
Nerve->AddElement(K, fractionmass=0.003);
/*		("H", 0.107), ("C", 0.145), ("N", 0.022), ("O", 0.712),  
		("Na", 0.002), ("P", 0.004), ("S", 0.002), ("Cl", 0.003),
		("K", 0.003) */

G4Material* Bone =
new G4Material("icru-44 bone", density= 1.92*gram/cm3, ncomponents=9);
Bone->AddElement(H, fractionmass=0.034);
Bone->AddElement(C, fractionmass=0.155);
Bone->AddElement(N, fractionmass=0.042);
Bone->AddElement(O, fractionmass=0.435);
Bone->AddElement(Na, fractionmass=0.001);
Bone->AddElement(Mg, fractionmass=0.002);
Bone->AddElement(P, fractionmass=0.103);
Bone->AddElement(S, fractionmass=0.003);
Bone->AddElement(Ca, fractionmass=0.225);
/*		("H", 0.034), ("C", 0.155), ("N", 0.042), ("O", 0.435),
		("Na", 0.001), ("Mg", 0.002), ("P", 0.103), ("S", 0.003),
		("Ca", 0.225) */

//
// define a material from elements and/or others materials (mixture of mixtures)
//

G4Material* Aerog = 
new G4Material("Aerogel", density= 0.200*g/cm3, ncomponents=3);
Aerog->AddMaterial(SiO2, fractionmass=62.5*perCent);
Aerog->AddMaterial(H2O , fractionmass=37.4*perCent);
Aerog->AddElement (C   , fractionmass= 0.1*perCent);

//
// examples of gas in non STP conditions
//

G4Material* CO2 = 
new G4Material("CarbonicGas", density= 1.842*mg/cm3, ncomponents=2,
                              kStateGas, 325.*kelvin, 50.*atmosphere);
CO2->AddElement(C, natoms=1);
CO2->AddElement(O, natoms=2);
 
G4Material* steam = 
new G4Material("WaterSteam", density= 0.3*mg/cm3, ncomponents=1,
                             kStateGas, 500.*kelvin, 2.*atmosphere);
steam->AddMaterial(H2O, fractionmass=1.);

//
// examples of vacuum
//

G4Material* Vacuum =
new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                           kStateGas, 2.73*kelvin, 3.e-18*pascal);

G4Material* beam = 
new G4Material("Beam", density= 1.e-5*g/cm3, ncomponents=1,
                       kStateGas, STP_Temperature, 2.e-2*bar);
beam->AddMaterial(Air, fractionmass=1.);

//
// or use G4-NIST materials data base
//
G4NistManager* man = G4NistManager::Instance();
man->FindOrBuildMaterial("G4_SODIUM_IODIDE");

// print table
//
G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//default materials of the World
defaultMaterial  = Vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{
  // Variables copied from fixedbeam30.py
  G4double worldFullLength = 3.0*m;

  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the Calor parameters definition
  //ComputeCalorParameters();
   
  //     
  // World
  //
  solidWorld = new G4Box("World",                             //its name
                   worldFullLength/2,worldFullLength/2,worldFullLength/2); //its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,           //its solid
                                   defaultMaterial,      //its material
                                   "World");             //its name
                                   
  physiWorld = new G4PVPlacement(0,                      //no rotation
                                 G4ThreeVector(),        //at (0,0,0)
                                 logicWorld,             //its logical volume                                 
                                 "World",                //its name
                                 0,                      //its mother  volume
                                 false,                  //no boolean operation
                                 0);                     //copy number

 // Pasted from fixedbeam30.py
		
		int firstlayer=0;
		int lastlayer=303; //#255
		//#dataset_shape=(128,256,256)
		std::vector<int> dataset_shape;
		dataset_shape.push_back(512);
		dataset_shape.push_back(512);
		dataset_shape.push_back(303);

		std::vector<int> keep_organs;
		keep_organs.push_back(4);
		//#keep_organs=[0, 2, 4]
				
		//uncomment later maybe
		//detailedLayers=None
		G4double visFraction=0.0001;
		G4double x_pixel_size=0.625*mm;
		G4double y_pixel_size=0.625*mm;
		G4double z_pixel_size=0.625*mm;
		G4bool useRLEbyOrgan=false;
		G4bool useRLEbymaterial=true;
		
		G4String voxelFile="ct_phantom.dat";

		
		G4String tissueMapFile="zubal_head_map.txt";
/*		
		quads, material_key=read_zubal_key(tissueMapFile)
		data=smoother_debrecen.voxel_file(filename=voxelFile, shape=dataset_shape,
			deltas=(x_pixel_size, y_pixel_size, z_pixel_size) )

		materials_info=create_zubal_materials()
*/
		
		G4double zcount=lastlayer-firstlayer;
		int xcount = dataset_shape[0];
		int ycount = dataset_shape[1];
		int junk = dataset_shape[2];
	
		G4Colour skin(1,0,0);
		G4Colour fat(0,1,0);
		G4Colour skull(0,0,1);
		G4Colour jawbone(1,1,0);
		G4Colour cerebral_fluid(0,1,1);

		// Uncomment later?
		//print xcount, ycount, zcount, Vector((x_pixel_size, y_pixel_size, z_pixel_size))
/*
		
		if not orbOnly:
			compressor=G4Core.mhmCompressedVoxels(
				xcount, ycount, zcount, 
				Vector((x_pixel_size, y_pixel_size, z_pixel_size)),
				visFraction)
					
			//colors={	'skin':G4Core.G4Colour(1,0,0, 0.5), 
			//			'fat':G4Core.G4Colour(0,1,0, 0.5), 
			//			'skull':G4Core.G4Colour(0,0,1,1), 
			//			'jaw bone':G4Core.G4Colour(1,1,0,1), 
			//			'cerebral fluid': G4Core.G4Colour(0,1,1,1)
			//}

			vis_attr_list=[]
			faded=G4Core.G4Colour(1,1,1,1)
			invisible=G4Core.G4VisAttributes(0, faded)

			for junk, orgidx, name, matidx in quads:
				compressor.SetOrganMaterial(orgidx, materials_info.tissuemap[matidx])

				if (orgidx == 0) {
					v=invisible //#organ 0 is always invisible
				} else {
					color=colors.get(name, None)
					if color is None:
						color=G4Core.G4Colour(
							0.5*random.random()+0.25, 
							0.5*random.random()+0.25,
							0.5*random.random()+0.25, 1)
					v=G4Core.G4VisAttributes(1, color)
					v.SetForceSolid(1)
					v.SetLineWidth(0.0)
				}

				vis_attr_list.push_back(v); //#keep ownership since we didn't disown these
				compressor.SetOrganVisAttributes(orgidx, v)
			
			logic_vol=compressor.AssembleGeometry(physiWorld, "zubal", materials_info.air)
						
			G4PVPlacement* phantomPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(0.,0.,0.), logic_vol, "phantom", parent=logicWorld, false, 0);

			for orgidx in keep_organs:
				compressor.AddKeptOrgan(orgidx)
				compressor.SetOrganVisSampling(orgidx,-1) //# < 0 forces visibility

			for (int iz = firstlayer; iz < lastlayer; ++iz) {
				detailed = (detailedLayers is None or iz in detailedLayers)
				compressor.CompressAndAppendVoxelPlane(
					data.get_plane(iz), detailed,  useRLEbyOrgan, useRLEbyMaterial
				)
			}
			
			compressor.BuildContainerSolid(phantomPhys);
			
		else:
			G4Orb* phantomSolid = new G4Orb("sphere", 100.0*mm)
			G4LogicalVolume* phantomLog = new G4LogicalVolume(phantomSolid, Tissue, "orb");
			phantomPhys=Placement(phantomLog, 'orb', parent=logicWorld, pos=(0,0,0))
			phantomLog.SetForceWireFrame(True)	
			phantomLog.vis.SetForceAuxEdgeVisible(True)	
*/
		// Uncomment later?				
		G4cout << "Phantom fully populated\n";
		

		//#########################
		//#### source construction
		//#########################
		

		//uncomment later
		//sourceAssembly=G4Support.AssemblyVolume()

		//####container
		G4double containerOuterDiam= 107.6352;
		G4double containerHeight = 50.0126 + 51.5*2 + 89 + 1.5;
		G4double containerPosX=0.0;
		G4double containerPosY=0.0;
		G4double containerPosZ=0.0;
		
		G4Tubs* containerSolid = new G4Tubs("containerSolid", 0, containerOuterDiam/2, containerHeight/2, 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* container = new G4LogicalVolume(containerSolid, /*material=*/Hydrogen, /*name=*/"container");
				//, /*color=*/jawbone);
		
		G4RotationMatrix* zero_rotation = new G4RotationMatrix();
		//uncomment later
		G4PVPlacement* containerPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(containerPosX,containerPosY,containerPosZ), container, "containerPhys", logicWorld, false, 0);

		//####largest tube

		G4double largestTubeInnerDiam=12.1412;
		G4double largestTubeOuterDiam=27.6351;
		G4double largestTubeHeight = 50.0126;
		G4double largestTubePosX=0.0;
		G4double largestTubePosY=0.0;
		G4double largestTubePosZ=largestTubeHeight/2 - containerHeight/2;

		G4Tubs* largestTube = new G4Tubs("largestTube", largestTubeInnerDiam/2, largestTubeOuterDiam/2, largestTubeHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* largestTubeLogic = new G4LogicalVolume(largestTube, /*material=*/Iron, /*name=*/"largestTubeLogic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* largestTubePhys = new G4PVPlacement(zero_rotation, G4ThreeVector(largestTubePosX,largestTubePosY,largestTubePosZ), largestTubeLogic, "largestTubePhys", container, false, 0);
		//#largestTubeLogic.SetForceWireFrame(1)	

		//#sourceAssembly.AddPlacedVolume( largestTubeLogic, pos=(largestTubePosX,largestTubePosY,largestTubePosZ) )  

		//####middle tube

		G4double middleTubeInnerDiam=7.112;
		G4double middleTubeOuterDiam=11.6078;
		G4double middleTubeHeight=36.4236;
		G4double middleTubePosX=0.0;
		G4double middleTubePosY=0.0;
		G4double middleTubePosZ= ( (largestTubeHeight-middleTubeHeight)/2 - 1.4224 + largestTubeHeight/2 ) - containerHeight/2;

		G4Tubs* middleTube = new G4Tubs("middleTube", middleTubeInnerDiam/2, middleTubeOuterDiam/2, middleTubeHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* middleTubeLogic = new G4LogicalVolume(middleTube, /*material=*/Iron, /*name=*/"middleTubeLogic");
				//, /*color=*/fat);
		//uncomment later
		G4PVPlacement* middleTubePhys = new G4PVPlacement(zero_rotation, G4ThreeVector(middleTubePosX,middleTubePosY,middleTubePosZ), middleTubeLogic, "middleTubePhys", container, false, 0);
		//middleTubeLogic.SetForceWireFrame(1)

		//#sourceAssembly.AddPlacedVolume( middleTubeLogic, pos=(middleTubePosX,middleTubePosY,middleTubePosZ) )  

		//####smallest tube

		G4double smallestTubeInnerDiam=2.5908;
		G4double smallestTubeOuterDiam=6.604;
		G4double smallestTubeHeight=32.9438;
		G4double smallestTubePosX=0.0;
		G4double smallestTubePosY=0.0;
		G4double smallestTubePosZ=( (largestTubeHeight-smallestTubeHeight) /2 - 1.4224 -0.4064 + largestTubeHeight/2 ) - containerHeight/2;

		G4Tubs* smallestTube = new G4Tubs("smallestTube", smallestTubeInnerDiam/2, smallestTubeOuterDiam/2, smallestTubeHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* smallestTubeLogic = new G4LogicalVolume(smallestTube, /*material=*/Iron, /*name=*/"smallestTubeLogic");
				//, /*color=*/skull);
		//uncomment later
		G4PVPlacement* smallestTubePhys = new G4PVPlacement(zero_rotation, G4ThreeVector(smallestTubePosX,smallestTubePosY, smallestTubePosZ), smallestTubeLogic, "smallestTubePhys", container, false, 0);
		//smallestTubeLogic.SetForceWireFrame(1)

		//#sourceAssembly.AddPlacedVolume( smallestTubeLogic, pos=(smallestTubePosX,smallestTubePosY, smallestTubePosZ) )  

		//####cap of the largest tube  
 
		G4double largestCapInnerDiam=8;
		//#largestCapOuterDiam=27.6351
		G4double largestCapOuterDiam=middleTubeOuterDiam;
		G4double largestCapHeight=1.4224;
		G4double largestCapPosX=0.0;
		G4double largestCapPosY=0.0;
		G4double largestCapPosZ=(largestTubeHeight/2 - 1.4224/2 + largestTubeHeight/2) - containerHeight/2;

		G4Tubs* largestCap = new G4Tubs("largestCap", largestCapInnerDiam/2, largestCapOuterDiam/2, largestCapHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* largestCapLogic = new G4LogicalVolume(largestCap, /*material=*/Iron, /*name=*/"largestCapLogic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* largestCapPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(largestCapPosX,largestCapPosY,largestCapPosZ), largestCapLogic, "largestCapPhys", container, false, 0);
		//largestCapLogic.SetForceWireFrame(1)	

		//#sourceAssembly.AddPlacedVolume( largestCapLogic, pos=(largestCapPosX,largestCapPosY,largestCapPosZ) )  


		//####cap of the middle tube 

		G4double middleCapInnerDiam=0;
		//#middleCapOuterDiam=11.6078;
		G4double middleCapOuterDiam=smallestTubeOuterDiam;
		G4double middleCapHeight=0.4064;
		G4double middleCapPosX=0.0;
		G4double middleCapPosY=0.0;
		G4double middleCapPosZ=( (largestTubeHeight -middleCapHeight)/2 - largestCapHeight + largestTubeHeight/2 ) - containerHeight/2;

		G4Tubs* middleCap = new G4Tubs("middleCap", middleCapInnerDiam/2, middleCapOuterDiam/2, middleCapHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* middleCapLogic = new G4LogicalVolume(middleCap, /*material=*/Iron, /*name=*/"middleCapLogic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* middleCapPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(middleCapPosX,middleCapPosY,middleCapPosZ), middleCapLogic, "middleCapPhys", container, false, 0);
		//middleCapLogic.SetForceWireFrame(1)	

		//#sourceAssembly.AddPlacedVolume( middleCapLogic, pos=(middleCapPosX,middleCapPosY,middleCapPosZ) )  

		//####cap of the smallest tube 
		
		G4double smallestCapInnerDiam=0;
		G4double smallestCapOuterDiam=2.5908;
		G4double smallestCapHeight=0.685;
		G4double smallestCapPosX=0.0;
		G4double smallestCapPosY=0.0;
		G4double smallestCapPosZ=( (50.0126 - 32.9438 ) /2 - 1.4224 -0.4064  + (smallestTubeHeight - smallestCapHeight )/2 + largestTubeHeight/2 ) - containerHeight/2;

		G4Tubs* smallestCap = new G4Tubs("smallestCap", smallestCapInnerDiam/2, smallestCapOuterDiam/2, smallestCapHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4EllipticalCone* smallestCapCone = new G4EllipticalCone( "smallestCapCone", smallestCapOuterDiam/2, smallestCapOuterDiam/2, 0.7/2, 0.7/2);

		//uncomment later
		G4SubtractionSolid* smallestCapSolid = new G4SubtractionSolid("smallestCap-smallestCapCone", smallestCap, smallestCapCone);

		//uncomment later
		G4LogicalVolume* smallestCapLogic = new G4LogicalVolume(smallestCapSolid, /*material=*/Iron, /*name=*/"smallestCapLogic");
				//, /*color=*/cerebral_fluid);
		//uncomment later
		G4PVPlacement* smallestCapPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(smallestCapPosX,smallestCapPosY,smallestCapPosZ), smallestCapLogic, "smallestCapPhys", container, false, 0);
		//smallestCapLogic.SetForceSolid(1)	

		//#sourceAssembly.AddPlacedVolume( smallestCapLogic, pos=(smallestCapPosX,smallestCapPosY,smallestCapPosZ) )  


		//####backcap of the largest tube 
		
		G4double largestBackCapInnerDiam=0;
		//#largestBackCapOuterDiam=27.6352;
		G4double largestBackCapOuterDiam=middleTubeOuterDiam;
		G4double largestBackCapHeight=2.2606;
		G4double largestBackCapPosX=0.0;
		G4double largestBackCapPosY=0.0;
		G4double largestBackCapPosZ= (-largestTubeHeight/2 + largestBackCapHeight/2  + largestTubeHeight/2) -containerHeight/2;

		G4Tubs* largestBackCap = new G4Tubs("largestBackCap", largestBackCapInnerDiam/2, largestBackCapOuterDiam/2, largestBackCapHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* largestBackCapLogic = new G4LogicalVolume(largestBackCap, /*material=*/Iron, /*name=*/"largestBackCapLogic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* largestBackCapPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(largestBackCapPosX,largestBackCapPosY,largestBackCapPosZ), largestBackCapLogic, "largestBackCapPhys", container, false, 0);
		//largestBackCapLogic.SetForceWireFrame(1)	

		//#sourceAssembly.AddPlacedVolume( largestBackCapLogic, pos=(largestBackCapPosX,largestBackCapPosY,largestBackCapPosZ) )  

		//####backcap of the middle tube 
		
		G4double middleBackCapInnerDiam=0;
		//#middleBackCapOuterDiam=11.6078
		G4double middleBackCapOuterDiam=smallestTubeOuterDiam;
		G4double middleBackCapHeight=2.0574;
		G4double middleBackCapPosX=0.0;
		G4double middleBackCapPosY=0.0;
		G4double middleBackCapPosZ=( middleTubePosZ-middleTubeHeight/2 + middleBackCapHeight/2 );

		G4Tubs* middleBackCap = new G4Tubs("middleBackCap", middleBackCapInnerDiam/2, middleBackCapOuterDiam/2, middleBackCapHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* middleBackCapLogic = new G4LogicalVolume(middleBackCap, /*material=*/Iron, /*name=*/"middleBackCapLogic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* middleBackCapPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(middleBackCapPosX,middleBackCapPosY,middleBackCapPosZ), middleBackCapLogic, "middleBackCapPhys", container, false, 0);
		//middleBackCapLogic.SetForceWireFrame(1)	

		//#sourceAssembly.AddPlacedVolume( middleBackCapLogic, pos=(middleBackCapPosX,middleBackCapPosY,middleBackCapPosZ) )  

		//####backcap of the smallest tube 
		
		G4double smallestBackCapInnerDiam=0;
		//#smallestBackCapOuterDiam=6.64
		G4double smallestBackCapOuterDiam=smallestTubeInnerDiam;
		G4double smallestBackCapHeight=2;
		G4double smallestBackCapPosX=0.0;
		G4double smallestBackCapPosY=0.0;
		G4double smallestBackCapPosZ=( smallestTubePosZ - smallestTubeHeight/2 + smallestBackCapHeight/2 );

		G4Tubs* smallestBackCap = new G4Tubs("smallestBackCap", smallestBackCapInnerDiam/2, smallestBackCapOuterDiam/2, smallestBackCapHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* smallestBackCapLogic = new G4LogicalVolume(smallestBackCap, /*material=*/Iron, /*name=*/"smallestBackCapLogic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* smallestBackCapPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(smallestBackCapPosX,smallestBackCapPosY,smallestBackCapPosZ), smallestBackCapLogic, "smallestBackCapPhys", container, false, 0);
		//uncomment later
		//smallestBackCapLogic.SetForceWireFrame(1)	

		//#sourceAssembly.AddPlacedVolume( smallestBackCapLogic, pos=(smallestBackCapPosX,smallestBackCapPosY,smallestBackCapPosZ) )  

		//####Cobalt filling
 		
		G4double cobaltInnerDiam=0;
		G4double cobaltOuterDiam=smallestTubeInnerDiam;
		G4double cobaltHeight=smallestTubeHeight - smallestBackCapHeight - smallestCapHeight;
		G4double cobaltPosX=0.0;
		G4double cobaltPosY=0.0;
		G4double cobaltPosZ=( smallestTubePosZ + (smallestTubeHeight - cobaltHeight )/2 - smallestCapHeight  );
		

		G4Tubs* cobalt = new G4Tubs("cobalt", cobaltInnerDiam/2, cobaltOuterDiam/2, cobaltHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* cobaltLogic = new G4LogicalVolume(cobalt, /*material=*/Iron, /*name=*/"cobaltLogic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* cobaltPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(cobaltPosX,cobaltPosY,cobaltPosZ), cobaltLogic, "cobaltPhys", container, false, 0);
		//cobaltLogic.SetForceWireFrame(1)	

		//#sourceAssembly.AddPlacedVolume( cobaltLogic, pos=(cobaltPosX,cobaltPosY,cobaltPosZ) )  


		//####Tungsten plug
		
		G4double tungstenPlugInnerDiam=0;
		G4double tungstenPlugOuterDiam=12.0142;
		G4double tungstenPlugHeight=9.271;
		G4double tungstenPlugPosX=0.0;
		G4double tungstenPlugPosY=0.0;
		G4double tungstenPlugPosZ=(middleTubePosZ-middleTubeHeight/2 -tungstenPlugHeight/2 );

		G4Tubs* tungstenPlug = new G4Tubs("tungstenPlug", tungstenPlugInnerDiam/2, tungstenPlugOuterDiam/2, tungstenPlugHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* tungstenPlugLogic = new G4LogicalVolume(tungstenPlug, /*material=*/Tungsten, /*name=*/"tungstenPlugLogic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* tungstenPlugPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(tungstenPlugPosX,tungstenPlugPosY,tungstenPlugPosZ), tungstenPlugLogic, "tungstenPlugPhys", container, false, 0);
		//tungstenPlugLogic.SetForceWireFrame(1)	

		//#sourceAssembly.AddPlacedVolume( tungstenPlugLogic, pos=(tungstenPlugPosX,tungstenPlugPosY,tungstenPlugPosZ) )  


		//####primary coll
		
		G4double primaryCollInnerDiam=6.5;
		G4double primaryCollOuterDiam=25;
		G4double primaryCollHeight=51.5;
		G4double primaryCollPosX=0.0;
		G4double primaryCollPosY=0.0;
		G4double primaryCollPosZ= (50.0126/2 + primaryCollHeight/2  + largestTubeHeight/2 ) - containerHeight/2;

		G4Tubs* primaryColl = new G4Tubs("primaryColl", primaryCollInnerDiam/2, primaryCollOuterDiam/2, primaryCollHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* primaryCollLogic = new G4LogicalVolume(primaryColl, /*material=*/Iron, /*name=*/"primaryCollLogic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* primaryCollPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(primaryCollPosX,primaryCollPosY,primaryCollPosZ), primaryCollLogic, "primaryCollPhys", container, false, 0);
		//primaryCollLogic.SetForceWireFrame(1)	

		//#sourceAssembly.AddPlacedVolume( primaryCollLogic, pos=(primaryCollPosX,primaryCollPosY,primaryCollPosZ) )  
		
		G4double primaryColl2InnerDiam=6.5;
		G4double primaryColl2OuterDiam=25;
		G4double primaryColl2Height=51.5;
		G4double primaryColl2PosX=0.0;
		G4double primaryColl2PosY=0.0;
		G4double primaryColl2PosZ=( primaryCollPosZ +  primaryColl2Height );

		G4Tubs* primaryColl2 = new G4Tubs("primaryColl2", primaryColl2InnerDiam/2, primaryColl2OuterDiam/2, primaryColl2Height/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* primaryColl2Logic = new G4LogicalVolume(primaryColl2, /*material=*/Iron, /*name=*/"primaryColl2Logic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* primaryColl2Phys = new G4PVPlacement(zero_rotation, G4ThreeVector(primaryColl2PosX,primaryColl2PosY,primaryColl2PosZ), primaryColl2Logic, "primaryColl2Phys", container, false, 0);
		//primaryColl2Logic.SetForceWireFrame(1)	

		//#sourceAssembly.AddPlacedVolume( primaryColl2Logic, pos=(primaryColl2PosX,primaryColl2PosY,primaryColl2PosZ) )  
		
		//####14mm coll
		
		G4double coll14OuterDiam=20;
		G4double coll14Height=89;
		G4double coll14PosX=0.0;
		G4double coll14PosY=0.0;
		G4double coll14PosZ=( primaryColl2PosZ + coll14Height/2 + primaryColl2Height/2 + 1.5 );
		G4Tubs* coll14body = new G4Tubs("coll14body", 0.0, coll14OuterDiam/2, coll14Height/2, 0.0, 2.0*pi);
//#bad dir
//#		G4Cons* coll14coneontheleft = new G4Cons("coll14coneontheleft",0.0, 5/2, 0.0, 6.15/2, 1.0/2, 0.0, 2.0*pi);
//#		G4Cons* coll14coneinthemiddle = new G4Cons("coll14coneinthemiddle",0.0, 7.85/2, 0.0, 5/2, 87.0/2, 0.0, 2.0*pi);
//#		G4Cons* coll14coneontheright = new G4Cons("coll14coneontheright",0.0, 9.0/2, 0.0, 7.85/2, 1.0/2, 0.0, 2.0*pi);
//#good dir
		G4Cons* coll14coneontheleft = new G4Cons("coll14coneontheleft",0.0, 7.85/2, 0.0, 9.0/2, 1.0/2, 0.0, 2.0*pi);
		G4Cons* coll14coneinthemiddle = new G4Cons("coll14coneinthemiddle",0.0, 5/2, 0.0, 7.85/2, 87.0/2, 0.0, 2.0*pi);
		G4Cons* coll14coneontheright = new G4Cons("coll14coneontheright",0.0, 6.15/2, 0.0, 5/2, 1.0/2, 0.0, 2.0*pi);
		//uncomment later
		G4RotationMatrix* rotationa = new G4RotationMatrix();
		G4SubtractionSolid* coll14a = new G4SubtractionSolid("coll14body-coll14coneontheleft", coll14body, coll14coneontheleft, rotationa, G4ThreeVector(0.0 ,0.0, 89.0/2));
		G4RotationMatrix* rotationb = new G4RotationMatrix();
		G4SubtractionSolid* coll14b = new G4SubtractionSolid("coll14body-coll14coneontheleft-coll14coneinthemiddle", coll14a, coll14coneinthemiddle,rotationa , G4ThreeVector(0.0 ,0.0, 0.0));
		G4RotationMatrix* rotation = new G4RotationMatrix();
		G4SubtractionSolid* coll14 = new G4SubtractionSolid("coll14body-coll14coneontheleft-coll14coneinthemiddle-coll14coneontheright", coll14b, coll14coneontheright,rotation , G4ThreeVector(0.0 ,0.0, -89.0/2));

		//uncomment later
		G4LogicalVolume* coll14Logic = new G4LogicalVolume(coll14, /*material=*/Iron, /*name=*/"coll14Logic");
				//, /*color=*/jawbone);
		//uncomment later
		G4PVPlacement* coll14Phys = new G4PVPlacement(zero_rotation, G4ThreeVector(coll14PosX,coll14PosY,coll14PosZ), coll14Logic, "coll14Phys", container, false, 0);
		//#coll14Logic.SetForceWireFrame(1)	
		//coll14Logic.SetForceSolid(1)

		//#sourceAssembly.AddPlacedVolume( coll4Logic, pos=(coll4PosX,coll4PosY,coll4PosZ) )  

		//#sourceAssembly.MakeImprint( logicWorld, pos=(0, 0, 0) )       

		
		//####target
		//#target=G4Core.G4Box("target", 1, 1, 1)
		//#G4LogicalVolume* targetLogic= new G4LogicalVolume(target, /*material=*/Iron, /*name=*/"targetLogic", /*color=*/(0,1,0));
		//#targetPhys=G4Support.Placement(targetLogic, "targetPhys", logicWorld, pos=(0, 0, coll4PosZ+coll4Height/2+100) )


		//####Source shielding  

		G4double sourceShieldingInnerDiam=27.6352;
		G4double sourceShieldingOuterDiam=107.6352;
		G4double sourceShieldingHeight=50.0126 + 103;
		G4double sourceShieldingPosX=0.0;
		G4double sourceShieldingPosY=0.0;
		G4double sourceShieldingPosZ= sourceShieldingHeight/2 - containerHeight/2;

		G4Tubs* sourceShielding = new G4Tubs("sourceShielding", sourceShieldingInnerDiam/2, sourceShieldingOuterDiam/2, sourceShieldingHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* sourceShieldingLogic = new G4LogicalVolume(sourceShielding, /*material=*/Iron, /*name=*/"sourceShieldingLogic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* sourceShieldingPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(sourceShieldingPosX,sourceShieldingPosY,sourceShieldingPosZ), sourceShieldingLogic, "sourceShieldingPhys", container, false, 0);
		//uncomment later
		//sourceShieldingLogic->SetForceWireFrame(true);
		//#sourceShieldingLogic->SetForceSolid(true);

		//####Collimator shielding  

		G4double collimatorShieldingInnerDiam=20;
		G4double collimatorShieldingOuterDiam=107.6352;
		G4double collimatorShieldingHeight=89;
		G4double collimatorShieldingPosX=0.0;
		G4double collimatorShieldingPosY=0.0;
		G4double collimatorShieldingPosZ= ( primaryColl2PosZ + coll14Height/2 + primaryColl2Height/2 + 1.5 );

		G4Tubs* collimatorShielding = new G4Tubs("collimatorShielding", collimatorShieldingInnerDiam/2, collimatorShieldingOuterDiam/2, collimatorShieldingHeight/2., 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* collimatorShieldingLogic = new G4LogicalVolume(collimatorShielding, /*material=*/Iron, /*name=*/"collimatorLogic");
				//, /*color=*/skin);
		//uncomment later
		G4PVPlacement* collimatorShieldingPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(collimatorShieldingPosX,collimatorShieldingPosY,collimatorShieldingPosZ), collimatorShieldingLogic, "collimatorShieldingPhys", container, false, 0);
		//collimatorShieldingLogic.SetForceWireFrame(1)
		//#collimatorShieldingLogic.SetForceSolid(1)
		
		
		//####plexi
		G4double plexiInnerDiam= 0.0;
		G4double plexiOuterDiam= 107.6352;
		G4double plexiHeight = 2.0;
		G4double plexiPosX=0.0;
		G4double plexiPosY=0.0;
		G4double plexiPosZ=0.0;
		
		G4Tubs* plexiSolid = new G4Tubs("plexiSolid", plexiInnerDiam/2, plexiOuterDiam/2, plexiHeight/2, 0.0, 2.0*pi);
		//uncomment later
		G4LogicalVolume* plexiLogic = new G4LogicalVolume(plexiSolid, /*material=*/Plexiglas, /*name=*/"plexi");
				//, /*color=*/fat);
		//uncomment later
		G4PVPlacement* plexiPhys = new G4PVPlacement(zero_rotation, G4ThreeVector(plexiPosX,plexiPosY,plexiPosZ), plexiLogic, "plexiPhys", logicWorld, false, 0);
		//plexiLogic.SetForceSolid(1)

 /* 
  //                               
  // Calorimeter
  //  
  solidCalor=0; logicCalor=0; physiCalor=0;
  solidLayer=0; logicLayer=0; physiLayer=0;
  solidLayer2=0; logicLayer2=0; physiLayer2=0;
  solidLayer3=0; logicLayer3=0; physiLayer3=0;
  
  if (CalorThickness > 0.)  
    { solidCalor = new G4Box("Calorimeter",              //its name
                           CalorThickness/2,CalorSizeYZ/2,CalorSizeYZ/2);//size
                                 
      logicCalor = new G4LogicalVolume(solidCalor,       //its solid
                                       defaultMaterial,  //its material
                                       "Calorimeter");   //its name
                                           
      physiCalor = new G4PVPlacement(0,                  //no rotation
                                     G4ThreeVector(),    //at (0,0,0)
                                     logicCalor,         //its logical volume
                                     "Calorimeter",      //its name
                                     logicWorld,         //its mother  volume
                                     false,              //no boolean operation
                                     0);                 //copy number
  
  //                                 
  // Layer
  //
      solidLayer = new G4Box("Layer",                    //its name
                       LayerThickness/2,CalorSizeYZ/2,CalorSizeYZ/2); //size
                       
      logicLayer = new G4LogicalVolume(solidLayer,       //its solid
                                       defaultMaterial,  //its material
                                       "Layer");         //its name
      if (NbOfLayers > 1)                                      
        physiLayer = new G4PVReplica("Layer",            //its name
                                      logicLayer,        //its logical volume
                                      logicCalor,        //its mother
                                      kXAxis,            //axis of replication
                                      NbOfLayers,        //number of replica
                                      LayerThickness);   //width of replica
      else
        physiLayer = new G4PVPlacement(0,                //no rotation
                                      G4ThreeVector(-CalorThickness/2 + LayerThickness/2,0.,0.),   //position
                                      logicLayer,        //its logical volume                                     
                                      "Layer",           //its name
                                      logicCalor,        //its mother  volume
                                      false,             //no boolean operation
                                      0);                //copy number     
    }     

  // Layer 2
	solidLayer2 = new G4Box ("Layer2",
		      Layer2Thickness/2,CalorSizeYZ/2,CalorSizeYZ/2);

	logicLayer2 = new G4LogicalVolume (solidLayer2,
		      			   defaultMaterial,
					   "Layer2");

	physiLayer2 = new G4PVPlacement(0,
					G4ThreeVector((-CalorThickness/2 + LayerThickness + Layer2Thickness/2),0,0),
					logicLayer2,

					"Layer2",
					logicCalor,
					false,
					0);
                    
  // Layer 3
	solidLayer3 = new G4Box ("Layer3",
		      Layer3Thickness/2, CalorSizeYZ/2,CalorSizeYZ/2);

	logicLayer3 = new G4LogicalVolume (solidLayer3,
					   defaultMaterial,
					   "Layer3");

	physiLayer3 = new G4PVPlacement(0,
					G4ThreeVector(LayerThickness + Layer2Thickness + Layer3Thickness/2 - CalorThickness/2,0,0),
					logicLayer3,

					"Layer3",
					logicCalor,
					false,
					0);  
  //                               
  // Absorber
  //
  solidAbsorber=0; logicAbsorber=0; physiAbsorber=0;  
  
  if (AbsorberThickness > 0.) 
    { solidAbsorber = new G4Box("Absorber",              //its name
                          AbsorberThickness/2,CalorSizeYZ/2,CalorSizeYZ/2); 
                          
      logicAbsorber = new G4LogicalVolume(solidAbsorber,    //its solid
                                          AbsorberMaterial, //its material
                                          AbsorberMaterial->GetName()); //name
                                                
      physiAbsorber = new G4PVPlacement(0,                 //no rotation
                          G4ThreeVector(-GapThickness/2,0.,0.),  //its position
                                        logicAbsorber,     //its logical volume                    
                                        AbsorberMaterial->GetName(), //its name
                                        logicLayer,        //its mother
                                        false,             //no boolean operat
                                        0);                //copy number
                                        
    }

  //Absorber2
  solidAbsorber2=0;  logicAbsorber2=0;  physiAbsorber2=0;
  if (Absorber2Thickness>0)
	{ solidAbsorber2 = new G4Box("Absorber2",
			      Absorber2Thickness/2,CalorSizeYZ/2,CalorSizeYZ/2);

	  logicAbsorber2 = new G4LogicalVolume(solidAbsorber2,
						Absorber2Material,
						Absorber2Material->GetName());

	  physiAbsorber2 = new G4PVPlacement(0,
				G4ThreeVector(-Gap2Thickness/2,0.,0.),
				logicAbsorber2,

				Absorber2Material->GetName(),
				logicLayer2,
				false,
				0);
  }

  //Absorber3
  solidAbsorber3=0;  logicAbsorber3=0;  physiAbsorber3=0;
  if (Absorber3Thickness>0)
	{ solidAbsorber3 = new G4Box("Absorber3",
			      Absorber3Thickness/2,CalorSizeYZ/2,CalorSizeYZ/2);

	  logicAbsorber3 = new G4LogicalVolume(solidAbsorber3,
						Absorber3Material,
						Absorber3Material->GetName());

	  physiAbsorber3 = new G4PVPlacement(0,
				G4ThreeVector(-Gap3Thickness/2,0.,0.),
				logicAbsorber3,

				Absorber3Material->GetName(),
				logicLayer3,
				false,
				0);
  }

  
  //                                 
  // Gap
  //
  solidGap=0; logicGap=0; physiGap=0; 
  
  if (GapThickness > 0.)
    { solidGap = new G4Box("Gap",
                               GapThickness/2,CalorSizeYZ/2,CalorSizeYZ/2);
                               
      logicGap = new G4LogicalVolume(solidGap,
                                           GapMaterial,
                                           GapMaterial->GetName());
                                           
      physiGap = new G4PVPlacement(0,                      //no rotation
               G4ThreeVector(AbsorberThickness/2,0.,0.),   //its position
                                   logicGap,               //its logical volume               
                                   GapMaterial->GetName(), //its name
                                   logicLayer,             //its mother
                                   false,                  //no boolean operat
                                   0);                     //copy number
    }

  //Gap2
  solidGap2=0;  logicGap2=0;  physiGap2=0;
  if (Gap2Thickness>0)
  {  solidGap2 = new G4Box("Gap2",
			Gap2Thickness/2,CalorSizeYZ/2,CalorSizeYZ/2);
    
     logicGap2 = new G4LogicalVolume(solidGap2,
					Gap2Material,
					Gap2Material->GetName());

     physiGap2 = new G4PVPlacement(0,
			G4ThreeVector((Absorber2Thickness/2),0.,0.),
			logicGap2,

			Gap2Material->GetName(),
			logicLayer2,
			false,
 			0);
  }

  //Gap3
  solidGap3=0;  logicGap3=0; physiGap3=0;
  if (Gap3Thickness>0)
  {  solidGap3 = new G4Box("Gap3",
			Gap3Thickness/2, CalorSizeYZ/2, CalorSizeYZ/2);

     logicGap3 = new G4LogicalVolume(solidGap3,
					Gap3Material,
					Gap3Material->GetName());

     physiGap3 = new G4PVPlacement(0,
			G4ThreeVector(Absorber3Thickness/2,0.,0.),
			logicGap3,

			Gap3Material->GetName(),
			logicLayer3,
			false,
			0);
  }
    
  //PrintCalorParameters();     
  
  //                                        
  // Visualization attributes
  //
  */
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

  //G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  //simpleBoxVisAtt->SetVisibility(true);
  //logicCalor->SetVisAttributes(simpleBoxVisAtt);

 /*
  // Below are vis attributes that permits someone to test / play 
  // with the interactive expansion / contraction geometry system of the
  // vis/OpenInventor driver :
 {G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  simpleBoxVisAtt->SetVisibility(true);
  delete logicCalor->GetVisAttributes();
  logicCalor->SetVisAttributes(simpleBoxVisAtt);}

 {G4VisAttributes* atb= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  logicLayer->SetVisAttributes(atb);}
  
 {G4VisAttributes* atb= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  atb->SetForceSolid(true);
  logicAbsorber->SetVisAttributes(atb);}
  
 {//Set opacity = 0.2 then transparency = 1 - 0.2 = 0.8
  G4VisAttributes* atb= new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.2));
  atb->SetForceSolid(true);
  logicGap->SetVisAttributes(atb);}
  */

  //
  //always return the physical World
  //
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n------------------------------------------------------------"
         << "\n---> The calorimeter is " << NbOfLayers << " layers of: [ "
         << AbsorberThickness/mm << "mm of " << AbsorberMaterial->GetName() 
         << " + "
         << GapThickness/mm << "mm of " << GapMaterial->GetName() << " ] " 
         << "\n------------------------------------------------------------\n";
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) AbsorberMaterial = pttoMaterial;
}

void DetectorConstruction::SetAbsorber2Material(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) Absorber2Material = pttoMaterial;
}

void DetectorConstruction::SetAbsorber3Material(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) Absorber3Material = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) GapMaterial = pttoMaterial;
}

void DetectorConstruction::SetGap2Material(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) Gap2Material = pttoMaterial;
}

void DetectorConstruction::SetGap3Material(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) Gap3Material = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  AbsorberThickness = val;
}

void DetectorConstruction::SetAbsorber2Thickness(G4double val)
{
  Absorber2Thickness = val;
}

void DetectorConstruction::SetAbsorber3Thickness(G4double val)
{
  Absorber3Thickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapThickness(G4double val)
{
  // change Gap thickness and recompute the calorimeter parameters
  GapThickness = val;
}

void DetectorConstruction::SetGap2Thickness(G4double val)
{
  Gap2Thickness = val;
}

void DetectorConstruction::SetGap3Thickness(G4double val)
{
  Gap3Thickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCalorSizeYZ(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  CalorSizeYZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfLayers(G4int val)
{
  NbOfLayers = val;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(magField) delete magField;                //delete the existing magn field

  if(fieldValue!=0.)                        // create a new one if non nul
  { magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } else {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
