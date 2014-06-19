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

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleSource  = new G4GeneralParticleSource();
  Detector = (DetectorConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction();  
  
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);
  rndmFlag = "on";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleSource;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 
  G4double x0 = -0.5*(Detector->GetWorldSizeX());
  G4double y0 = 0.*cm, z0 = 0.*cm;
  if (rndmFlag == "on")
     {y0 = (Detector->GetCalorSizeYZ())*(G4UniformRand()-0.5);
      z0 = (Detector->GetCalorSizeYZ())*(G4UniformRand()-0.5);
     } 
  particleSource->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  particleSource->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
#source coordinates
sources= [[0]*4 for _ in [0]*30]

  	
sources[0][0] = 14.3	
sources[0][1] = 52.5	
sources[0][2] = 407.4	
sources[0][3] = 103.8	
  	
sources[1][0] = 15.6	
sources[1][1] = 43.75	
sources[1][2] = 405.5	
sources[1][3] = 113.2	
  	
sources[2][0] = 16.9	
sources[2][1] = 35.	
sources[2][2] = 403.5	
sources[2][3] = 122.6	
  	
sources[3][0] = 18.2	
sources[3][1] = 26.25	
sources[3][2] = 401.3	
sources[3][3] = 131.9	
  	
sources[4][0] = 19.5	
sources[4][1] = 17.5	
sources[4][2] = 398.8	
sources[4][3] = 141.2	
  	
sources[5][0] = 20.8	
sources[5][1] = 8.75	
sources[5][2] = 396.1	
sources[5][3] = 150.5	
 	
sources[6][0] = 22.1	
sources[6][1] = 0.	
sources[6][2] = 393.2	
sources[6][3] = 159.7	
  	
sources[7][0] = 23.4	
sources[7][1] = 52.5	
sources[7][2] = 390.	
sources[7][3] = 168.8	

sources[8][0] = 24.7	
sources[8][1] = 42.	
sources[8][2] = 386.7	
sources[8][3] = 177.9	
  	
sources[9][0] = 26.	
sources[9][1] = 31.5	
sources[9][2] = 383.1	
sources[9][3] = 186.9	
  	
sources[10][0] = 27.3	
sources[10][1] = 21.	
sources[10][2] = 379.4	
sources[10][3] = 195.8	


sources[11][0] = 28.6	
sources[11][1] = 10.5	
sources[11][2] = 375.4	
sources[11][3] = 204.7	
  	  	
sources[12][0] = 29.9	
sources[12][1] = 0.	
sources[12][2] = 371.2	
sources[12][3] = 213.4
  	
sources[13][0] = 31.2	
sources[13][1] = 52.5	
sources[13][2] = 366.7	
sources[13][3] = 222.1	
  	
sources[14][0] = 32.5	
sources[14][1] = 42.	
sources[14][2] = 362.1	
sources[14][3] = 230.7	
  	
sources[15][0] = 33.8	
sources[15][1] = 31.5	
sources[15][2] = 357.4	
sources[15][3] = 239.2	
 	
sources[16][0] = 35.1	
sources[16][1] = 21.	
sources[16][2] = 352.2	
sources[16][3] = 247.6	
  	
sources[17][0] = 36.4	
sources[17][1] = 10.5	
sources[17][2] = 347.	
sources[17][3] = 255.8	
  	
sources[18][0] = 37.7	
sources[18][1] = 0.	
sources[18][2] = 341.5	
sources[18][3] = 264.	
  	  	
sources[19][0] = 39.	
sources[19][1] = 52.5	
sources[19][2] = 335.9	
sources[19][3] = 272.	
  	  	
sources[20][0] = 40.3	
sources[20][1] = 42.	
sources[20][2] = 330.1	
sources[20][3] = 279.9	
 	  	
sources[21][0] = 41.6	
sources[21][1] = 31.5	
sources[21][2] = 324.	
sources[21][3] = 287.7	
  	  	
sources[22][0] = 42.9	
sources[22][1] = 21.	
sources[22][2] = 317.8	
sources[22][3] = 295.3	
  	  	
sources[23][0] = 44.2	
sources[23][1] = 10.5	
sources[23][2] = 311.4	
sources[23][3] = 302.8	
  	  	
sources[24][0] = 45.5	
sources[24][1] = 0.	
sources[24][2] = 304.8	
sources[24][3] = 310.2	
  	  
sources[25][0] = 46.8	
sources[25][1] = 52.	
sources[25][2] = 298.	
sources[25][3] = 317.4	
  	  	
sources[26][0] = 48.1	
sources[26][1] = 40.	
sources[26][2] = 291.1	
sources[26][3] = 324.4	
  	  	
sources[27][0] = 49.4	
sources[27][1] = 28.	
sources[27][2] = 283.9	
sources[27][3] = 331.3	
	  	  	
sources[28][0] = 50.7	
sources[28][1] = 16.	
sources[28][2] = 276.6	
sources[28][3] = 338.	
  	  	
sources[29][0] = 52.	# T
sources[29][1] = 4.	    # Q
sources[29][2] = 269.2	# R
sources[29][3] = 344.5	# Z+30
#source coordinates end


import sys
pV=getFlag('pV', 0, int) #paralell value
nIons=getFlag('nIons', 0, int)
sDist=getFlag('sDist', 10, int)

outFile=getFlag('outFile', 
	default='notused.'+time.strftime("%Y%m%d.%H%M%S")+'.aida',
	datatype=str)
gunEnergy=getFlag('gunEnergy', 180.0, float)

sourceDistance=sDist	 # sDist means the zubal distance from secondary collimator in cm
length=80+10*sourceDistance
coll14Height=89
containerHeight = 50.0126 + 51.5*2 + 89 + 1.5
plexiHeight=10.0

print "********* Leftover arguments: ", sys.argv

mainVars=setupGeant(detector=construction, detector_args=dict(orbOnly=False)
)

G4Sys=mainVars.G4Sys

print G4Core.FormatMaterialTable()
	
mainVars.myTrackingAction.SetShowPrimaries(1)
mainVars.myTrackingAction.SetShowGammas(1)

RadiationPhysics.ShowParticles()

ion=G4Core.FindParticle('gamma')
mainVars.gun.SetParticleDefinition(ion)
mainVars.G4Sys.SetWorldCuts(1*millimeter)

detinfo=mainVars.detectorInfo
try:
	compressor=detinfo.compressor
	compressor.SetZVisSampling( (detinfo.lastlayer+detinfo.firstlayer)//2,1000)	
	compressor.SetXVisSampling(detinfo.xcount//2,10000)	
	compressor.SetYVisSampling(detinfo.ycount//2,1000)	
	for i in range(256): 
		compressor.SetOrganSensitiveDetector(i, mainVars.detectorInfo.PhantomDet)
except AttributeError:
	pass #compressor is probably an orb!

def set_source(new_point, rotationX, rotationZ, thetaDeg, gunAxis):
	"""set up the source and collimators for a given angle and hypocenter"""
	theta=thetaDeg*3.14159265/180.0
	new_point=Vector(new_point)

	gen=mainVars.myGenerator
	mainVars.myGenerator.SetBeamSize(2.5908/2)    #smallestTubeInnerDiam=2.5908
	mainVars.myGenerator.SetAngularHalfSpread(math.pi) #distribute beam over whole target

#	coords for Gun position
	containerHeight = 50.0126 + 51.5*2 + 89 + 1.5
	largestTubeHeight = 50.0126
	largestTubePosZ=largestTubeHeight/2 - containerHeight/2
	smallestTubeHeight = 32.9438 * mm
	smallestBackCapHeight = 2 * mm
	smallestCapHeight = 0.685 * mm
	cobaltHeight = smallestTubeHeight - smallestBackCapHeight - smallestCapHeight
	largestTubeHeight = 50.0126 * mm
	smallestTubePosZ = ((largestTubeHeight/mm-smallestTubeHeight/mm)/2-1.4224 -0.4064 )*mm
	cobaltPosZ = (smallestTubePosZ/mm+(smallestTubeHeight/mm-cobaltHeight/mm)/2-smallestCapHeight/mm)*mm 
	insideFillingX = 0.0 
	insideFillingY = 0.0 
	fillingPosZ= cobaltPosZ+largestTubeHeight/2 - containerHeight/2
	print new_point,rotationX,rotationZ,theta
	#lenti kettot meg kell rotalni thetaval, kulonben a forrasok a helyukon maradnak, csak a vas megy
	rotgunAxis=gunAxis.rotate(theta,0,0)
	rotnew_point=new_point.rotate(theta,0,0)
	gen.SetGunAxis(rotgunAxis)
	# veget transzfotmaljuk a source[][]-szal
	gen.SetGunPosition(rotnew_point+containerHeight/2*rotgunAxis+fillingPosZ*rotgunAxis)  #
	mainVars.detectorInfo.containerPhys.relocate( pos=new_point+containerHeight/2*rotgunAxis, rot=rotationX*rotationZ*Matrix(RotationZ(theta)))
	mainVars.detectorInfo.plexiPhys.relocate( pos=new_point+(containerHeight+5)*rotgunAxis, rot=rotationX*rotationZ*Matrix(RotationZ(theta)))  # ez nem biztos, hogy a helyen van

	#kozepe trafo
	#gen.SetGunPosition(rotnew_point+fillingPosZ*rotgunAxis)  #
	#mainVars.detectorInfo.containerPhys.relocate( pos=new_point, rot=rotationX*rotationZ*Matrix(RotationZ(theta)))
	#mainVars.detectorInfo.plexiPhys.relocate( pos=new_point+(28.5+89+5)*gunAxis, rot=rotationX*rotationZ*Matrix(RotationZ(theta)))  # ez nem biztos, hogy a helyen van



if nIons>0:		
	def createTree(name=None, id=0, memory_only=True):
		import AIDASupport #only import AIDA if we really need it!
		if name is None:
			runcode=time.strftime("%Y%m%d.%H%M%S")
			name="4pi-run-beam30_%d-%.0e.%s"%(pV,id,runcode)
			if memory_only:
				dataTree=AIDASupport.AIDATree('')
			else:
				dataTree=AIDASupport.AIDATree(name+'.hdf5', createNew=1, storeType='hdf5', options="")
		return name, dataTree
		

	###
	if 0:
		SetupViewer(name="OGLIX 800", delayed=True)
		view_angle(90,180)
		view_accumulate()
		view_on()
		view_zoom(0.5)
	
	beamZ = 0

	mainVars.myEvent.SetPrintInterval(10000)
	totalIons=0
	#mainVars.G4Sys.SetupRun()
	rm=mainVars.G4Sys.runMgr

	runName, dataTree=createTree(id=nIons)
	if 0:
		hist=dataTree.GetHistogramFactory().createHistogram3D("energy","energy", 150, -85., 85.,
									  150, -85., 85., 150, -85, 85*mm)
		dose_map=smoother_debrecen.AIDA_dose_map(hist)
	else:
		import numpy
		hist_data=numpy.zeros((150,150,150), numpy.float32)
		hist=G4Core.Numpy3DCalorimeter("energy", hist_data, -85, 85, -85, 85, -85, 85)
		dose_map=smoother_debrecen.numpy_dose_map(hist)

	mainVars.detectorInfo.PhantomDet.SetXYZHistogram(hist)

	gen=mainVars.myGenerator
	
	mat_nx, mat_ny, mat_nz=detinfo.dataset_shape
	#read in a tissue map which has indices which correspond to mu values in the muvals array

	#materials in test data set are 0=air, 2=tissue, 4=bone
	#these correspond to elements 0, 4, and 9 of the zuba materials data set from 
	#_CVSVers="Id: irradiate_zubal_compressed.py,v 1.33 2009/12/20 20:13:45 marcus Exp"
	zubal_matatten=[
	0.001, 1.00278697894*1.1, 1.26951357429, 1.22682864055, 1.2615128051, 0.3,
		1.05145316682*1.2, 1.0,  1.26392344202, 1.75568066692*1.1*1.05
	]
	muvals_map=[zubal_matatten[0], 0, zubal_matatten[4], 0, zubal_matatten[9]]
		
	muvals=smoother_debrecen.muvals_manager(dose_map, 
		detinfo.data, muvals_map, cache_file_name=detinfo.data.filename+'_cache.pickle')

	#view_zoom(0.2)
	#view_on()
	histinfo=dose_map
	
	smooth=None
	import threading
	
	sigma_long=10*mm
	sigma_transverse=1.4*mm
	
	theta_start=-180
	theta_finish=180
	deltatheta=360
		
	#this function will get run in a thread.  It, in turn, spins off the actual G4ant4 sim in another thread
	#and releases control, so smoothing can be done in parallel
	def runsweep(n_events, thetaDeg, deltatheta):
#		for conv_point in conv_point_list:
		for i in range(30):

		### translation and rotation
			x = math.sin( sources[i][1] *deg) * ( - (sources[i][2]-1.5 ) )
			y = math.sqrt( ( sources[i][2]-1.5 )*(sources[i][2]-1.5 ) - x*x )
			z = -sources[i][3]
	
			rotationX=G4Core.HepRotation()
			rotationZ=G4Core.HepRotation()
			rotationX.rotateX( (-90+sources[i][0])*deg )
			rotationZ.rotateZ( -sources[i][1] *deg )
			source_distance=math.sqrt(x*x+y*y+z*z)
			new_point=Vector((x,y,z))
			#print "HAHO",x,y,z, x*x+y*y+z*z
			gunAxis=Vector((-x/source_distance, -y/source_distance, -z/source_distance ))
		#### rotation end
			for thetaStep in range(thetaDeg, thetaDeg+deltatheta, 360):
				set_source(new_point, rotationX, rotationZ, thetaStep, gunAxis)				
				#randomize((189554476, 1193620722))
				#print "*****seeds: ", thetaStep, (mainVars.myGenerator.GetSeedWord(0), mainVars.myGenerator.GetSeedWord(1))
				th=G4Support.RunAsync(n_events) #this is qthreaded!
				th.join() #let it finish
				#print "*****final seeds: ", thetaStep,  get_seeds(mainVars)
			
	for thetaDeg in range(theta_start, theta_finish+1, deltatheta): #take final step at end, by setting the upper bound to end+1, to get last smoothing done
		#create a thread with the runsweep function in it, and let it run
		if thetaDeg != theta_finish: #not final pass, go ahead and start a sim
			histinfo.reset()
			th=threading.Thread(target=runsweep, args=(nIons, thetaDeg, deltatheta))
			th.start()
				
		if thetaDeg != theta_start: #only do this after first pass, when data will be available
			kern=smoother_debrecen.construct_kernel(histinfo, (0, math.sin(sim_center_theta), math.cos(sim_center_theta) ) , sigma_long, sigma_transverse )
			if smooth is None:
				smooth=smoother_debrecen.do_smooth(kern, histinfo.dose, muvals.muvals, destructive=False)
			else:
				smooth +=smoother_debrecen.do_smooth(kern, histinfo.dose, muvals.muvals, destructive=False)

		if thetaDeg != theta_finish: 
			th.join() #make sure runsweep is finished by joining its thread BU-join: wait until thread is terminated
			histinfo.reload() #transfer data from AIDA hist to usable form before we start next run and overwrite
			sim_center_theta=(thetaDeg+0.5*deltatheta)*3.14159265/180.0 #center of this angular range, for smoothing kernel, saved for next pass!
		

	resamp=muvals.convert_flux_to_dose(smooth)
	smoother_debrecen.create_dx_output(runName, muvals, resamp)
*/
