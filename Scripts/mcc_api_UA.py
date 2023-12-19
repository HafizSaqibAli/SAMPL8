"""
This script uses the CodeEntropy API to compute MCC 
entropy of a molecule at its united atom (UA) level.
CodeEntropy package should be installed for this script to work.
Download the package from: https://github.com/arghya90/CodeEntropy.git

The instructions to install it are in the README.md file.

NOTE: READ THE CODE CAREFULLY TO UNDERSTAND WHAT IS HAPPENING.

"""

import os, sys, re
import numpy as nmp
from CodeEntropy.Reader import GromacsReader
from CodeEntropy.Reader import Constants as CONST
from CodeEntropy.IO import Writer
from CodeEntropy.FunctionCollection import Utils, EntropyFunctions, GeometricFunctions, CustomFunctions
from CodeEntropy.FunctionCollection import UnitsAndConversions as UAC
from CodeEntropy.ClassCollection import DataContainer, BeadClasses, ModeClasses, CustomDataTypes, ConformationEntity


## Function definitions
def get_svib_UA(arg_hostDataContainer, \
				arg_baseMolecule, \
				arg_outFile, \
				arg_fScale, \
				arg_tScale, \
				arg_verbose, \
				arg_moutFile, \
				arg_nmdFile, \
				arg_temper):
	Utils.printflush('-'*60)
	Utils.printflush("{:^60}".format("Hierarchy level. --> United Atom <--"))
	Utils.printflush('-'*60)

	Utils.printOut(arg_outFile,'-'*60)
	Utils.printOut(arg_outFile,"{:^60}".format("Hierarchy level. --> United Atom <--"))
	Utils.printOut(arg_outFile,'-'*60)

	# preparing header for output file
	Utils.printOut(arg_outFile,'	  {:<10s}{:>5s}{:>12s}{:>12s}'.format('RESNAME', 'RESID', 'FF_ENTROPY', 'TT_ENTROPY'))

	# initialize total entropy values
	totalUAEntropyFF = 0.
	totalUAEntropyTT = 0.

	# number of frames
	numFrames = len(arg_hostDataContainer.trajSnapshots)

	#reset
	arg_hostDataContainer.reset_rotationAxesArray()
	arg_hostDataContainer.reset_translationAxesArray()

	# for each residue:
	for rid in range(arg_baseMolecule.numResidues):
		resLabel = "{}{}".format(arg_baseMolecule.resnameArray[rid], arg_baseMolecule.residArray[rid])
		Utils.printflush('Working on resid : {}'.format(resLabel))

		# create a bead collection 
		ridBeadCollection = BeadClasses.BeadCollection("{}_bead".format(resLabel),arg_hostDataContainer)
		ridBeadCollection.listOfBeads = []

		# add UA beads to it
		atoms_in_rid = []
		iAtom_in_rid = int(arg_baseMolecule.residueHead[rid])
		while iAtom_in_rid != -1:

			atoms_in_rid.append(iAtom_in_rid)

			# if it is not a hydrogen
			if not arg_baseMolecule.isHydrogenArray[iAtom_in_rid]:
				atomList = [iAtom_in_rid]

				numBondedHydrogenAtoms = len(arg_baseMolecule.bondedHydrogenTable[iAtom_in_rid]) 
				if numBondedHydrogenAtoms != 0:
					# add all the bonded hydrogens to its list
					for iPriorityH,iH in arg_baseMolecule.bondedHydrogenTable[iAtom_in_rid].list_in_order():
						if iH != -1:
							atomList.append(iH)

				# create a bead
				newBead = BeadClasses.Bead(arg_atomList=atomList,
								  arg_hostDataContainer=arg_hostDataContainer,
								  arg_numFrames=numFrames,
									arg_beadName = arg_baseMolecule.atomNameArray[iAtom_in_rid],
									arg_beadResi = arg_baseMolecule.residArray[rid],
									arg_beadResn = arg_baseMolecule.resnameArray[rid],
									arg_beadChid = "X")
				newBead.position = arg_hostDataContainer._labCoords[0, iAtom_in_rid]
				ridBeadCollection.listOfBeads.append(newBead)
			else:
				pass

			iAtom_in_rid = int(arg_baseMolecule.atomArray[iAtom_in_rid])
		#END WHILE

		# by this point, the UA beads for that residue have been created
		Utils.printflush('Total number of UA beads in residue {} : {}'\
					  .format(resLabel, len(ridBeadCollection.listOfBeads)))

		# reset weighted vectors for each bead
		for iBead in ridBeadCollection.listOfBeads:
			iBead.reset_totalWeightedVectors( (numFrames,3) )

		# reseting all the F-T combo matrices to zero
		ridBeadCollection.reinitialize_matrices()

		# setup Translation and Rotation axes
		# Translation axes : each atom is in the principal axes of its host residue
		Utils.printflush("Assigning Translation Axes at the UA level->", end = ' ')
		for iFrame in range(numFrames):
			selMOI, selAxes = arg_baseMolecule\
					  .get_principal_axes(arg_atomList = atoms_in_rid, arg_frame = iFrame, arg_sorted=False)
			selCOM = arg_baseMolecule\
					 .get_center_of_mass(arg_atomList = atoms_in_rid, arg_frame = iFrame)
			arg_hostDataContainer.update_translationAxesArray_at(arg_frame = iFrame, arg_atomList = atoms_in_rid, arg_pAxes = selAxes, arg_orig = selCOM)

		Utils.printflush('Done')

		Utils.printflush("Assigning Rotational Axes at the UA level->", end = ' ')
		for iAtom_in_rid in atoms_in_rid:
			for iFrame in range(numFrames):
				# Rotation axes : 
				# the axes will have the geometry of a 
				# local sphereical-polar coordinate system. 
				if arg_baseMolecule.isHydrogenArray[iAtom_in_rid]:
					# do nothing except do a check 
					# on number of heavy atoms it is bonded to 
					# should be exactly 1
					assert(len(arg_baseMolecule.bondedHeavyAtomTable[iAtom_in_rid]) == 1)

				else:
					# it is a heavy atom
					# from each of the hydrogen atoms bonded to it 
					# get the average position lab coordinate
					avgHydrogenPosition = EntropyFunctions.get_avg_hpos(arg_atom= iAtom_in_rid, \
						arg_frame = iFrame, \
						arg_baseMolecule = arg_baseMolecule, \
						arg_hostDataContainer = arg_hostDataContainer)

					# use the resultant vector to generate an 
					# orthogonal local coordinate axes system
					# with origin at the heavy atom position
					heavyOrigin = arg_hostDataContainer._labCoords[iFrame, iAtom_in_rid]
					iAtomBasis = nmp.vstack((GeometricFunctions.get_sphCoord_axes(arg_r=avgHydrogenPosition), heavyOrigin))

					# heavy atom
					arg_hostDataContainer.localCoords[iFrame, iAtom_in_rid] = arg_hostDataContainer._labCoords[iFrame, iAtom_in_rid] - iAtomBasis[-1]
					arg_hostDataContainer.localCoords[iFrame, iAtom_in_rid] = iAtomBasis[0:3,] @ arg_hostDataContainer.localCoords[iFrame, iAtom_in_rid]
					arg_hostDataContainer.rotationAxesArray[iFrame, iAtom_in_rid] = iAtomBasis

					# its hydrogens
					for iPriorityH, iH in arg_baseMolecule.bondedHydrogenTable[iAtom_in_rid].list_in_order():
						arg_hostDataContainer.localCoords[iFrame, iH] = arg_hostDataContainer._labCoords[iFrame, iH] - iAtomBasis[-1]
						arg_hostDataContainer.localCoords[iFrame, iH] = iAtomBasis[0:3,] @ arg_hostDataContainer.localCoords[iFrame, iH]
						arg_hostDataContainer.rotationAxesArray[iFrame, iH] = iAtomBasis
		Utils.printflush('Done')

		# cast lab forces onto these axes simulateneously 
		Utils.printflush('Casting forces at the UA level->',end=' ')
		for iAtom_in_rid in atoms_in_rid:
			for iFrame in range(numFrames):
				arg_hostDataContainer.localForces[iFrame][iAtom_in_rid] = arg_hostDataContainer.translationAxesArray[iFrame,iAtom_in_rid, 0:3] \
																		  @ arg_hostDataContainer._labForces[iFrame][iAtom_in_rid]
		Utils.printflush('Done')


		# update torques in the arg_hostDataContainer.
		# This is redundant because local forces are not the ones that 
		# will be used to get the torques. Local coordinates however are.
		Utils.printflush('Casting torques at the UA level->', end = ' ')
		for iAtom_in_rid in atoms_in_rid:
			for iFrame in range(numFrames):
				coords_i = arg_hostDataContainer.localCoords[iFrame, iAtom_in_rid]
				forces_i = arg_hostDataContainer.rotationAxesArray[iFrame, iAtom_in_rid][0:3,]@arg_hostDataContainer._labForces[iFrame,iAtom_in_rid]
				arg_hostDataContainer.localTorques[iFrame,iAtom_in_rid,:] = CustomFunctions.cross_product(coords_i,forces_i)
		Utils.printflush('Done')


		# mass weighting the forces and torque
		Utils.printflush('Computing total weighted forces and torques on each bead->', end = ' ')
		for iBead in ridBeadCollection.listOfBeads:
			for iFrame in range(numFrames):

				# mass weighting the forces for each bead (iBead) in each direction (j) 
				# inertia weighting the torques for each bead (iBead) in each direction (j)

				# define local basis as the rotationalAxes of the first atom in the atomList of iBead 
				# doesnt matter because they all have the same R and T axes
				iLocalBasis = arg_hostDataContainer.rotationAxesArray[iFrame][iBead.atomList[0]]

				#get the moment of inertia tensor for ibead in thid local basis
				beadMOITensor = iBead.get_moment_of_inertia_tensor_local(arg_localBasis = iLocalBasis, arg_frame = iFrame)

				# get total torque and force in each direction and weight them by √beadMOITensor[jj]
				for j in range(3):
					iBead.totalWeightedForces[iFrame,j] = nmp.sum(nmp.asarray([arg_hostDataContainer.localForces[iFrame,aid,j] for aid in iBead.atomList]))
					iBead.totalWeightedTorques[iFrame,j] = nmp.sum(nmp.asarray([arg_hostDataContainer.localTorques[iFrame,aid,j] for aid in iBead.atomList]))

					try:
						assert(iBead.get_total_mass() != 0)
					except:
						raise AssertionError("Bead with zero mass found. Clearly an error!")

					# mass weight the total force component
					iBead.totalWeightedForces[iFrame,j] /= nmp.sqrt(iBead.get_total_mass())

					try:
						if nmp.isclose(iBead.totalWeightedTorques[iFrame,j] , 0.0):
							# then the beadMOITensor[j,j] must be 0 as well
							# ensure that
							assert(nmp.isclose(beadMOITensor[j,j] , 0.0))
						else:
							# inertia weight the total torque component
							iBead.totalWeightedTorques[iFrame,j] /= nmp.sqrt(beadMOITensor[j,j])
					except:
						raise AssertionError("Moment of Intertia is non-zero for a bead lying on axis {}".format(j))

		Utils.printflush('Done')

		# now fill in the matrices
		Utils.printflush("Updating the submatrices ... ")
		ridBeadCollection.update_subMatrix(arg_pairString="FF",arg_verbose=arg_verbose)
		ridBeadCollection.update_subMatrix(arg_pairString="TT",arg_verbose=arg_verbose)
		Utils.printflush('Done')

		#make quadrant from subMatrices
		Utils.printflush("Generating Quadrants->",end = ' ')
		ffQuadrant = ridBeadCollection.generate_quadrant(arg_pairString="FF",arg_filterZeros=0)
		ttQuadrant = ridBeadCollection.generate_quadrant(arg_pairString="TT",arg_filterZeros=0)
		Utils.printflush("Done")

		# scale forces/torques of these quadrants
		ffQuadrant = nmp.multiply(arg_fScale**2, ffQuadrant)
		ttQuadrant = nmp.multiply(arg_tScale**2, ttQuadrant)

		# remove any row or column with zero axis
		# this could have been done while generating quadrants. Can be merged if wished for
		ffQuadrant = ridBeadCollection.filter_zero_rows_columns(ffQuadrant)
		ttQuadrant = ridBeadCollection.filter_zero_rows_columns(ttQuadrant)


		# print matrices if asked
		if arg_moutFile:
			Writer.write_a_matrix(arg_matrix = ffQuadrant\
								  , arg_descriptor = "FF COV AT UNITED ATOM LEVEL FOR RES {}".format(resLabel)\
								  , arg_outFile = arg_moutFile)
			Writer.write_a_matrix(arg_matrix = ttQuadrant\
								  , arg_descriptor = "TT COV AT UNITED ATOM LEVEL FOR RES {}".format(resLabel)\
								  , arg_outFile = arg_moutFile)

		#diagnolaize
		Utils.printflush("Diagonalizing->", end = ' ')
		lambdasFF, eigVectorsFF  = Utils.diagonalize(ffQuadrant)
		lambdasTT, eigVectorsTT  = Utils.diagonalize(ttQuadrant)
		Utils.printflush('Done')

		# since eigen values can be complex numbers 
		# but with imag parts very close to zero
		# use numpy's real_if_close with some tolerance to mask the imag parts
		# Utils.printflush('Checking the nature of eigen values and conditioning them ...', end = ' ')
		# tol = 1e+5
		# lambdasFF = nmp.real_if_close(lambdasFF/1e+5, tol= tol)
		# lambdasTT = nmp.real_if_close(lambdasTT/1e+5, tol= tol)
		# Utils.printflush('Done')

		# filter real zero values
		lambdasFF = nmp.asarray([lm for lm in lambdasFF if not nmp.isclose(lm, 0.0)])
		lambdasTT = nmp.asarray([lm for lm in lambdasTT if not nmp.isclose(lm, 0.0)])


		# change to SI units
		Utils.printflush('Changing the units of eigen values to SI units->', end = ' ')
		lambdasFF = UAC.change_lambda_units(lambdasFF)
		lambdasTT = UAC.change_lambda_units(lambdasTT)
		Utils.printflush('Done')

		# Create a spectrum to store these modes for 
		# proper output and analyses.
		modeSpectraFF = []
		for midx, mcombo in enumerate(zip(lambdasFF, eigVectorsFF)):
			fflmb, evec = mcombo
			# compute mode frequencies
			# nu = sqrt(lambda/kT)*(1/2pi)
			# Units: 1/s
			mfreq = EntropyFunctions.compute_frequency_from_lambda(fflmb, arg_temper)
			newMode = ModeClasses.Mode(arg_modeIdx = midx + 1, \
				arg_modeEval = fflmb, \
				arg_modeEvec = evec, \
				arg_modeFreq = mfreq)
			newMode.modeAmpl = EntropyFunctions.compute_ampfac_from_lambda(fflmb, arg_temper)
			modeSpectraFF.append(newMode)

		ridBeadCollection.assign_attribute("modeSpectraFF", modeSpectraFF)


		modeSpectraTT = []
		for midx, mcombo in enumerate(zip(lambdasTT, eigVectorsTT)):
			ttlmb, evec = mcombo
			# compute mode frequencies
			# nu = sqrt(lambda/kT)*(1/2pi)
			# Units: 1/s
			mfreq = EntropyFunctions.compute_frequency_from_lambda(ttlmb, arg_temper)
			newMode = ModeClasses.Mode(arg_modeIdx = midx + 1, \
				arg_modeEval = ttlmb, \
				arg_modeEvec = evec, \
				arg_modeFreq = mfreq)
			newMode.modeAmpl = EntropyFunctions.compute_ampfac_from_lambda(ttlmb, arg_temper)
			modeSpectraTT.append(newMode)

		ridBeadCollection.assign_attribute("modeSpectraTT", modeSpectraTT)

		# sorting the spectrum
		Utils.printflush('Sorting spectrum in ascending order of frequencies->', end = ' ')
		ridBeadCollection.modeSpectraFF = ModeClasses.sort_modes(ridBeadCollection.modeSpectraFF)
		ridBeadCollection.modeSpectraTT = ModeClasses.sort_modes(ridBeadCollection.modeSpectraTT)
		Utils.printflush('Done')

		# Print modes if asked
		if arg_nmdFile:
			Writer.append_file(arg_nmdFile)
			ridBeadCollection.write_nmd_file(arg_nmdfile = arg_nmdFile, \
										   arg_spectrum = ridBeadCollection.modeSpectraFF, \
										   arg_wfac = [iBead.get_total_mass() for iBead in ridBeadCollection.listOfBeads])

		# compute entropy
		# 1. remove the smallest 6 freqs from FF sprectrum 
		#	 because they may be overlapping with residue level motions
		# 2. DO NOT remove any freq from TT spectrum because 
		#	they are uncoupled to any TT freq in any other hierarchy
		entropyFF = [EntropyFunctions.calculate_entropy_per_dof(m.modeFreq, arg_temper) for m in ridBeadCollection.modeSpectraFF[6:]]
		entropyTT = [EntropyFunctions.calculate_entropy_per_dof(m.modeFreq, arg_temper) for m in ridBeadCollection.modeSpectraTT[0:]]

		ridTotalEntropyFF = nmp.sum(entropyFF)
		ridTotalEntropyTT = nmp.sum(entropyTT)

		# print final outputs
		Utils.printflush("Entropy values:")

		Utils.printflush('{:<40s} : {:.4f} J/mol/K'.format('FF Entropy (UA for {})'.format(resLabel), ridTotalEntropyFF))
		Utils.printflush('{:<40s} : {:.4f} J/mol/K'.format('TT Entropy (UA for {})'.format(resLabel), ridTotalEntropyTT))
		Utils.printOut(arg_outFile,'UATOM {:<10}{:>5}{:>12.3f}{:>12.3f}'\
								.format(arg_baseMolecule.resnameArray[rid]\
								, arg_baseMolecule.residArray[rid]\
								, ridTotalEntropyFF\
								, ridTotalEntropyTT))

		totalUAEntropyFF += ridTotalEntropyFF
		totalUAEntropyTT += ridTotalEntropyTT

	# Final information	
	Utils.printflush('_'*60)
	Utils.printflush('{:<25} : {:>15.3f} J/mol/K'.format('Total Entropy FF (UA level)', totalUAEntropyFF))
	Utils.printflush('{:<25} : {:>15.3f} J/mol/K'.format('Total Entropy TT (UA level)', totalUAEntropyTT))
	Utils.printflush('-'*60)

	Utils.printOut(arg_outFile,'_'*60)
	Utils.printOut(arg_outFile,'{:<25} : {:>15.3f} J/mol/K'.format('Total Entropy FF (UA level)', totalUAEntropyFF))
	Utils.printOut(arg_outFile,'{:<25} : {:>15.3f} J/mol/K'.format('Total Entropy TT (UA level)', totalUAEntropyTT))
	Utils.printOut(arg_outFile,'-'*60)
	
	return
#END

def get_stopo_UA(arg_baseMolecule, arg_hostDataContainer, arg_outFile, arg_verbose):
	"""
	Compute entropy by Adaptive Enumeration Method (AEM).
	This method deals with each dihedral in a conformational entity on an individual basis. After that it coalesces
	the state vectors of each dihedral in the conformational entity to help compute entropy using p-logP formulation. 
	This function computes the total entropy from all residue in the base molecule.
	"""

	Utils.printflush('-'*60)
	Utils.printflush("{:^60}".format("Topographical entropy of residue side chains \ncomputed using all the dihedrals with AEM method"))
	Utils.printflush('-'*60)

	Utils.printOut(arg_outFile,'-'*60)
	Utils.printOut(arg_outFile,"{:^60}".format("Topographical entropy of residue side chains \ncomputed using all the dihedrals with AEM method"))
	Utils.printOut(arg_outFile,'-'*60)

	# number of frames
	numFrames = len(arg_hostDataContainer.trajSnapshots)

	# log of number of frames (a constant)
	logNumFrames = nmp.log(numFrames)

	# total SC entropy
	totalTopogEntropySC = 0.
	

	# browse through each residue in the system and get their dihedrals
	for rid in range(arg_baseMolecule.numResidues):
		Utils.printflush('-'*10,end='')
		Utils.printflush('Working on resid : {} ({})'.format(rid, arg_baseMolecule.resnameArray[rid]), end='')
		Utils.printflush('-'*10)


		# build a binary tree that will hold unique dihedrals 
		# uniqueness is defined based on 2-3 atom indexes
		diheds_in_rid = CustomDataTypes.BinaryTree()
		iAtom_in_rid = int(arg_baseMolecule.residueHead[rid])
		while iAtom_in_rid != -1:

			for iDih in arg_baseMolecule.dihedralTable[iAtom_in_rid]:
				# see if it is a side chain dihedral exclusive to this resid 
				if iDih.is_from_same_residue() == rid and iDih.is_heavy_dihedral():
					dihNode = CustomDataTypes.TreeNode(None, None, iDih)
					diheds_in_rid.add_node(dihNode)

			iAtom_in_rid = int(arg_baseMolecule.atomArray[iAtom_in_rid])

		Utils.printflush('Found {} exclusive dihedrals in residue {}{}'.\
			format(len(diheds_in_rid), arg_baseMolecule.resnameArray[rid], arg_baseMolecule.residArray[rid]))

		# create an object of Class  ormationEntity corresponding to this residue
		newEntity = ConformationEntity.ConformationEntity(arg_order = len(diheds_in_rid), arg_numFrames = numFrames)

		# also initialize a string array that will store the state in each frame as a distinct string
		# made from coalesced character cast of numeric arrays
		ridDecimalReprArray = []

		# for each dihedral identified, get the state vector
		for i, iDih in enumerate(diheds_in_rid.list_in_order()):
			stateTS = iDih.get_state_ts(arg_dataContainer = arg_hostDataContainer, arg_verbose = arg_verbose)
			newEntity.timeSeries[i,:] = stateTS

		# Now coalesce integer labels of the constituent dihedrals in each time point to get 
		# an expression of the conformation at that time.
		for iFrame in range(numFrames):
			ridDecimalReprArray.append(Utils.coalesce_numeric_array(newEntity.timeSeries[:,iFrame]))


		# for each of the unique state get their count and compute the topographical entropy for this residue
		setOfstates = set(ridDecimalReprArray)
		Utils.printflush('Found {} dihedrals which collectively acquire {} unique conformers'.format(len(diheds_in_rid), len(setOfstates)))

		# print(ridDecimalReprArray)

		# total SC entropy at the topographical level of this residue
		ridTopogEntropy = 0.

		for iState in setOfstates:
			iCount = ridDecimalReprArray.count(iState) 

			# p Log(p) for this state
			iPlogP = iCount * (nmp.log(iCount) - logNumFrames)
			ridTopogEntropy += iPlogP;

		ridTopogEntropy /= numFrames;
		ridTopogEntropy *= -CONST.GAS_CONST  #(R)

		# Final residue SC information	
		Utils.printflush('{:<40s} : {:.4f}'.format('Residue Topographical Entropy from AEM ({} {})'.format(arg_baseMolecule.resnameArray[rid], arg_baseMolecule.residArray[rid]), ridTopogEntropy))

		Utils.printOut(arg_outFile, '{:<40s} : {:.4f}'.format('Residue Topographical Entropy from AEM ({} {})'.format(arg_baseMolecule.resnameArray[rid], arg_baseMolecule.residArray[rid]), ridTopogEntropy))

		# add this residue's SC entropy to the total SC entropy
		totalTopogEntropySC += ridTopogEntropy
			
	# total SC topographical entropy
	Utils.printflush('_'*60)
	Utils.printflush('{:<40} : {:>15.3f}'.format('Total Topog. Entropy (AEM) ', totalTopogEntropySC))
	Utils.printflush('-'*60)

	Utils.printOut(arg_outFile, '_'*60)
	Utils.printOut(arg_outFile, '{:<40} : {:>15.3f}'.format('Total Topog. Entropy (AEM)', totalTopogEntropySC))
	Utils.printOut(arg_outFile, '-'*60)

	return
#END

if __name__ == "__main__":

	helpm = "This script uses the CodeEntropy API to compute MCC \n \
				entropy of a molecule at its united atom (UA) level.\n\
				CodeEntropy package should be installed for this script to work.\n\
				Download the package from: https://github.com/arghya90/CodeEntropy.git\n\n\
				The instructions to install it are in the README.md file.\n\n\
				NOTE: READ THE CODE CAREFULLY TO UNDERSTAND WHAT IS HAPPENING.\n\n"

	syntax = "Provide inputs to `python /PATH/TO/mcc_api_UA.py` in this exact order:\n\n\
			  <trajfile with forces (.trr)> \n\
			  <topol file (.tpr)> \n\
			  <output file name> \n\
			  <beginTime (ps)>\n\
			  <endtime (ps)>\n\
			  <stride (integer; skip frames every so often)>\n\
			  <verbosity (0-5; how much noisy do you want the output to be>\n\
			  <fscale (e.g. 0.5 for force halving)>\n\
			  <tscale (e.g. 0.5 for torque halving)>\n\
			  <temperature (K; for entropy calculations)>"

	if len(sys.argv) != 11:
		Utils.printflush("Enough inputs not provided ({}). See below:\n\n".format(len(sys.argv)))
		Utils.printflush(helpm)
		Utils.printflush(syntax)
		sys.exit("Exiting...")

	# assemble inputs
	trajFile = sys.argv[1]
	topolFile = sys.argv[2]
	outFile = sys.argv[3]
	beginTime = float(sys.argv[4])
	endTime = float(sys.argv[5])
	stride=int(sys.argv[6])
	verbose = int(sys.argv[7])
	fScale = float(sys.argv[8])
	tScale = float(sys.argv[9])
	temper= float(sys.argv[10])

	# prepare outfile for writing
	Writer.write_file(outFile)

	# read gromacs traj and topol files to populate data structures
	mainMolecule, mainContainer = GromacsReader.read_gromacs_input(arg_beginTime=beginTime, \
														  arg_endTime=endTime, \
														  arg_outFile=outFile, \
														  arg_stride=stride, \
														  arg_tprFile=topolFile, \
														  arg_trajFile=trajFile, \
														  arg_verbose=verbose)

	# compute vib entropy
	get_svib_UA(arg_hostDataContainer = mainContainer, \
					arg_baseMolecule = mainMolecule, \
					arg_outFile = outFile, \
					arg_fScale = fScale, \
					arg_tScale = tScale, \
					arg_verbose = verbose, \
					arg_moutFile = None, \
					arg_nmdFile = None, \
					arg_temper = temper)

	# compute topographical entropy
	get_stopo_UA(arg_baseMolecule = mainMolecule, \
		                              arg_hostDataContainer = mainContainer,\
		                              arg_outFile = outFile, \
		                              arg_verbose = verbose)


