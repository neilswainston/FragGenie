package uk.ac.liverpool.metfraglib.fragmenter;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;

import de.ipbhalle.metfraglib.FastBitArray;
import uk.ac.liverpool.metfraglib.fragment.BitArrayNeutralLoss;
import de.ipbhalle.metfraglib.parameter.Constants;
import uk.ac.liverpool.metfraglib.precursor.Precursor;
import uk.ac.liverpool.metfraglib.fragment.Fragment;

public class Fragmenter {

	private Double minimumMassDeviationForFragmentGeneration = Constants.DEFAULT_MIN_MASS_DEV_FOR_FRAGMENT_GENERATION;
	private byte maximumNumberOfAFragmentAddedToQueue = 2;
	private boolean ringBondsInitialised = false;
	private FastBitArray ringBondFastBitArray;
	private Precursor precursor;
	private double minimumFragmentMassLimit = 0.0;
	private BitArrayNeutralLoss[] detectedNeutralLosses;
	private List<Short> brokenBondToNeutralLossIndex = new ArrayList<>();
	private List<Integer> neutralLossIndex = new ArrayList<>();
	
	private final String[] smartPatterns = {"O", "C(=O)O", "N", "C[Si](C)(C)O", "C[Si](C)C", "CO", "CN"};
	private final short[] minimumNumberImplicitHydrogens = {1, 1, 2, 9, 9, 1, 0};
	

	/**
	 * 
	 * @param precursor
	 * @param maximumTreeDepth
	 * @throws Exception
	 */
	public Fragmenter(final Precursor precursor) throws Exception {
		this.precursor = precursor;
		this.ringBondFastBitArray = new FastBitArray(this.precursor.getNonHydrogenBondCount(), false);
		this.detectedNeutralLosses = getMatchingAtoms(this.precursor);
	}

	/**
	 * return true if neutral loss has been detected before true is returned mass as
	 * well as molecular formula is modified
	 * 
	 * @param newGeneratedTopDownFragments
	 * @return
	 */
	private boolean checkForNeutralLossesAdaptMolecularFormulas(Fragment[] newGeneratedTopDownFragments,
			short removedBondIndex) {
		if (newGeneratedTopDownFragments.length != 2) {
			System.err.println("Error: Cannot check for neutral losses for these fragments."); //$NON-NLS-1$
			return false;
		}
		byte neutralLossFragment = -1;
		for (int i = 0; i < this.detectedNeutralLosses.length; i++) {
			for (int ii = 0; ii < this.detectedNeutralLosses[i].getNumberNeutralLosses(); ii++) {
				if (newGeneratedTopDownFragments[0].getAtomsFastBitArray()
						.equals(this.detectedNeutralLosses[i].getNeutralLossAtomFastBitArray(ii))) {
					newGeneratedTopDownFragments[1].getMolecularFormula(this.precursor)
							.setNumberHydrogens((short) (newGeneratedTopDownFragments[1]
									.getMolecularFormula(this.precursor).getNumberHydrogens()
									+ this.detectedNeutralLosses[i].getHydrogenDifference()));
					/*
					 * check for previous broken bonds caused by neutral loss
					 */
					int[] brokenBondIndeces = newGeneratedTopDownFragments[1].getBrokenBondIndeces();
					for (int s = 0; s < brokenBondIndeces.length; s++) {
						int index = this.brokenBondToNeutralLossIndex.indexOf((short) brokenBondIndeces[s]);
						if ((short) brokenBondIndeces[s] == removedBondIndex) {
							if (index == -1) {
								this.brokenBondToNeutralLossIndex.add(removedBondIndex);
								this.neutralLossIndex.add(i);
							}
							continue;
						}
						if (index != -1) {
							newGeneratedTopDownFragments[1].getMolecularFormula(this.precursor)
									.setNumberHydrogens((short) (newGeneratedTopDownFragments[1]
											.getMolecularFormula(this.precursor).getNumberHydrogens()
											+ this.detectedNeutralLosses[this.neutralLossIndex.get(index)]
													.getHydrogenDifference()));
						}
					}
					return true;
				} else if (newGeneratedTopDownFragments[1].getAtomsFastBitArray()
						.equals(this.detectedNeutralLosses[i].getNeutralLossAtomFastBitArray(ii))) {
					newGeneratedTopDownFragments[0].getMolecularFormula(this.precursor)
							.setNumberHydrogens((short) (newGeneratedTopDownFragments[0]
									.getMolecularFormula(this.precursor).getNumberHydrogens()
									+ this.detectedNeutralLosses[i].getHydrogenDifference()));
					// newGeneratedTopDownFragments[0].setTreeDepth((byte)(newGeneratedTopDownFragments[0].getTreeDepth()
					// - 1));
					/*
					 * check for previous broken bonds caused by neutral loss
					 */
					int[] brokenBondIndeces = newGeneratedTopDownFragments[0].getBrokenBondIndeces();
					for (int s = 0; s < brokenBondIndeces.length; s++) {
						int index = this.brokenBondToNeutralLossIndex.indexOf((short) brokenBondIndeces[s]);
						if ((short) brokenBondIndeces[s] == removedBondIndex) {
							if (index == -1) {
								this.brokenBondToNeutralLossIndex.add(removedBondIndex);
								this.neutralLossIndex.add(i);
							}
							continue;
						}
						if (index != -1) {
							newGeneratedTopDownFragments[0].getMolecularFormula(this.precursor)
									.setNumberHydrogens((short) (newGeneratedTopDownFragments[0]
											.getMolecularFormula(this.precursor).getNumberHydrogens()
											+ this.detectedNeutralLosses[this.neutralLossIndex.get(index)]
													.getHydrogenDifference()));
						}
					}
					return true;
				}
			}
		}
		if (neutralLossFragment == -1)
			return false;
		return true;
	}

	/**
	 * generates all fragments of the given precursor fragment to reach the new tree
	 * depth
	 * 
	 * @throws Exception
	 */
	public ArrayList<Fragment> getFragmentsOfNextTreeDepth(Fragment precursorFragment) throws Exception {
		FastBitArray ringBonds = new FastBitArray(precursorFragment.getBondsFastBitArray().getSize(), false);
		java.util.Queue<Fragment> ringBondCuttedFragments = new java.util.LinkedList<>();
		java.util.Queue<Short> lastCuttedBondOfRing = new java.util.LinkedList<>();
		ArrayList<Fragment> fragmentsOfNextTreeDepth = new ArrayList<>();
		/*
		 * generate fragments of skipped bonds
		 */
		if (this.ringBondsInitialised)
			this.generateFragmentsOfSkippedBonds(fragmentsOfNextTreeDepth, precursorFragment);
		/*
		 * get the last bond index that was removed; from there on the next bonds will
		 * be removed
		 */
		short nextBrokenIndexBondIndexToRemove = (short) (precursorFragment.getMaximalIndexOfRemovedBond() + 1);
		/*
		 * start from the last broken bond index
		 */
		for (short i = nextBrokenIndexBondIndexToRemove; i < precursorFragment.getBondsFastBitArray().getSize(); i++) {
			if (!precursorFragment.getBondsFastBitArray().get(i))
				continue;
			short[] indecesOfBondConnectedAtoms = this.precursor.getConnectedAtomIndecesOfBondIndex(i);
			/*
			 * try to generate at most two fragments by the removal of the given bond
			 */
			Fragment[] newGeneratedTopDownFragments = precursorFragment.traverseMolecule(this.precursor, i,
					indecesOfBondConnectedAtoms);
			/*
			 * in case the precursor wasn't splitted try to cleave an additional bond until
			 * 
			 * 1. two fragments are generated or 2. the maximum number of trials have been
			 * reached 3. no further bond can be removed
			 */
			if (newGeneratedTopDownFragments.length == 1) {
				ringBonds.set(i, true);
				newGeneratedTopDownFragments[0].setLastSkippedBond((short) (i + 1));
				ringBondCuttedFragments.add(newGeneratedTopDownFragments[0]);
				lastCuttedBondOfRing.add(i);
				if (!this.ringBondsInitialised)
					this.ringBondFastBitArray.set(i);
			}
			/*
			 * pre-processing of the generated fragment/s
			 */
			this.processGeneratedFragments(newGeneratedTopDownFragments);
			/*
			 * if two new fragments have been generated set them as valid
			 */
			if (newGeneratedTopDownFragments.length == 2) {
				this.checkForNeutralLossesAdaptMolecularFormulas(newGeneratedTopDownFragments, i);
				newGeneratedTopDownFragments[0].setAsValidFragment();
				newGeneratedTopDownFragments[1].setAsValidFragment();
			}
			/*
			 * add fragment/s to vector after setting the proper precursor
			 */
			for (int k = 0; k < newGeneratedTopDownFragments.length; k++) {
				// precursorFragment.addChild(newGeneratedTopDownFragments[k]);
				if (newGeneratedTopDownFragments.length == 2)
					fragmentsOfNextTreeDepth.add(newGeneratedTopDownFragments[k]);
				/*
				 * if (precursorFragment.isValidFragment()) {
				 * newGeneratedTopDownFragments[k].setPrecursorFragment(precursorFragment); }
				 * else {
				 * newGeneratedTopDownFragments[k].setPrecursorFragment(precursorFragment.
				 * hasPrecursorFragment() ? precursorFragment.getPrecursorFragment() :
				 * precursorFragment); }
				 */

			}
		}
		/*
		 * create fragments by ring bond cleavage and store them in the given vector
		 */
		this.createRingBondCleavedFragments(fragmentsOfNextTreeDepth, precursorFragment, ringBondCuttedFragments,
				ringBonds, lastCuttedBondOfRing);
		this.ringBondsInitialised = true;

		return fragmentsOfNextTreeDepth;
	}

	/**
	 * @throws Exception
	 * 
	 */

	/*
	 * generate fragments by removing bonds that were skipped due to ring bond
	 * cleavage
	 */
	private void generateFragmentsOfSkippedBonds(ArrayList<Fragment> newGeneratedTopDownFragments,
			Fragment precursorFragment) throws Exception {
		short lastSkippedBonds = precursorFragment.getLastSkippedBond();
		short lastCuttedBond = (short) (precursorFragment.getMaximalIndexOfRemovedBond());
		if (lastSkippedBonds == -1)
			return;
		for (short currentBond = lastSkippedBonds; currentBond < lastCuttedBond; currentBond++) {

			if (this.ringBondFastBitArray.get(currentBond))
				continue;
			if (!precursorFragment.getBondsFastBitArray().get(currentBond))
				continue;

			short[] connectedAtomIndeces = this.precursor.getConnectedAtomIndecesOfBondIndex(currentBond);

			Fragment[] newFragments = precursorFragment.traverseMolecule(this.precursor, currentBond,
					connectedAtomIndeces);

			this.processGeneratedFragments(newFragments);
			if (newFragments.length == 2) {
				this.checkForNeutralLossesAdaptMolecularFormulas(newFragments, currentBond);
				newFragments[0].setAsValidFragment();
				newFragments[1].setAsValidFragment();
			} else {
				System.err.println("problem generating fragments"); //$NON-NLS-1$
				System.exit(1);
			}

			for (int k = 0; k < newFragments.length; k++) {
				// precursorFragment.addChild(newFragments[k]);
				/*
				 * if (precursorFragment.isValidFragment())
				 * newFragments[k].setPrecursorFragment(precursorFragment); else
				 * newFragments[k].setPrecursorFragment(precursorFragment.hasPrecursorFragment()
				 * ? precursorFragment.getPrecursorFragment() : precursorFragment);
				 */
				if (newFragments.length == 2) {
					newGeneratedTopDownFragments.add(newFragments[k]);
				}
			}
		}
	}

	/**
	 * 
	 * @param newGeneratedTopDownFragments
	 */
	private void processGeneratedFragments(Fragment[] newGeneratedTopDownFragments) {
		if (newGeneratedTopDownFragments.length == 2) {
			newGeneratedTopDownFragments[0].setAddedToQueueCounts((byte) 1);
			newGeneratedTopDownFragments[1].setAddedToQueueCounts((byte) 1);
			if (newGeneratedTopDownFragments[0].getMonoisotopicMass(this.precursor) <= this.minimumFragmentMassLimit
					- this.minimumMassDeviationForFragmentGeneration) {
				newGeneratedTopDownFragments[0].setAsDiscardedForFragmentation();
			}
			if (newGeneratedTopDownFragments[1].getMonoisotopicMass(this.precursor) <= this.minimumFragmentMassLimit
					- this.minimumMassDeviationForFragmentGeneration) {
				newGeneratedTopDownFragments[1].setAsDiscardedForFragmentation();
			}
		}
	}

	/**
	 * 
	 * @param newGeneratedTopDownFragments
	 * @param precursorFragment
	 * @param toProcess
	 * @param ringBondFastBitArray
	 * @param lastCuttedRingBond
	 * @return
	 * @throws Exception
	 */
	private ArrayList<Fragment> createRingBondCleavedFragments(ArrayList<Fragment> newGeneratedTopDownFragments,
			Fragment precursorFragment, java.util.Queue<Fragment> toProcess, FastBitArray ringBondFastBitArray,
			java.util.Queue<Short> lastCuttedRingBond) throws Exception {
		/*
		 * process all fragments that have been cutted in a ring without generating a
		 * new one
		 */
		while (!toProcess.isEmpty() && lastCuttedRingBond.size() != 0) {
			/*
			 * 
			 */
			Fragment currentFragment = toProcess.poll();
			short nextRingBondToCut = (short) (lastCuttedRingBond.poll() + 1);
			/*
			 * 
			 */
			for (short currentBond = nextRingBondToCut; currentBond < ringBondFastBitArray.getSize(); currentBond++) {
				if (!ringBondFastBitArray.get(currentBond))
					continue;
				if (currentFragment.getBrokenBondsFastBitArray().get(currentBond))
					continue;
				Fragment[] newFragments = { currentFragment };
				short[] connectedAtomIndeces = this.precursor.getConnectedAtomIndecesOfBondIndex(currentBond);
				newFragments = currentFragment.traverseMolecule(this.precursor, currentBond, connectedAtomIndeces);

				//
				// pre-processing of the generated fragment/s
				//
				this.processGeneratedFragments(newFragments);
				//
				// if two new fragments have been generated set them as valid
				//
				if (newFragments.length == 2) {
					newFragments[0].setAsValidFragment();
					newFragments[1].setAsValidFragment();
					newFragments[0].setLastSkippedBond(currentFragment.getLastSkippedBond());
					newFragments[1].setLastSkippedBond(currentFragment.getLastSkippedBond());
				}
				//
				// set precursor fragment of generated fragment(s) and the child(ren) of
				// precursor fragments
				//
				for (int k = 0; k < newFragments.length; k++) {
					if (newFragments.length == 2) {
						newGeneratedTopDownFragments.add(newFragments[k]);
					}
				}

				if (newFragments.length == 1) {
					if (newFragments[0].getAddedToQueueCounts() < this.maximumNumberOfAFragmentAddedToQueue) {
						toProcess.add(newFragments[0]);
						lastCuttedRingBond.add(currentBond);
					} else {
						newGeneratedTopDownFragments.add(newFragments[0]);
					}
				}
			}
		}

		return newGeneratedTopDownFragments;
	}
	
	/**
	 * 
	 * @param precursorMolecule
	 * @return
	 */
	public BitArrayNeutralLoss[] getMatchingAtoms(Precursor precursorMolecule) {
		SMARTSQueryTool[] smartsQuerytools = new SMARTSQueryTool[smartPatterns.length];
		for(int i = 0; i < smartsQuerytools.length; i++) {
			smartsQuerytools[i] = new SMARTSQueryTool(smartPatterns[i], DefaultChemObjectBuilder.getInstance());
		}
		java.util.ArrayList<BitArrayNeutralLoss> matchedNeutralLossTypes = new java.util.ArrayList<>();
		for(byte i = 0; i < smartsQuerytools.length; i++) {
			try {
				if(smartsQuerytools[i].matches(precursorMolecule.getStructureAsIAtomContainer())) {
					/*
					 * get atom indeces containing to a neutral loss
					 */
					java.util.List<java.util.List<Integer>> matchingAtoms = smartsQuerytools[i].getMatchingAtoms();
					/*
					 * store which is a valid loss based on the number of hydrogens
					 */
					boolean[] validMatches = new boolean[matchingAtoms.size()];
					FastBitArray[] allMatches = new FastBitArray[matchingAtoms.size()];
					int numberOfValidNeutralLosses = 0;
					/*
					 * check each part that is marked as neutral loss
					 */
					for(int ii = 0; ii < matchingAtoms.size(); ii++) {
						java.util.List<Integer> part = matchingAtoms.get(ii);
						/*
						 * count number of implicit hydrogens of this neutral loss
						 */
						int numberImplicitHydrogens = 0;
						allMatches[ii] = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount());
						/*
						 * check all atoms 
						 */
						for(int iii = 0; iii < part.size(); iii++) {
							allMatches[ii].set(part.get(iii));
							/*
							 * count number of implicit hydrogens of this neutral loss
							 */
							numberImplicitHydrogens += precursorMolecule.getNumberHydrogensConnectedToAtomIndex(part.get(iii));
						}
						/*
						 * valid neutral loss match if number implicit hydrogens are at least the number of hydrogens
						 * needed for the certain neutral loss
						 */
						if(numberImplicitHydrogens >= minimumNumberImplicitHydrogens[i]) {
							validMatches[ii] = true;
							numberOfValidNeutralLosses++;
						}
					}
					/*
					 * create BitArrayNeutralLosses of valid neutral loss part detections
					 */
					if(numberOfValidNeutralLosses != 0) {
						BitArrayNeutralLoss newDetectedNeutralLoss = 
							new BitArrayNeutralLoss(numberOfValidNeutralLosses, i);
						int neutralLossIndexOfBitArrayNeutralLoss = 0;
						for(int k = 0; k < validMatches.length; k++) {
							if(validMatches[k]) {
								newDetectedNeutralLoss.setNeutralLoss(neutralLossIndexOfBitArrayNeutralLoss, allMatches[k]);
								neutralLossIndexOfBitArrayNeutralLoss++;
							}
						}
						/*
						 * store them in vector
						 */
						matchedNeutralLossTypes.add(newDetectedNeutralLoss);
					}
				}
			} catch (CDKException e) {
				e.printStackTrace();
			}
		}
		BitArrayNeutralLoss[] matchedNeutralLossTypesArray = new BitArrayNeutralLoss[matchedNeutralLossTypes.size()];
		for(int i = 0; i < matchedNeutralLossTypes.size(); i++) {
			matchedNeutralLossTypesArray[i] = matchedNeutralLossTypes.get(i);
		}
		return matchedNeutralLossTypesArray;
	}
}
