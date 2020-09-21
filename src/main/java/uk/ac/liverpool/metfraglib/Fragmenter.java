package uk.ac.liverpool.metfraglib;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.TreeMap;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;

import uk.ac.liverpool.metfraglib.FastBitArray;
import uk.ac.liverpool.metfraglib.Fragment;
import uk.ac.liverpool.metfraglib.Precursor;

public class Fragmenter {

	/**
	 * 
	 * @param newGeneratedTopDownFragments
	 */
	private static void processGeneratedFragments(Fragment[] newGeneratedTopDownFragments) {
		if (newGeneratedTopDownFragments.length == 2) {
			newGeneratedTopDownFragments[0].setAddedToQueueCounts((byte) 1);
			newGeneratedTopDownFragments[1].setAddedToQueueCounts((byte) 1);
		}
	}

	private List<Short> brokenBondToNeutralLossIndex = new ArrayList<>();
	private FastBitArray[][] detectedNeutralLosses;
	private byte maximumNumberOfAFragmentAddedToQueue = 2;
	private final short[] minimumNumberImplicitHydrogens = { 1, 1, 2, 9, 9, 1, 0 };
	private List<Integer> neutralLossIndex = new ArrayList<>();
	private Precursor precursor;

	private FastBitArray ringBondFastBitArray;
	private boolean ringBondsInitialised = false;

	private final String[] smartPatterns = { "O", "C(=O)O", "N", "C[Si](C)(C)O", "C[Si](C)C", "CO", "CN" }; //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$ //$NON-NLS-4$ //$NON-NLS-5$ //$NON-NLS-6$ //$NON-NLS-7$

	/**
	 * 
	 * @param precursor
	 * @param maximumTreeDepth
	 * @throws Exception
	 */
	public Fragmenter(final String smiles) throws Exception {
		this.precursor = Precursor.fromSmiles(smiles);
		this.ringBondFastBitArray = new FastBitArray(this.precursor.getNonHydrogenBondCount(), false);
		this.detectedNeutralLosses = getMatchingAtoms(this.precursor);
	}
	
	/**
	 * 
	 * @return Map<String, Float>
	 * @throws Exception
	 */
	public Map<String, Float> getFormulaToMasses(final int maxTreeDepth) throws Exception {
		final Map<String, Float> formulaToMasses = new TreeMap<>();
		Queue<Fragment> fragments = new LinkedList<>();
		
		final Fragment precursorFragment = new Fragment(this.precursor);
		fragments.add(precursorFragment);
		
		formulaToMasses.put(precursorFragment.getFormula(), Float.valueOf(precursorFragment.getMonoisotopicMass()));
		
		for (int k = 1; k <= maxTreeDepth; k++) {
			Queue<Fragment> newFragments = new LinkedList<>();

			while (!fragments.isEmpty()) {
				final Fragment fragment = fragments.poll();

				for (final Fragment childFragment : getFragmentsOfNextTreeDepth(fragment)) {
					formulaToMasses.put(childFragment.getFormula(), Float.valueOf(childFragment.getMonoisotopicMass()));

					if (maxTreeDepth > 0) {
						newFragments.add(childFragment);
					}
				}
			}

			fragments = newFragments;
		}
		
		return formulaToMasses;
	}

	/**
	 * return true if neutral loss has been detected before true is returned mass as
	 * well as molecular formula is modified
	 * 
	 * @param newGeneratedTopDownFragments
	 * @return
	 * @throws Exception
	 */
	private boolean checkForNeutralLossesAdaptMolecularFormulas(Fragment[] newGeneratedTopDownFragments,
			short removedBondIndex) throws Exception {
		if (newGeneratedTopDownFragments.length != 2) {
			System.err.println("Error: Cannot check for neutral losses for these fragments."); //$NON-NLS-1$
			return false;
		}
		byte neutralLossFragment = -1;
		for (int i = 0; i < this.detectedNeutralLosses.length; i++) {
			for (int ii = 0; ii < this.detectedNeutralLosses[i].length; ii++) {
				if (newGeneratedTopDownFragments[0].getAtomsFastBitArray()
						.equals(this.detectedNeutralLosses[i][ii])) {
					/*
					 * check for previous broken bonds caused by neutral loss
					 */
					int[] brokenBondIndeces = newGeneratedTopDownFragments[1].getBrokenBondIndeces();
					for (int s = 0; s < brokenBondIndeces.length; s++) {
						int index = this.brokenBondToNeutralLossIndex.indexOf(Short.valueOf((short)brokenBondIndeces[s]));
						if ((short) brokenBondIndeces[s] == removedBondIndex) {
							if (index == -1) {
								this.brokenBondToNeutralLossIndex.add(Short.valueOf(removedBondIndex));
								this.neutralLossIndex.add(Integer.valueOf(i));
							}
							continue;
						}
					}
					return true;
				} else if (newGeneratedTopDownFragments[1].getAtomsFastBitArray()
						.equals(this.detectedNeutralLosses[i][ii])) {

					/*
					 * check for previous broken bonds caused by neutral loss
					 */
					int[] brokenBondIndeces = newGeneratedTopDownFragments[0].getBrokenBondIndeces();
					for (int s = 0; s < brokenBondIndeces.length; s++) {
						int index = this.brokenBondToNeutralLossIndex.indexOf(Short.valueOf((short)brokenBondIndeces[s]));
						if ((short) brokenBondIndeces[s] == removedBondIndex) {
							if (index == -1) {
								this.brokenBondToNeutralLossIndex.add(Short.valueOf(removedBondIndex));
								this.neutralLossIndex.add(Integer.valueOf(i));
							}
							continue;
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
			short nextRingBondToCut = (short) (lastCuttedRingBond.poll().intValue() + 1);
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
				Fragmenter.processGeneratedFragments(newFragments);
				//
				// if two new fragments have been generated set them as valid
				//
				if (newFragments.length == 2) {
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
						lastCuttedRingBond.add(Short.valueOf(currentBond));
					} else {
						newGeneratedTopDownFragments.add(newFragments[0]);
					}
				}
			}
		}

		return newGeneratedTopDownFragments;
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

			Fragmenter.processGeneratedFragments(newFragments);
			if (newFragments.length == 2) {
				this.checkForNeutralLossesAdaptMolecularFormulas(newFragments, currentBond);
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
	 * generates all fragments of the given precursor fragment to reach the new tree
	 * depth
	 * 
	 * @throws Exception
	 */
	private ArrayList<Fragment> getFragmentsOfNextTreeDepth(Fragment precursorFragment) throws Exception {
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
			final short[] indecesOfBondConnectedAtoms = this.precursor.getConnectedAtomIndecesOfBondIndex(i);

			// Check bond strength:
			final List<Object> bond = this.precursor.getBond(i);
			System.out.println(bond);
			
			// try to generate at most two fragments by the removal of the given bond
			Fragment[] newGeneratedTopDownFragments = precursorFragment.traverseMolecule(this.precursor, i,
					indecesOfBondConnectedAtoms);
			/*
			 * in case the precursor wasn't splitted try to cleave an additional bond until
			 * 
			 * 1. two fragments are generated or 2. the maximum number of trials have been
			 * reached 3. no further bond can be removed
			 */
			if (newGeneratedTopDownFragments.length == 1) {
				ringBonds.set(i);
				newGeneratedTopDownFragments[0].setLastSkippedBond((short) (i + 1));
				ringBondCuttedFragments.add(newGeneratedTopDownFragments[0]);
				lastCuttedBondOfRing.add(Short.valueOf(i));
				if (!this.ringBondsInitialised)
					this.ringBondFastBitArray.set(i);
			}
			/*
			 * pre-processing of the generated fragment/s
			 */
			Fragmenter.processGeneratedFragments(newGeneratedTopDownFragments);
			/*
			 * if two new fragments have been generated set them as valid
			 */
			if (newGeneratedTopDownFragments.length == 2) {
				this.checkForNeutralLossesAdaptMolecularFormulas(newGeneratedTopDownFragments, i);
			}
			/*
			 * add fragment/s to vector after setting the proper precursor
			 */
			for (int k = 0; k < newGeneratedTopDownFragments.length; k++) {
				// precursorFragment.addChild(newGeneratedTopDownFragments[k]);
				if (newGeneratedTopDownFragments.length == 2)
					fragmentsOfNextTreeDepth.add(newGeneratedTopDownFragments[k]);
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
	 * 
	 * @param precursorMolecule
	 * @return
	 * @throws CDKException 
	 */
	private FastBitArray[][] getMatchingAtoms(Precursor precursorMolecule) throws CDKException {
		final SMARTSQueryTool[] smartsQuerytools = new SMARTSQueryTool[this.smartPatterns.length];
		
		for (int i = 0; i < smartsQuerytools.length; i++) {
			smartsQuerytools[i] = new SMARTSQueryTool(this.smartPatterns[i], DefaultChemObjectBuilder.getInstance());
		}
		
		final List<FastBitArray[]> matchedNeutralLossTypes = new ArrayList<>();
		
		for (byte i = 0; i < smartsQuerytools.length; i++) {
			if (smartsQuerytools[i].matches(precursorMolecule.getStructureAsIAtomContainer())) {
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
				for (int ii = 0; ii < matchingAtoms.size(); ii++) {
					java.util.List<Integer> part = matchingAtoms.get(ii);
					/*
					 * count number of implicit hydrogens of this neutral loss
					 */
					int numberImplicitHydrogens = 0;
					allMatches[ii] = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount());
					/*
					 * check all atoms
					 */
					for (int iii = 0; iii < part.size(); iii++) {
						allMatches[ii].set(part.get(iii).intValue());
						/*
						 * count number of implicit hydrogens of this neutral loss
						 */
						numberImplicitHydrogens += precursorMolecule
								.getNumberHydrogensConnectedToAtomIndex(part.get(iii).intValue());
					}
					/*
					 * valid neutral loss match if number implicit hydrogens are at least the number
					 * of hydrogens needed for the certain neutral loss
					 */
					if (numberImplicitHydrogens >= this.minimumNumberImplicitHydrogens[i]) {
						validMatches[ii] = true;
						numberOfValidNeutralLosses++;
					}
				}
				/*
				 * create BitArrayNeutralLosses of valid neutral loss part detections
				 */
				if (numberOfValidNeutralLosses != 0) {
					FastBitArray[] newDetectedNeutralLoss = new FastBitArray[numberOfValidNeutralLosses];
					int neutralLossIndexOfBitArrayNeutralLoss = 0;
					for (int k = 0; k < validMatches.length; k++) {
						if (validMatches[k]) {
							newDetectedNeutralLoss[neutralLossIndexOfBitArrayNeutralLoss] = allMatches[k];
							neutralLossIndexOfBitArrayNeutralLoss++;
						}
					}
					/*
					 * store them in vector
					 */
					matchedNeutralLossTypes.add(newDetectedNeutralLoss);
				}
			}
		}
		
		final FastBitArray[][] matchedNeutralLossTypesArray = new FastBitArray[matchedNeutralLossTypes.size()][];
		
		for (int i = 0; i < matchedNeutralLossTypes.size(); i++) {
			matchedNeutralLossTypesArray[i] = matchedNeutralLossTypes.get(i);
		}
		
		return matchedNeutralLossTypesArray;
	}
}