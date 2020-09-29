package uk.ac.liverpool.metfrag;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.TreeSet;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.smarts.SmartsPattern;

@SuppressWarnings("deprecation")
public class Fragmenter {
	
	/**
	 * 
	 */
	private final static String[] SMART_PATTERNS = { "O", "C(=O)O", "N", "C[Si](C)(C)O", "C[Si](C)C", "CO", "CN" }; //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$ //$NON-NLS-4$ //$NON-NLS-5$ //$NON-NLS-6$ //$NON-NLS-7$

	/**
	 * 
	 */
	private final SmartsPattern[] smartsPatterns;
	
	/**
	 * 
	 */
	private List<Integer> brokenBondToNeutralLossIndex = new ArrayList<>();
	
	/**
	 * 
	 */
	private FastBitArray[][] detectedNeutralLosses;
	
	/**
	 * 
	 */
	private byte maximumNumberOfAFragmentAddedToQueue = 2;
	
	/**
	 * 
	 */
	private final short[] minimumNumberImplicitHydrogens = { 1, 1, 2, 9, 9, 1, 0 };
	
	/**
	 * 
	 */
	private List<Integer> neutralLossIndex = new ArrayList<>();
	
	/**
	 * 
	 */
	private Precursor precursor;

	/**
	 * 
	 */
	private FastBitArray ringBondFastBitArray;
	
	/**
	 * 
	 */
	private boolean ringBondsInitialised = false;
	
	/**
	 * 
	 * @param smiles
	 * @throws CDKException 
	 * @throws IOException 
	 * @throws Exception
	 */
	public Fragmenter(final String smiles) throws CDKException, IOException {
		this.precursor = new Precursor(smiles);
		this.ringBondFastBitArray = new FastBitArray(this.precursor.getAtomContainer().getBondCount(), false);
		this.smartsPatterns = new SmartsPattern[SMART_PATTERNS.length];

		for (int i = 0; i < this.smartsPatterns.length; i++) {
			this.smartsPatterns[i] = SmartsPattern.create(SMART_PATTERNS[i]);
		}
		
		this.detectedNeutralLosses = getMatchingAtoms(this.precursor);
	}

	/**
	 * 
	 * @return Collection<List<Object>>
	 * @throws Exception
	 */
	public Collection<Fragment> getFragments(final int maxTreeDepth) throws Exception {
		final Collection<Fragment> fragments = new TreeSet<>();
		Queue<Fragment> fragmentsQueue = new LinkedList<>();

		final Fragment precursorFragment = new Fragment(this.precursor);
		fragmentsQueue.add(precursorFragment);
		fragments.add(precursorFragment);

		for (int k = 1; k <= maxTreeDepth; k++) {
			final Queue<Fragment> newFragmentsQueue = new LinkedList<>();

			while (!fragmentsQueue.isEmpty()) {
				final Fragment fragment = fragmentsQueue.poll();

				for (final Fragment childFragment : getFragmentsOfNextTreeDepth(fragment)) {
					fragments.add(childFragment);

					if (maxTreeDepth > 0) {
						newFragmentsQueue.add(childFragment);
					}
				}
			}

			fragmentsQueue = newFragmentsQueue;
		}

		return fragments;
	}
	
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

	/**
	 * return true if neutral loss has been detected before true is returned mass as
	 * well as molecular formula is modified
	 * 
	 * @param newGeneratedTopDownFragments
	 * @return
	 * @throws Exception
	 */
	private boolean checkForNeutralLossesAdaptMolecularFormulas(Fragment[] newGeneratedTopDownFragments,
			int removedBondIndex) throws Exception {

		byte neutralLossFragment = -1;
		for (int i = 0; i < this.detectedNeutralLosses.length; i++) {
			for (int ii = 0; ii < this.detectedNeutralLosses[i].length; ii++) {
				if (newGeneratedTopDownFragments[0].getAtomsArray().equals(this.detectedNeutralLosses[i][ii])) {
					/*
					 * check for previous broken bonds caused by neutral loss
					 */
					int[] brokenBondIndeces = newGeneratedTopDownFragments[1].getBrokenBondsArray().getSetIndices();
					for (int s = 0; s < brokenBondIndeces.length; s++) {
						int index = this.brokenBondToNeutralLossIndex.indexOf(Integer.valueOf(brokenBondIndeces[s]));
						if ((short) brokenBondIndeces[s] == removedBondIndex) {
							if (index == -1) {
								this.brokenBondToNeutralLossIndex.add(Integer.valueOf(removedBondIndex));
								this.neutralLossIndex.add(Integer.valueOf(i));
							}
							continue;
						}
					}
					return true;
				} else if (newGeneratedTopDownFragments[1].getAtomsArray().equals(this.detectedNeutralLosses[i][ii])) {

					/*
					 * check for previous broken bonds caused by neutral loss
					 */
					int[] brokenBondIndeces = newGeneratedTopDownFragments[0].getBrokenBondsArray().getSetIndices();
					for (int s = 0; s < brokenBondIndeces.length; s++) {
						int index = this.brokenBondToNeutralLossIndex.indexOf(Integer.valueOf(brokenBondIndeces[s]));
						if ((short) brokenBondIndeces[s] == removedBondIndex) {
							if (index == -1) {
								this.brokenBondToNeutralLossIndex.add(Integer.valueOf(removedBondIndex));
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
			Fragment precursorFragment, java.util.Queue<Fragment> toProcess, FastBitArray ringBondArray,
			java.util.Queue<Integer> lastCuttedRingBond) throws Exception {
		/*
		 * process all fragments that have been cutted in a ring without generating a
		 * new one
		 */
		while (!toProcess.isEmpty() && lastCuttedRingBond.size() != 0) {
			/*
			 * 
			 */
			Fragment currentFragment = toProcess.poll();
			int nextRingBondToCut = lastCuttedRingBond.poll().intValue() + 1;
			/*
			 * 
			 */
			for (int currentBond = nextRingBondToCut; currentBond < ringBondArray.getSize(); currentBond++) {
				if (!ringBondArray.get(currentBond))
					continue;
				if (currentFragment.getBrokenBondsArray().get(currentBond))
					continue;
				Fragment[] newFragments = { currentFragment };
				int[] connectedAtomIndeces = this.getConnectedAtomIndicesFromBondIndex(currentBond);
				newFragments = currentFragment.fragment(this.precursor, currentBond, connectedAtomIndeces);

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
						lastCuttedRingBond.add(Integer.valueOf(currentBond));
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
		int lastSkippedBonds = precursorFragment.getLastSkippedBond();
		short lastCuttedBond = (short) (precursorFragment.getBrokenBondsArray().getLastSetBit());
		if (lastSkippedBonds == -1)
			return;
		for (int currentBond = lastSkippedBonds; currentBond < lastCuttedBond; currentBond++) {

			if (this.ringBondFastBitArray.get(currentBond))
				continue;
			if (!precursorFragment.getBondsArray().get(currentBond))
				continue;

			int[] connectedAtomIndeces = this.getConnectedAtomIndicesFromBondIndex(currentBond);

			final Fragment[] newFragments = precursorFragment.fragment(this.precursor, currentBond,
					connectedAtomIndeces);

			Fragmenter.processGeneratedFragments(newFragments);
			if (newFragments.length == 2) {
				this.checkForNeutralLossesAdaptMolecularFormulas(newFragments, currentBond);
			} else {
				System.err.println("problem generating fragments"); //$NON-NLS-1$
				System.exit(1);
			}

			for (int k = 0; k < newFragments.length; k++) {
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
		FastBitArray ringBonds = new FastBitArray(precursorFragment.getBondsArray().getSize(), false);
		java.util.Queue<Fragment> ringBondCuttedFragments = new java.util.LinkedList<>();
		java.util.Queue<Integer> lastCuttedBondOfRing = new java.util.LinkedList<>();
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
		int nextBrokenIndexBondIndexToRemove = precursorFragment.getBrokenBondsArray().getLastSetBit() + 1;
		/*
		 * start from the last broken bond index
		 */
		for (int i = nextBrokenIndexBondIndexToRemove; i < precursorFragment.getBondsArray().getSize(); i++) {
			if (!precursorFragment.getBondsArray().get(i))
				continue;
			final int[] indecesOfBondConnectedAtoms = this.getConnectedAtomIndicesFromBondIndex(i);

			// try to generate at most two fragments by the removal of the given bond
			Fragment[] newGeneratedTopDownFragments = precursorFragment.fragment(this.precursor, i,
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
				lastCuttedBondOfRing.add(Integer.valueOf(i));
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
	 * @throws IOException 
	 */
	private FastBitArray[][] getMatchingAtoms(Precursor precursorMolecule) throws CDKException, IOException {
		final List<FastBitArray[]> matchedNeutralLossTypes = new ArrayList<>();

		for (byte i = 0; i < this.smartsPatterns.length; i++) {
			if (this.smartsPatterns[i].matches(precursorMolecule.getAtomContainer())) {
				/*
				 * get atom indeces containing to a neutral loss
				 */
				final int[][] matchingAtoms = this.smartsPatterns[i].matchAll(precursorMolecule.getAtomContainer()).toArray();
				/*
				 * store which is a valid loss based on the number of hydrogens
				 */
				boolean[] validMatches = new boolean[matchingAtoms.length];
				FastBitArray[] allMatches = new FastBitArray[matchingAtoms.length];
				int numberOfValidNeutralLosses = 0;
				/*
				 * check each part that is marked as neutral loss
				 */
				for (int ii = 0; ii < matchingAtoms.length; ii++) {
					int[] part = matchingAtoms[ii];
					/*
					 * count number of implicit hydrogens of this neutral loss
					 */
					int numberImplicitHydrogens = 0;
					allMatches[ii] = new FastBitArray(precursorMolecule.getAtomContainer().getAtomCount(), false);
					/*
					 * check all atoms
					 */
					for (int iii = 0; iii < part.length; iii++) {
						allMatches[ii].set(part[iii]);
						/*
						 * count number of implicit hydrogens of this neutral loss
						 */
						numberImplicitHydrogens += this.getImplicitHydrogenCount(part[iii]);
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
	
	/**
	 * Returns atom indices that are connected to bond with bondIdx.
	 * 
	 * @param bondIdx
	 * @return int[]
	 */
	private int[] getConnectedAtomIndicesFromBondIndex(final int bondIdx) {
		final IAtomContainer atomContainer = this.precursor.getAtomContainer();
		final IBond bond = atomContainer.getBond(bondIdx);
		return new int[] { atomContainer.indexOf(bond.getAtom(0)), atomContainer.indexOf(bond.getAtom(1)) };
	}
	
	/**
	 * 
	 * @param atomIdx
	 * @return int
	 */
	private int getImplicitHydrogenCount(final int atomIdx) {
		final IAtomContainer atomContainer = this.precursor.getAtomContainer();
		return atomContainer.getAtom(atomIdx).getImplicitHydrogenCount().intValue();
	}
}