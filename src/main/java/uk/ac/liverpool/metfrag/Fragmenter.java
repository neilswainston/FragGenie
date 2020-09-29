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
import org.openscience.cdk.smiles.smarts.SmartsPattern;

@SuppressWarnings("deprecation")
public class Fragmenter {
	
	/**
	 * 
	 */
	private final static String[] SMARTS_PATTERNS_STRINGS = { "O", "C(=O)O", "N", "C[Si](C)(C)O", "C[Si](C)C", "CO", "CN" }; //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$ //$NON-NLS-4$ //$NON-NLS-5$ //$NON-NLS-6$ //$NON-NLS-7$

	/**
	 * 
	 */
	private final static int[] MIN_NUM_IMPLICIT_HYDROGENS = { 1, 1, 2, 9, 9, 1, 0 };
	
	/**
	 * 
	 */
	private final static SmartsPattern[] SMARTS_PATTERNS;
	
	/**
	 * 
	 */
	private List<List<boolean[]>> neutralLosses = new ArrayList<>();
	
	/**
	 * 
	 */
	private List<Integer> brokenBondToNeutralLossIndex = new ArrayList<>();
	

	
	/**
	 * 
	 */
	private byte maximumNumberOfAFragmentAddedToQueue = 2;
	
	/**
	 * 
	 */
	private List<Integer> neutralLossIndex = new ArrayList<>();
	
	/**
	 * 
	 */
	private IAtomContainer prec;

	/**
	 * 
	 */
	private boolean[] ringBondArray;
	
	/**
	 * 
	 */
	private boolean ringBondsInitialised = false;
	
	
	static {
		SMARTS_PATTERNS = new SmartsPattern[SMARTS_PATTERNS_STRINGS.length];

		for (int i = 0; i < SMARTS_PATTERNS_STRINGS.length; i++) {
			try {
				SMARTS_PATTERNS[i] = SmartsPattern.create(SMARTS_PATTERNS_STRINGS[i]);
			}
			catch (IOException e) {
				throw new ExceptionInInitializerError(e);
			}
		}
	}
	
	/**
	 * 
	 * @param precursor
	 * @throws CDKException
	 * @throws IOException
	 */
	public Fragmenter(final IAtomContainer precursor) throws IOException, CDKException {
		this.prec = precursor;
		this.ringBondArray = new boolean[this.prec.getBondCount()];
		initialiseNeutralLosses();
	}
	
	/**
	 * 
	 * @throws CDKException
	 * @throws IOException
	 */
	private void initialiseNeutralLosses() throws CDKException, IOException {
		int i = 0;
		
		for (SmartsPattern pattern : SMARTS_PATTERNS) {
			// Get atom indices corresponding to a neutral loss:
			final int[][] matchingAtoms = pattern.matchAll(this.prec).toArray();
			
			if (matchingAtoms.length > 0) {
				// Store valid losses based on the number of hydrogens:
				final boolean[] validMatches = new boolean[matchingAtoms.length];
				final boolean[][] allMatches = new boolean[matchingAtoms.length][];
				
				// Check each part that is marked as neutral loss:
				for (int ii = 0; ii < matchingAtoms.length; ii++) {
					final int[] part = matchingAtoms[ii];
					
					// Count number of implicit hydrogens of this neutral loss:
					int numberImplicitHydrogens = 0;
					allMatches[ii] = new boolean[this.prec.getAtomCount()];
					
					// Check all atoms:
					for (int atomIdx : part) {
						allMatches[ii][atomIdx] = true;
						// Count number of implicit hydrogens of this neutral loss:
						numberImplicitHydrogens += this.getImplicitHydrogenCount(atomIdx);
					}
					
					// Valid neutral loss match if number implicit hydrogens are at least the number
					// of hydrogens needed for the certain neutral loss:
					if (numberImplicitHydrogens >= MIN_NUM_IMPLICIT_HYDROGENS[i]) {
						validMatches[ii] = true;
					}
				}
				
				// Create new neutral loss of valid neutral loss part detections:
				final List<boolean[]> newNeutralLoss = new ArrayList<>();
				
				for (int k = 0; k < validMatches.length; k++) {
					if (validMatches[k]) {
						newNeutralLoss.add(allMatches[k]);
					}
				}
				
				this.neutralLosses.add(newNeutralLoss);
			}
			
			i++;
		}
	}

	/**
	 * 
	 * @return Collection<List<Object>>
	 */
	public Collection<Fragment> getFragments(final int maxTreeDepth) {
		final Collection<Fragment> fragments = new TreeSet<>();
		final Queue<Fragment> fragmentsQueue = new LinkedList<>();

		final Fragment precursorFragment = new Fragment(this.prec);
		fragmentsQueue.add(precursorFragment);
		fragments.add(precursorFragment);

		for (int k = 0; k < maxTreeDepth; k++) {
			final Collection<Fragment> newFragmentsQueue = new LinkedList<>();

			while (!fragmentsQueue.isEmpty()) {
				final Fragment fragment = fragmentsQueue.poll();
				final Collection<Fragment> childFragments = getFragmentsOfNextTreeDepth(fragment);
				fragments.addAll(childFragments);
				newFragmentsQueue.addAll(childFragments);
			}

			fragmentsQueue.addAll(newFragmentsQueue);
		}

		return fragments;
	}
	
	/**
	 * Generates all fragments of the given precursor fragment to reach the new tree depth.
	 * 
	 * @param parentFragment
	 * @return Collection<Fragment>
	 */
	private Collection<Fragment> getFragmentsOfNextTreeDepth(final Fragment parentFragment) {
		final Collection<Fragment> childFragments = new ArrayList<>();
		
		final boolean[] ringBonds = new boolean[parentFragment.getBondsArray().length];
		final Queue<Fragment> ringBondCuttedFragments = new LinkedList<>();
		final Queue<Integer> lastCuttedBondOfRing = new LinkedList<>();
		
		// Generate fragments of skipped bonds:
		if (this.ringBondsInitialised) {
			this.generateFragmentsOfSkippedBonds(childFragments, parentFragment);
		}
			
		// Get the last bond index that was removed; from there on the next bonds will be removed:
		final int nextBondIndexToRemove = getLastSetBit(parentFragment.getBrokenBondsArray()) + 1;
		
		/*
		 * start from the last broken bond index
		 */
		for (int bondIdx = nextBondIndexToRemove; bondIdx < parentFragment.getBondsArray().length; bondIdx++) {
			if (!parentFragment.getBondsArray()[bondIdx]) {
				continue;
			}
			
			// Try to generate at most two fragments by the breaking of the given bond:
			final Fragment[] newGeneratedTopDownFragments = parentFragment.fragment(bondIdx);
			
			/*
			 * in case the precursor wasn't splitted try to cleave an additional bond until
			 * 
			 * 1. two fragments are generated or 2. the maximum number of trials have been
			 * reached 3. no further bond can be removed
			 */
			if (newGeneratedTopDownFragments.length == 1) {
				ringBonds[bondIdx] = true;
				newGeneratedTopDownFragments[0].setLastSkippedBond(bondIdx + 1);
				ringBondCuttedFragments.add(newGeneratedTopDownFragments[0]);
				lastCuttedBondOfRing.add(Integer.valueOf(bondIdx));
				if (!this.ringBondsInitialised)
					this.ringBondArray[bondIdx] = true;
			}
			/*
			 * pre-processing of the generated fragment/s
			 */
			Fragmenter.processGeneratedFragments(newGeneratedTopDownFragments);
			/*
			 * if two new fragments have been generated set them as valid
			 */
			if (newGeneratedTopDownFragments.length == 2) {
				this.checkForNeutralLossesAdaptMolecularFormulas(newGeneratedTopDownFragments, bondIdx);
			}
			/*
			 * add fragment/s to vector after setting the proper precursor
			 */
			for (int k = 0; k < newGeneratedTopDownFragments.length; k++) {
				// precursorFragment.addChild(newGeneratedTopDownFragments[k]);
				if (newGeneratedTopDownFragments.length == 2)
					childFragments.add(newGeneratedTopDownFragments[k]);
			}
		}
		/*
		 * create fragments by ring bond cleavage and store them in the given vector
		 */
		this.createRingBondCleavedFragments(childFragments, parentFragment, ringBondCuttedFragments, ringBonds, lastCuttedBondOfRing);
		this.ringBondsInitialised = true;

		return childFragments;
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
	private boolean checkForNeutralLossesAdaptMolecularFormulas(Fragment[] newGeneratedTopDownFragments, int removedBondIndex) {

		byte neutralLossFragment = -1;
		for (int i = 0; i < this.neutralLosses.size(); i++) {
			for (int ii = 0; ii < this.neutralLosses.get(i).size(); ii++) {
				if (newGeneratedTopDownFragments[0].getAtomsArray().equals(this.neutralLosses.get(i).get(ii))) {
					/*
					 * check for previous broken bonds caused by neutral loss
					 */
					boolean[] brokenBondsArray = newGeneratedTopDownFragments[1].getBrokenBondsArray();
					
					for (int s = 0; s < brokenBondsArray.length; s++) {
						if(brokenBondsArray[s]) {
							int index = this.brokenBondToNeutralLossIndex.indexOf(Integer.valueOf(s));
							if (s == removedBondIndex) {
								if (index == -1) {
									this.brokenBondToNeutralLossIndex.add(Integer.valueOf(removedBondIndex));
									this.neutralLossIndex.add(Integer.valueOf(i));
								}
								continue;
							}
						}
						
					}
					return true;
				}
				else if (newGeneratedTopDownFragments[1].getAtomsArray().equals(this.neutralLosses.get(i).get(ii))) {

					/*
					 * check for previous broken bonds caused by neutral loss
					 */
					boolean[] brokenBondsArray = newGeneratedTopDownFragments[0].getBrokenBondsArray();
					
					for (int s = 0; s < brokenBondsArray.length; s++) {
						if(brokenBondsArray[s]) {
							int index = this.brokenBondToNeutralLossIndex.indexOf(Integer.valueOf(s));
							if (s == removedBondIndex) {
								if (index == -1) {
									this.brokenBondToNeutralLossIndex.add(Integer.valueOf(removedBondIndex));
									this.neutralLossIndex.add(Integer.valueOf(i));
								}
								continue;
							}
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
	 * @param childFragments
	 * @param precursorFragment
	 * @param toProcess
	 * @param ringBondArray
	 * @param lastCuttedRingBond
	 * @return
	 */
	private void createRingBondCleavedFragments(Collection<Fragment> childFragments,
			Fragment precursorFragment, Queue<Fragment> toProcess, boolean[] ringBondArray,
			Queue<Integer> lastCuttedRingBond) {
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
			for (int currentBond = nextRingBondToCut; currentBond < ringBondArray.length; currentBond++) {
				if (!ringBondArray[currentBond])
					continue;
				if (currentFragment.getBrokenBondsArray()[currentBond])
					continue;
				Fragment[] newFragments = currentFragment.fragment(currentBond);

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
						childFragments.add(newFragments[k]);
					}
				}

				if (newFragments.length == 1) {
					if (newFragments[0].getAddedToQueueCounts() < this.maximumNumberOfAFragmentAddedToQueue) {
						toProcess.add(newFragments[0]);
						lastCuttedRingBond.add(Integer.valueOf(currentBond));
					} else {
						childFragments.add(newFragments[0]);
					}
				}
			}
		}
	}

	/**
	 * Generate fragments by removing bonds that were skipped due to ring bond cleavage.
	 * 
	 * @throws Exception
	 */
	private void generateFragmentsOfSkippedBonds(Collection<Fragment> childFragments, Fragment precursorFragment) {
		int lastSkippedBonds = precursorFragment.getLastSkippedBond();
		int lastCuttedBond = getLastSetBit(precursorFragment.getBrokenBondsArray());
		if (lastSkippedBonds == -1)
			return;
		for (int currentBond = lastSkippedBonds; currentBond < lastCuttedBond; currentBond++) {

			if (this.ringBondArray[currentBond])
				continue;
			if (!precursorFragment.getBondsArray()[currentBond])
				continue;

			final Fragment[] newFragments = precursorFragment.fragment(currentBond);

			Fragmenter.processGeneratedFragments(newFragments);
			if (newFragments.length == 2) {
				this.checkForNeutralLossesAdaptMolecularFormulas(newFragments, currentBond);
			} else {
				System.err.println("problem generating fragments"); //$NON-NLS-1$
				System.exit(1);
			}

			for (int k = 0; k < newFragments.length; k++) {
				if (newFragments.length == 2) {
					childFragments.add(newFragments[k]);
				}
			}
		}
	}
	
	/**
	 * 
	 * @param atomIdx
	 * @return int
	 */
	private int getImplicitHydrogenCount(final int atomIdx) {
		return this.prec.getAtom(atomIdx).getImplicitHydrogenCount().intValue();
	}
	
	/**
	 * 
	 * @param array
	 * @return int
	 */
	private static int getLastSetBit(final boolean[] array) {
		for(int i = array.length - 1; i > -1; i--) {
			if(array[i]) {
				return i;
			}
		}
		
		return -1;
	}
}