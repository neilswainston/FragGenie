package uk.ac.liverpool.metfraglib.fragmenter;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;

import de.ipbhalle.metfraglib.FastBitArray;
import de.ipbhalle.metfraglib.additionals.NeutralLosses;
import de.ipbhalle.metfraglib.fragment.BitArrayNeutralLoss;
import de.ipbhalle.metfraglib.list.FragmentList;
import de.ipbhalle.metfraglib.parameter.Constants;
import uk.ac.liverpool.metfraglib.precursor.Precursor;
import uk.ac.liverpool.metfraglib.candidate.PrecursorCandidate;
import uk.ac.liverpool.metfraglib.fragment.Fragment;

public class Fragmenter {

	private Double minimumMassDeviationForFragmentGeneration = Constants.DEFAULT_MIN_MASS_DEV_FOR_FRAGMENT_GENERATION;
	private Byte maximumNumberOfAFragmentAddedToQueue;
	private boolean ringBondsInitialised;
	private FastBitArray ringBondFastBitArray;
	private PrecursorCandidate candidate;
	private Double minimumFragmentMassLimit;
	private BitArrayNeutralLoss[] detectedNeutralLosses;
	private List<Short> brokenBondToNeutralLossIndex = new ArrayList<>();
	private List<Integer> neutralLossIndex = new ArrayList<>();

	/**
	 * 
	 * @param candidate
	 * @param maximumTreeDepth
	 * @throws Exception
	 */
	public Fragmenter(final PrecursorCandidate candidate) throws Exception {
		this.candidate = candidate;
		this.minimumFragmentMassLimit = 0.0;
		this.maximumNumberOfAFragmentAddedToQueue = 2;
		this.ringBondsInitialised = false;
		this.ringBondFastBitArray = new FastBitArray(this.candidate.getPrecursorMolecule().getNonHydrogenBondCount(),
				false);
		this.detectedNeutralLosses = new NeutralLosses().getMatchingAtoms(this.candidate.getPrecursorMolecule());
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
					newGeneratedTopDownFragments[1].getMolecularFormula(this.candidate.getPrecursorMolecule())
							.setNumberHydrogens((short) (newGeneratedTopDownFragments[1]
									.getMolecularFormula(this.candidate.getPrecursorMolecule()).getNumberHydrogens()
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
							newGeneratedTopDownFragments[1].getMolecularFormula(this.candidate.getPrecursorMolecule())
									.setNumberHydrogens((short) (newGeneratedTopDownFragments[1]
											.getMolecularFormula(this.candidate.getPrecursorMolecule())
											.getNumberHydrogens()
											+ this.detectedNeutralLosses[this.neutralLossIndex.get(index)]
													.getHydrogenDifference()));
						}
					}
					return true;
				} else if (newGeneratedTopDownFragments[1].getAtomsFastBitArray()
						.equals(this.detectedNeutralLosses[i].getNeutralLossAtomFastBitArray(ii))) {
					newGeneratedTopDownFragments[0].getMolecularFormula(this.candidate.getPrecursorMolecule())
							.setNumberHydrogens((short) (newGeneratedTopDownFragments[0]
									.getMolecularFormula(this.candidate.getPrecursorMolecule()).getNumberHydrogens()
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
							newGeneratedTopDownFragments[0].getMolecularFormula(this.candidate.getPrecursorMolecule())
									.setNumberHydrogens((short) (newGeneratedTopDownFragments[0]
											.getMolecularFormula(this.candidate.getPrecursorMolecule())
											.getNumberHydrogens()
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
			short[] indecesOfBondConnectedAtoms = this.candidate.getPrecursorMolecule()
					.getConnectedAtomIndecesOfBondIndex(i);
			/*
			 * try to generate at most two fragments by the removal of the given bond
			 */
			Fragment[] newGeneratedTopDownFragments = precursorFragment
					.traverseMolecule(this.candidate.getPrecursorMolecule(), i, indecesOfBondConnectedAtoms);
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

			short[] connectedAtomIndeces = this.candidate.getPrecursorMolecule()
					.getConnectedAtomIndecesOfBondIndex(currentBond);

			Fragment[] newFragments = precursorFragment.traverseMolecule(this.candidate.getPrecursorMolecule(),
					currentBond, connectedAtomIndeces);

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
			if (newGeneratedTopDownFragments[0]
					.getMonoisotopicMass(this.candidate.getPrecursorMolecule()) <= this.minimumFragmentMassLimit
							- this.minimumMassDeviationForFragmentGeneration) {
				newGeneratedTopDownFragments[0].setAsDiscardedForFragmentation();
			}
			if (newGeneratedTopDownFragments[1]
					.getMonoisotopicMass(this.candidate.getPrecursorMolecule()) <= this.minimumFragmentMassLimit
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
				short[] connectedAtomIndeces = this.candidate.getPrecursorMolecule()
						.getConnectedAtomIndecesOfBondIndex(currentBond);
				newFragments = currentFragment.traverseMolecule(this.candidate.getPrecursorMolecule(), currentBond,
						connectedAtomIndeces);

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
}
