package uk.ac.liverpool.metfrag;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.Queue;
import java.util.TreeSet;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

public class Fragmenter {
	
	/**
	 * 
	 */
	private boolean ringBondsInitialised = false;
	
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
	 * @param precursor
	 * @throws CDKException
	 * @throws IOException
	 */
	public Fragmenter(final IAtomContainer precursor) throws IOException, CDKException {
		this.prec = precursor;
		this.ringBondArray = new boolean[this.prec.getBondCount()];
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
				final Collection<Fragment> childFragments = getFragments(fragment);
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
	private Collection<Fragment> getFragments(final Fragment parentFragment) {
		final Collection<Fragment> childFragments = new ArrayList<>();
		final boolean[] bondsArray = parentFragment.getBondsArray();
		
		final boolean[] ringBonds = new boolean[bondsArray.length];
		final Queue<Fragment> ringBondCuttedFragments = new LinkedList<>();
		final Queue<Integer> lastCuttedBondOfRing = new LinkedList<>();
		
		// Generate fragments of skipped bonds:
		if (this.ringBondsInitialised) {
			this.generateFragmentsOfSkippedBonds(childFragments, parentFragment);
		}
			
		// Get the last bond index that was removed; from there on the next bonds will be removed:
		final int nextBondBreakIdx = getLastSetBit(parentFragment.getBrokenBondsArray()) + 1;
		
		// Start from the last broken bond index:
		for (int bondIdx = nextBondBreakIdx; bondIdx < bondsArray.length; bondIdx++) {
			if (!bondsArray[bondIdx]) {
				continue;
			}
			
			// Try to generate at most two fragments by the breaking of the given bond:
			final Fragment[] currentChildFragments = parentFragment.fragment(bondIdx);
			
			/*
			 * In case the precursor wasn't split, try to cleave an additional bond until:
			 * 
			 * 1. two fragments are generated; or
			 * 2. the maximum number of trials have been reached; or
			 * 3. no further bond can be removed.
			 */
			if (currentChildFragments.length == 1) {
				ringBonds[bondIdx] = true;
				currentChildFragments[0].setLastSkippedBond(bondIdx + 1);
				ringBondCuttedFragments.add(currentChildFragments[0]);
				lastCuttedBondOfRing.add(Integer.valueOf(bondIdx));
				
				if (!this.ringBondsInitialised) {
					this.ringBondArray[bondIdx] = true;
				}
			}
			else if (currentChildFragments.length == 2) {
				// If two new fragments have been generated set them as valid:
				childFragments.addAll(Arrays.asList(currentChildFragments));
			}
		}
		/*
		 * create fragments by ring bond cleavage and store them in the given vector
		 */
		createRingBondCleavedFragments(childFragments, parentFragment, ringBondCuttedFragments, ringBonds, lastCuttedBondOfRing);
		this.ringBondsInitialised = true;

		return childFragments;
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
	private static void createRingBondCleavedFragments(Collection<Fragment> childFragments,
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
					toProcess.add(newFragments[0]);
					lastCuttedRingBond.add(Integer.valueOf(currentBond));
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

			for (int k = 0; k < newFragments.length; k++) {
				if (newFragments.length == 2) {
					childFragments.add(newFragments[k]);
				}
			}
		}
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