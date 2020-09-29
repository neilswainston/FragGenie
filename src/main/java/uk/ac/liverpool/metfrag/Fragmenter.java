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

/**
 * 
 * @author neilswainston
 */
public class Fragmenter {
	
	/**
	 * 
	 */
	private IAtomContainer prec;
	
	/**
	 * 
	 * @param precursor
	 * @throws CDKException
	 * @throws IOException
	 */
	public Fragmenter(final IAtomContainer precursor) throws IOException, CDKException {
		this.prec = precursor;
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
	private static Collection<Fragment> getFragments(final Fragment parentFragment) {
		final Collection<Fragment> childFragments = new ArrayList<>();
		final boolean[] bondsArray = parentFragment.getBondsArray();
		final boolean[] ringBonds = new boolean[bondsArray.length];
		final Queue<Fragment> ringBondFragments = new LinkedList<>();
		final Queue<Integer> lastBrokenBonds = new LinkedList<>();
			
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
				ringBondFragments.add(currentChildFragments[0]);
				lastBrokenBonds.add(Integer.valueOf(bondIdx));
			}
			else if (currentChildFragments.length == 2) {
				// If two new fragments have been generated set them as valid:
				childFragments.addAll(Arrays.asList(currentChildFragments));
			}
		}
		
		// Create fragments by ring bond cleavage:
		childFragments.addAll(createRingBondCleavedFragments(ringBondFragments, ringBonds, lastBrokenBonds));

		return childFragments;
	}

	/**
	 * 
	 * @param fragments
	 * @param ringBondArray
	 * @param lastBrokenBonds
	 * @return Collection<Fragment>
	 */
	private static Collection<Fragment> createRingBondCleavedFragments(final Queue<Fragment> fragments, final boolean[] ringBonds, final Queue<Integer> lastBrokenBonds) {
		final Collection<Fragment> childFragments = new ArrayList<>();
		
		// Process all fragments that have been cut in a ring without generating a new one:
		while (!fragments.isEmpty() && lastBrokenBonds.size() != 0) {
			final Fragment fragment = fragments.poll();
			final int nextRingBondToCut = lastBrokenBonds.poll().intValue() + 1;

			for (int currentBond = nextRingBondToCut; currentBond < ringBonds.length; currentBond++) {
				if (ringBonds[currentBond] && !fragment.getBrokenBondsArray()[currentBond]) {
					final Fragment[] newFragments = fragment.fragment(currentBond);
					
					childFragments.addAll(Arrays.asList(newFragments));
					
					if (newFragments.length == 1) {
						lastBrokenBonds.add(Integer.valueOf(currentBond));
					}
				}
			}
		}
		
		return childFragments;
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