package uk.ac.liverpool.metfrag;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * 
 * @author neilswainston
 */
public class Fragment implements Comparable<Fragment> {
	
	/**
	 * 
	 */
	private final Precursor prec;

	/**
	 * 
	 */
	private boolean[] atomsArray;

	/**
	 * 
	 */
	private boolean[] bondsArray;

	/**
	 * 
	 */
	private boolean[] brokenBondsArray;

	/**
	 * 
	 */
	private float monoisotopicMass = Float.NaN;

	/**
	 * Constructor.
	 * 
	 * @param precursor
	 */
	Fragment(final Precursor precursor) {
		this(precursor,
				new boolean[precursor.getAtomCount()],
				new boolean[precursor.getBondCount()]);
		
		Arrays.fill(this.atomsArray, true);
	}

	/**
	 * Constructor.
	 * 
	 * @param precursor
	 * @param atoms
	 * @param brokenBonds
	 */
	private Fragment(final Precursor precursor, final boolean[] atoms, final boolean[] brokenBonds) {
		this.prec = precursor;
		this.atomsArray = atoms;
		this.brokenBondsArray = brokenBonds;
	}
	
	/**
	 * 
	 * @return Collection<List<Object>>
	 */
	public Collection<Fragment> fragment(final int maxBrokenBonds, final float minMass) {
		final Collection<Fragment> fragments = new TreeSet<>();
		fragment(this, fragments, maxBrokenBonds, minMass);
		return fragments;
	}
	
	@Override
	public String toString() {
		return getFormula();
	}
	
	/**
	 * 
	 * @return float
	 */
	public float getMonoisotopicMass() {
		if(Float.isNaN(this.monoisotopicMass)) {
			float mass = 0.0f;

			for (int i = 0; i < this.atomsArray.length; i++) {
				if (this.atomsArray[i]) {
					mass += this.prec.getAtomMass(i);
				}
			}
			
			this.monoisotopicMass = mass;
		}
		
		return this.monoisotopicMass;
	}

	/**
	 * 
	 * @return String
	 */
	public String getFormula() {
		final String HYDROGEN = "H"; //$NON-NLS-1$
		final Map<String, Integer> elementCount = new TreeMap<>();

		for (int atomIdx = 0; atomIdx < this.atomsArray.length; atomIdx++) {
			if (this.atomsArray[atomIdx]) {
				final String element = this.prec.getAtomSymbol(atomIdx);

				if (elementCount.get(element) == null) {
					elementCount.put(element, Integer.valueOf(1));
				} else {
					elementCount.put(element, Integer.valueOf(elementCount.get(element).intValue() + 1));
				}

				final int hCount = this.prec.getAtomImplicitHydrogenCount(atomIdx);

				if (elementCount.get(HYDROGEN) == null) {
					if(hCount > 0) {
						elementCount.put(HYDROGEN, Integer.valueOf(hCount));
					}
				} else {
					elementCount.put(HYDROGEN, Integer.valueOf(elementCount.get(HYDROGEN).intValue() + hCount));
				}
			}
		}

		final StringBuilder builder = new StringBuilder();

		for (Map.Entry<String, Integer> entry : elementCount.entrySet()) {
			builder.append(entry.getKey());
			final int count = entry.getValue().intValue();

			if (count > 1) {
				builder.append(count);
			}
		}

		return builder.toString();
	}
	
	/**
	 * 
	 * @return Collection<Bond>
	 */
	public Collection<Bond> getBrokenBonds() {
		final List<Bond> bonds = new ArrayList<>();
		
		for (int i = 0; i < this.brokenBondsArray.length; i++) {
			if(this.brokenBondsArray[i]) {
				bonds.add(this.prec.getBond(i));
			}
		}
		
		Collections.sort(bonds);
		
		return bonds;
	}
	
	/**
	 * 
	 * @return int
	 */
	int getNumBrokenBonds() {
		int numBrokenBonds = 0;
		
		for (int i = 0; i < this.brokenBondsArray.length; i++) {
			if(this.brokenBondsArray[i]) {
				numBrokenBonds++;
			}
		}
		
		return numBrokenBonds;
	}
	
	/**
	 * 
	 * @param fragment
	 * @param fragments
	 * @param maxBrokenBonds
	 */
	private static void fragment(final Fragment fragment, final Collection<Fragment> fragments, final int maxBrokenBonds, final float minMass) {
		if(fragment.getNumBrokenBonds() <= maxBrokenBonds && fragment.getMonoisotopicMass() > minMass) {
			fragments.add(fragment);
			
			if(fragment.getNumBrokenBonds() < maxBrokenBonds) {
				for(int bondIdx = 0; bondIdx < fragment.brokenBondsArray.length; bondIdx++) {
					if(!fragment.brokenBondsArray[bondIdx]) {
						for(final Fragment childFragment : fragment.fragmentBond(bondIdx)) {
							fragment(childFragment, fragments, maxBrokenBonds, minMass);
						}
					}
				}
			}
		}
	}

	/**
	 * 
	 * @param bondIdx
	 * @return Fragment[]
	 */
	Fragment[] fragmentBond(final int bondIdx) {
		final int[] bondConnectedAtoms = this.prec.getBondConnectedAtoms(bondIdx);
		
		// Generate first fragment:
		// Traverse to first direction from atomIndex connected by broken bond:
		final Object[] result1 = this.traverse(bondConnectedAtoms[0], bondConnectedAtoms[1], bondIdx);

		final Fragment fragment1 = (Fragment)result1[1];

		// Only one fragment is generated when a ring bond was broken:
		if (((Boolean)result1[0]).booleanValue()) {
			return new Fragment[] { fragment1 };
		}
		
		// Generate second fragment:
		// Traverse the second direction from atomIndex connected by broken bond:
		final Object[] result2 = this.traverse(bondConnectedAtoms[1], bondConnectedAtoms[0], bondIdx);

		final Fragment fragment2 = (Fragment)result2[1];
		
		return new Fragment[] { fragment1, fragment2 };
	}
	


	/**
	 * Traverse the fragment to one direction starting from startAtom.
	 * 
	 * @param startAtomIdx
	 * @param endAtomIdx
	 * @param removeBondIdx
	 * @return Object[]
	 */
	private Object[] traverse(final int startAtomIdx, final int endAtomIdx, final int removeBondIdx) {
		final boolean[] newAtomArray = new boolean[this.prec.getAtomCount()];
		final boolean[] currentBrokenBondArray = this.brokenBondsArray;
		final boolean[] newBrokenBondArray = this.brokenBondsArray.clone();
		
		// When traversing the fragment graph, we want to know if we already visited
		// a node (atom) to check for ringed structures.
		// If traversed an already visited atom, then no new fragment was generated.
		final boolean[] visited = new boolean[this.prec.getAtomCount()];

		// Traverse molecule in the first direction:
		final Stack<int[]> toProcessConnectedAtoms = new Stack<>();
		final Stack<Integer> toProcessAtom = new Stack<>();
		
		toProcessConnectedAtoms.push(this.prec.getConnectedAtomIdxs(startAtomIdx));
		toProcessAtom.push(Integer.valueOf(startAtomIdx));
		visited[startAtomIdx] = true;
		
		boolean singleFragment = false;
		
		// Set the first atom of possible new fragment atom as the one direction of cut bond.
		newAtomArray[startAtomIdx] = true;
		
		while (!toProcessConnectedAtoms.isEmpty()) {
			final int midAtomIdx = toProcessAtom.pop().intValue();
			
			for (int nextAtomIdx : toProcessConnectedAtoms.pop()) {
				// Did we visit the current atom already?
				final int currentBondIdx = this.prec.getBondIdx(midAtomIdx, nextAtomIdx);

				if (currentBrokenBondArray[currentBondIdx] || currentBondIdx == removeBondIdx) {
					continue;
				}
				
				// If we visited the current atom already, then we do not have to check it again:
				if (visited[nextAtomIdx]) {
					continue;
				}
				
				// If we reach the second atom of the cleaved bond then still one fragment is present:
				if (nextAtomIdx == endAtomIdx) {
					singleFragment = true;
				}

				visited[nextAtomIdx] = true;
				newAtomArray[nextAtomIdx] = true;

				toProcessConnectedAtoms.push(this.prec.getConnectedAtomIdxs(nextAtomIdx));
				toProcessAtom.push(Integer.valueOf(nextAtomIdx));
			}
		}

		newBrokenBondArray[removeBondIdx] = true;
		
		final Fragment newFragment = new Fragment(this.prec, newAtomArray, newBrokenBondArray);

		return new Object[] {Boolean.valueOf(singleFragment), newFragment};
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((this.atomsArray == null) ? 0 : this.atomsArray.hashCode());
		result = prime * result + ((this.bondsArray == null) ? 0 : this.bondsArray.hashCode());
		result = prime * result + ((this.brokenBondsArray == null) ? 0 : this.brokenBondsArray.hashCode());
		result = prime * result + ((this.prec == null) ? 0 : this.prec.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
			
		if (obj == null) {
			return false;
		}
			
		if (getClass() != obj.getClass()) {
			return false;
		}
			
		
		final Fragment other = (Fragment) obj;
		
		if (this.atomsArray == null) {
			if (other.atomsArray != null) {
				return false;
			}	
		}
		else if (!this.atomsArray.equals(other.atomsArray)) {
			return false;
		}
			
		
		if (this.bondsArray == null) {
			if (other.bondsArray != null) {
				return false;
			}	
		}
		else if (!this.bondsArray.equals(other.bondsArray)) {
			return false;
		}
			
		if (this.brokenBondsArray == null) {
			if (other.brokenBondsArray != null) {
				return false;
			}
		}
		else if (!this.brokenBondsArray.equals(other.brokenBondsArray)) {
			return false;
		}
			
		if (this.prec == null) {
			if (other.prec != null) {
				return false;
			}
				
		}
		else if (!this.prec.equals(other.prec)) {
			return false;
		}
			
		return true;
	}

	@Override
	public int compareTo(final Fragment obj) {
		return (int)((this.getMonoisotopicMass() - obj.getMonoisotopicMass()) * 1000000);
	}
}