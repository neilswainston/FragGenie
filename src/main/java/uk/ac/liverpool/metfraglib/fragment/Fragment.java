package uk.ac.liverpool.metfraglib.fragment;

import java.util.LinkedList;
import java.util.List;

import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Element;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;

import de.ipbhalle.metfraglib.FastBitArray;
import de.ipbhalle.metfraglib.additionals.MoleculeFunctions;
import de.ipbhalle.metfraglib.exceptions.AtomTypeNotKnownFromInputListException;
import de.ipbhalle.metfraglib.interfaces.IMolecularStructure;
import uk.ac.liverpool.metfraglib.precursor.Precursor;

public class Fragment {

	private int addedToQueueCounts;
	private short lastSkippedBond;
	private FastBitArray atomsFastBitArray;
	private FastBitArray bondsFastBitArray;
	private FastBitArray brokenBondsFastBitArray;
	private byte treeDepth;
	private int numberHydrogens;
	private boolean isValidFragment;
	private boolean discardedForFragmentation;

	/**
	 * constructor setting all bits of atomsFastBitArray and bondsFastBitArray to
	 * true entire structure is represented
	 * 
	 * @param precursor
	 * @throws AtomTypeNotKnownFromInputListException
	 */
	public Fragment(Precursor precursor) {
		this.addedToQueueCounts = 0;
		this.lastSkippedBond = -1;
		this.atomsFastBitArray = new FastBitArray(precursor.getNonHydrogenAtomCount(), true);
		this.bondsFastBitArray = new FastBitArray(precursor.getNonHydrogenBondCount(), true);
		this.brokenBondsFastBitArray = new FastBitArray(precursor.getNonHydrogenBondCount());
		this.treeDepth = 0;
	}

	/**
	 * constructor setting bits of atomsFastBitArray and bondsFastBitArray by given
	 * ones
	 * 
	 * @param precursor
	 * @param atomsFastBitArray
	 * @param bondsFastBitArray
	 * @param brokenBondsFastBitArray
	 * @param numberHydrogens
	 * @throws Exception
	 */
	private Fragment(de.ipbhalle.metfraglib.FastBitArray atomsFastBitArray,
			de.ipbhalle.metfraglib.FastBitArray bondsFastBitArray,
			de.ipbhalle.metfraglib.FastBitArray brokenBondsFastBitArray, int numberHydrogens) throws Exception {
		this.atomsFastBitArray = atomsFastBitArray;
		this.bondsFastBitArray = bondsFastBitArray;
		this.brokenBondsFastBitArray = brokenBondsFastBitArray;
		this.treeDepth = 0;
		this.numberHydrogens = numberHydrogens;
		this.lastSkippedBond = -1;
	}

	/**
	 * main function of fragment generation traverse the given fragment and return
	 * two new fragments by removing bond with bondIndexToRemove
	 * 
	 * @param fragment
	 * @param bondNumber
	 * @param bondAtoms
	 * @return Fragment[]
	 * @throws Exception
	 */
	public Fragment[] traverseMolecule(Precursor precursorMolecule, short bondIndexToRemove,
			short[] indecesOfBondConnectedAtoms) throws Exception {

		/*
		 * generate first fragment
		 */
		de.ipbhalle.metfraglib.FastBitArray atomArrayOfNewFragment_1 = new de.ipbhalle.metfraglib.FastBitArray(
				precursorMolecule.getNonHydrogenAtomCount());
		de.ipbhalle.metfraglib.FastBitArray bondArrayOfNewFragment_1 = new de.ipbhalle.metfraglib.FastBitArray(
				precursorMolecule.getNonHydrogenBondCount());
		de.ipbhalle.metfraglib.FastBitArray brokenBondArrayOfNewFragment_1 = this.getBrokenBondsFastBitArray().clone();
		int[] numberHydrogensOfNewFragment = new int[1];

		/*
		 * traverse to first direction from atomIndex connected by broken bond
		 */
		boolean stillOneFragment = this.traverseSingleDirection(precursorMolecule, indecesOfBondConnectedAtoms[0],
				indecesOfBondConnectedAtoms[1], bondIndexToRemove, atomArrayOfNewFragment_1, bondArrayOfNewFragment_1,
				brokenBondArrayOfNewFragment_1, numberHydrogensOfNewFragment);

		Fragment firstNewGeneratedFragment = new Fragment(atomArrayOfNewFragment_1, bondArrayOfNewFragment_1,
				brokenBondArrayOfNewFragment_1, numberHydrogensOfNewFragment[0]);

		/*
		 * only one fragment is generated when a ring bond was broken
		 */
		if (stillOneFragment) {
			firstNewGeneratedFragment.setTreeDepth(this.treeDepth);
			firstNewGeneratedFragment.setAddedToQueueCounts((byte) (this.getAddedToQueueCounts() + 1));
			Fragment[] newFrags = { firstNewGeneratedFragment };
			return newFrags;
		}
		/*
		 * generate second fragment
		 */
		de.ipbhalle.metfraglib.FastBitArray atomArrayOfNewFragment_2 = new de.ipbhalle.metfraglib.FastBitArray(
				precursorMolecule.getNonHydrogenAtomCount());
		de.ipbhalle.metfraglib.FastBitArray bondArrayOfNewFragment_2 = new de.ipbhalle.metfraglib.FastBitArray(
				precursorMolecule.getNonHydrogenBondCount());
		de.ipbhalle.metfraglib.FastBitArray brokenBondArrayOfNewFragment_2 = this.getBrokenBondsFastBitArray().clone();
		numberHydrogensOfNewFragment[0] = 0;

		/*
		 * traverse the second direction from atomIndex connected by broken bond
		 */
		this.traverseSingleDirection(precursorMolecule, indecesOfBondConnectedAtoms[1], indecesOfBondConnectedAtoms[0],
				bondIndexToRemove, atomArrayOfNewFragment_2, bondArrayOfNewFragment_2, brokenBondArrayOfNewFragment_2,
				numberHydrogensOfNewFragment);

		Fragment secondNewGeneratedFragment = new Fragment(atomArrayOfNewFragment_2, bondArrayOfNewFragment_2,
				brokenBondArrayOfNewFragment_2, numberHydrogensOfNewFragment[0]);

		firstNewGeneratedFragment.setTreeDepth((byte) (this.treeDepth + 1));

		secondNewGeneratedFragment.setTreeDepth((byte) (this.treeDepth + 1));

		Fragment[] newFrags = { firstNewGeneratedFragment, secondNewGeneratedFragment };

		return newFrags;

	}

	public int getAddedToQueueCounts() {
		return this.addedToQueueCounts;
	}

	public void setAddedToQueueCounts(byte addedToQueueCounts) {
		this.addedToQueueCounts = addedToQueueCounts;
	}

	public de.ipbhalle.metfraglib.FastBitArray getAtomsFastBitArray() {
		return this.atomsFastBitArray;
	}

	public de.ipbhalle.metfraglib.FastBitArray getBondsFastBitArray() {
		return this.bondsFastBitArray;
	}

	public de.ipbhalle.metfraglib.FastBitArray getBrokenBondsFastBitArray() {
		return this.brokenBondsFastBitArray;
	}

	public int getMaximalIndexOfRemovedBond() {
		return this.brokenBondsFastBitArray.getLastSetBit();
	}

	public short getLastSkippedBond() {
		return this.lastSkippedBond;
	}

	public void setLastSkippedBond(short lastSkippedBond) {
		this.lastSkippedBond = lastSkippedBond;
	}

	public int[] getBrokenBondIndeces() {
		return this.brokenBondsFastBitArray.getSetIndeces();
	}

	public double getMonoisotopicMass(Precursor precursorMolecule) {
		double mass = 0.0;
		for (int i = 0; i < this.atomsFastBitArray.getSize(); i++) {
			if (this.atomsFastBitArray.get(i)) {
				mass += precursorMolecule.getMassOfAtom(i);
			}
		}
		return mass;
	}

	public boolean isValidFragment() {
		return this.isValidFragment;
	}

	public void setAsValidFragment() {
		this.isValidFragment = true;
	}

	public boolean isDiscardedForFragmentation() {
		return this.discardedForFragmentation;
	}

	public void setAsDiscardedForFragmentation() {
		this.discardedForFragmentation = true;
	}

	@Override
	public Object clone() {
		try {
			Fragment clone = new Fragment(this.atomsFastBitArray.clone(), this.bondsFastBitArray.clone(),
					this.brokenBondsFastBitArray.clone(), this.numberHydrogens);
			clone.setTreeDepth(this.treeDepth);
			return clone;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	private void setTreeDepth(byte treeDepth) {
		this.treeDepth = treeDepth;
	}

	/**
	 * traverse the fragment to one direction starting from startAtomIndex to set
	 * FastBitArrays of new fragment
	 * 
	 * @param startAtomIndex
	 * @param fragment
	 * @param bondIndexToRemove
	 * @param atomArrayOfNewFragment
	 * @param bondArrayOfNewFragment
	 * @param brokenBondArrayOfNewFragment
	 * @param numberHydrogensOfNewFragment
	 * @return
	 */
	private boolean traverseSingleDirection(Precursor precursorMolecule, short startAtomIndex, int endAtomIndex,
			int bondIndexToRemove, de.ipbhalle.metfraglib.FastBitArray atomArrayOfNewFragment,
			de.ipbhalle.metfraglib.FastBitArray bondArrayOfNewFragment,
			de.ipbhalle.metfraglib.FastBitArray brokenBondArrayOfNewFragment, int[] numberHydrogensOfNewFragment) {
		de.ipbhalle.metfraglib.FastBitArray bondFastBitArrayOfCurrentFragment = this.getBondsFastBitArray();
		/*
		 * when traversing the fragment graph then we want to know if we already visited
		 * a node (atom) need to be done for checking of ringed structures if traversed
		 * an already visited atom, then no new fragment was generated
		 */
		de.ipbhalle.metfraglib.FastBitArray visited = new de.ipbhalle.metfraglib.FastBitArray(
				precursorMolecule.getNonHydrogenAtomCount());
		numberHydrogensOfNewFragment[0] = 0;

		/*
		 * traverse molecule in the first direction
		 */
		java.util.Stack<short[]> toProcessConnectedAtoms = new java.util.Stack<>();
		java.util.Stack<Short> toProcessAtom = new java.util.Stack<>();
		toProcessConnectedAtoms.push(precursorMolecule.getConnectedAtomIndecesOfAtomIndex(startAtomIndex));
		toProcessAtom.push(startAtomIndex);
		visited.set(startAtomIndex);
		boolean stillOneFragment = false;
		/*
		 * set the first atom of possible new fragment atom is of the one direction of
		 * cutted bond
		 */
		atomArrayOfNewFragment.set(startAtomIndex);
		numberHydrogensOfNewFragment[0] += precursorMolecule.getNumberHydrogensConnectedToAtomIndex(startAtomIndex);
		while (!toProcessConnectedAtoms.isEmpty()) {
			short[] nextAtoms = toProcessConnectedAtoms.pop();
			short midAtom = toProcessAtom.pop();
			for (int i = 0; i < nextAtoms.length; i++) {
				/*
				 * did we visit the current atom already?
				 */
				int currentBondNumber = precursorMolecule.getBondIndexFromAtomAdjacencyList(nextAtoms[i], midAtom) - 1;

				if (!bondFastBitArrayOfCurrentFragment.get(currentBondNumber)
						|| currentBondNumber == bondIndexToRemove) {
					continue;
				}
				/*
				 * if we visited the current atom already then we do not have to check it again
				 */
				if (visited.get(nextAtoms[i])) {
					bondArrayOfNewFragment.set(currentBondNumber);
					continue;
				}
				/*
				 * if we reach the second atom of the cleaved bond then still one fragment is
				 * present
				 */
				if (nextAtoms[i] == endAtomIndex) {
					stillOneFragment = true;
				}

				visited.set(nextAtoms[i]);
				atomArrayOfNewFragment.set(nextAtoms[i]);
				/*
				 * add number of hydrogens of current atom
				 */
				numberHydrogensOfNewFragment[0] += precursorMolecule
						.getNumberHydrogensConnectedToAtomIndex(nextAtoms[i]);
				bondArrayOfNewFragment
						.set(precursorMolecule.getBondIndexFromAtomAdjacencyList(midAtom, nextAtoms[i]) - 1);
				toProcessConnectedAtoms.push(precursorMolecule.getConnectedAtomIndecesOfAtomIndex(nextAtoms[i]));
				toProcessAtom.push(nextAtoms[i]);
			}

		}

		brokenBondArrayOfNewFragment.set(bondIndexToRemove);
		bondArrayOfNewFragment.set(bondIndexToRemove, false);

		return stillOneFragment;
	}

	public IAtomContainer getStructureAsIAtomContainer(IMolecularStructure precursorMolecule) {
		return this.getStructureAsIAtomContainer(precursorMolecule, null, null);
	}

	private IAtomContainer getStructureAsIAtomContainer(IMolecularStructure precursorMolecule, IAtom atomToAdd,
			Integer atomPositionToAdd) {
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		if (this.atomsFastBitArray.cardinality() == 1) {
			IAtom curAtom = precursorMolecule.getStructureAsIAtomContainer()
					.getAtom(this.atomsFastBitArray.getFirstSetBit());
			if (atomToAdd != null
					&& precursorMolecule.getStructureAsIAtomContainer().indexOf(curAtom) == atomPositionToAdd) {
				IBond newBond = new Bond(curAtom, atomToAdd);
				fragmentStructure.addAtom(atomToAdd);
				fragmentStructure.addBond(newBond);
			}
			fragmentStructure.addAtom(curAtom);
			return fragmentStructure;
		}
		for (int i = 0; i < this.bondsFastBitArray.getSize(); i++) {
			if (this.bondsFastBitArray.get(i)) {
				IBond curBond = precursorMolecule.getStructureAsIAtomContainer().getBond(i);
				for (IAtom atom : curBond.atoms()) {
					if (atomToAdd != null
							&& precursorMolecule.getStructureAsIAtomContainer().indexOf(atom) == atomPositionToAdd) {
						IBond newBond = new Bond(atom, atomToAdd);
						fragmentStructure.addAtom(atomToAdd);
						fragmentStructure.addBond(newBond);
					}
					fragmentStructure.addAtom(atom);
				}
				fragmentStructure.addBond(curBond);
			}
		}
		// loss of hydrogens
		// MoleculeFunctions.prepareAtomContainer(fragmentStructure);

		List<IAtom> atomsToAdd = new LinkedList<>();
		List<IBond> bondsToAdd = new LinkedList<>();

		for (int i = 0; i < fragmentStructure.getAtomCount(); i++) {
			int bondOrderSum = 0;
			java.util.List<IBond> bonds = fragmentStructure.getConnectedBondsList(fragmentStructure.getAtom(i));
			for (int ii = 0; ii < bonds.size(); ii++) {
				bondOrderSum += bonds.get(ii).getOrder().numeric();
			}
			if (fragmentStructure.getAtom(i).getBondOrderSum() != null
					&& fragmentStructure.getAtom(i).getBondOrderSum() > bondOrderSum) {
				for (int k = 0; k < (fragmentStructure.getAtom(i).getBondOrderSum() - bondOrderSum); k++) {
					org.openscience.cdk.Atom hydrogenAtom = new org.openscience.cdk.Atom(new Element("H")); //$NON-NLS-1$
					IBond bond = new Bond(hydrogenAtom, fragmentStructure.getAtom(i));
					atomsToAdd.add(hydrogenAtom);
					bondsToAdd.add(bond);
				}
				fragmentStructure.getAtom(i).setImplicitHydrogenCount(0);
			}
		}
		for (IBond bond : bondsToAdd)
			fragmentStructure.addBond(bond);
		for (IAtom atom : atomsToAdd)
			fragmentStructure.addAtom(atom);
		atomsToAdd = new LinkedList<>();
		bondsToAdd = new LinkedList<>();
		for (int i = 0; i < fragmentStructure.getAtomCount(); i++) {
			if (fragmentStructure.getAtom(i).getImplicitHydrogenCount() != null
					&& fragmentStructure.getAtom(i).getImplicitHydrogenCount() > 0) {
				for (int k = 0; k < fragmentStructure.getAtom(i).getImplicitHydrogenCount(); k++) {
					org.openscience.cdk.Atom hydrogenAtom = new org.openscience.cdk.Atom(new Element("H")); //$NON-NLS-1$
					IBond bond = new Bond(hydrogenAtom, fragmentStructure.getAtom(i));
					atomsToAdd.add(hydrogenAtom);
					bondsToAdd.add(bond);
				}
				fragmentStructure.getAtom(i).setImplicitHydrogenCount(0);
			}
		}
		for (IBond bond : bondsToAdd)
			fragmentStructure.addBond(bond);
		for (IAtom atom : atomsToAdd)
			fragmentStructure.addAtom(atom);
		// MoleculeFunctions.removeHydrogens(fragmentStructure);

		fragmentStructure = MoleculeFunctions.convertExplicitToImplicitHydrogens(fragmentStructure);
		return fragmentStructure;
	}
}