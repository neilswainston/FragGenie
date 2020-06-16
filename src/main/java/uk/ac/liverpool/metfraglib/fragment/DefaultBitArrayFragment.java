package uk.ac.liverpool.metfraglib.fragment;

import de.ipbhalle.metfraglib.FastBitArray;
import de.ipbhalle.metfraglib.exceptions.AtomTypeNotKnownFromInputListException;
import de.ipbhalle.metfraglib.fragment.AbstractFragment;
import de.ipbhalle.metfraglib.interfaces.IMolecularFormula;
import de.ipbhalle.metfraglib.interfaces.IMolecularStructure;
import de.ipbhalle.metfraglib.molecularformula.BitArrayFragmentMolecularFormula;
import de.ipbhalle.metfraglib.parameter.Constants;
import de.ipbhalle.metfraglib.precursor.DefaultPrecursor;

/**
 * FastBitArrayFragment is an memory efficient way to store a fragment. It is
 * always related to a CDK AtomContainer object.
 * 
 * @author c-ruttkies
 *
 */
public abstract class DefaultBitArrayFragment extends AbstractFragment {

	protected short numberHydrogens;

	/**
	 * atoms represented as FastBitArray object
	 * 
	 */
	protected de.ipbhalle.metfraglib.FastBitArray atomsFastBitArray;
	protected de.ipbhalle.metfraglib.FastBitArray bondsFastBitArray;
	protected de.ipbhalle.metfraglib.FastBitArray brokenBondsFastBitArray;

	/**
	 * constructor setting precursor molecule of fragment
	 * 
	 * @param precursor
	 */
	public DefaultBitArrayFragment(IMolecularStructure precursorMolecule) {
		// super(precursorMolecule);
		this.atomsFastBitArray = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount(), true);
		this.bondsFastBitArray = new FastBitArray(precursorMolecule.getNonHydrogenBondCount(), true);
		this.brokenBondsFastBitArray = new FastBitArray(precursorMolecule.getNonHydrogenBondCount());
		this.treeDepth = 0;
	}

	/**
	 * 
	 * @param atomsFastBitArray
	 * @param bondsFastBitArray
	 */
	public DefaultBitArrayFragment(de.ipbhalle.metfraglib.FastBitArray atomsFastBitArray,
			de.ipbhalle.metfraglib.FastBitArray bondsFastBitArray,
			de.ipbhalle.metfraglib.FastBitArray brokenBondsFastBitArray) {
		// super(precursorMolecule);
		this.atomsFastBitArray = atomsFastBitArray;
		this.bondsFastBitArray = bondsFastBitArray;
		this.brokenBondsFastBitArray = brokenBondsFastBitArray;
		this.treeDepth = 0;
	}

	@Override
	public double getMonoisotopicMass(IMolecularStructure precursorMolecule) {
		// return this.molecularFormula.getMonoisotopicMass();
		double mass = 0.0;
		for (int i = 0; i < this.atomsFastBitArray.getSize(); i++) {
			if (this.atomsFastBitArray.get(i)) {
				mass += precursorMolecule.getMassOfAtom(i);
			}
		}
		return mass;
	}

	public byte matchToPeak(IMolecularStructure precursorMolecule, int precursorIonTypeIndex, boolean isPositive) {

		double[] ionisationTypeMassCorrection = new double[] {
				Constants.getIonisationTypeMassCorrection(precursorIonTypeIndex, isPositive),
				Constants.getIonisationTypeMassCorrection(0, isPositive) };

		for (int i = 0; i < ionisationTypeMassCorrection.length; i++) {
			double currentFragmentMass = this.getMonoisotopicMass(precursorMolecule) + ionisationTypeMassCorrection[i];

			System.out.println(currentFragmentMass);
		}

		return -1;
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

	public void setNumberHydrogens(int numberHydrogens) {
		this.numberHydrogens = (short) numberHydrogens;
	}

	@Override
	public IMolecularFormula getMolecularFormula(IMolecularStructure precursorMolecule) {
		try {
			BitArrayFragmentMolecularFormula form = new BitArrayFragmentMolecularFormula(
					(DefaultPrecursor) precursorMolecule, this.atomsFastBitArray);
			return form;
		} catch (AtomTypeNotKnownFromInputListException e) {
			e.printStackTrace();
		}
		return null;
	}

	public void setTreeDepth(byte treeDepth) {
		this.treeDepth = treeDepth;
	}

	@Override
	public byte getTreeDepth() {
		return this.treeDepth;
	}

	@Override
	public int[] getBrokenBondIndeces() {
		return this.brokenBondsFastBitArray.getSetIndeces();
	}
}