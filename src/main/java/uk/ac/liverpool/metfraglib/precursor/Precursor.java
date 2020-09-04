package uk.ac.liverpool.metfraglib.precursor;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;

import de.ipbhalle.metfraglib.FastBitArray;
import de.ipbhalle.metfraglib.additionals.MoleculeFunctions;
import de.ipbhalle.metfraglib.parameter.Constants;

public class Precursor {

	/**
	 * 
	 * @param smiles
	 * @return Precursor
	 * @throws Exception
	 */
	public static Precursor fromSmiles(final String smiles) throws Exception {
		final IAtomContainer molecule = MoleculeFunctions.getAtomContainerFromSMILES(smiles);
		MoleculeFunctions.prepareAtomContainer(molecule, true);
		MoleculeFunctions.convertExplicitToImplicitHydrogens(molecule);
		return new Precursor(molecule);
	}

	/**
	 * 
	 */
	private final FastBitArray aromaticBonds;

	/**
	 * 
	 */
	private final short[] atomAdjacencyList;

	/**
	 * 
	 */
	private final List<short[]> atomIndexToConnectedAtomIndeces = new ArrayList<>();

	/**
	 * 
	 */
	private final short[][] bondIndexToConnectedAtomIndeces;
	/**
	 * 
	 */
	private final double[] massesOfAtoms;

	/**
	 * 
	 */
	private final Integer[] numberHydrogensConnectedToAtom;

	/**
	 * 
	 */
	private final IAtomContainer precursorMolecule;

	/**
	 * 
	 */
	private final FastBitArray[] ringBondToBelongingRingBondIndeces;

	/**
	 * 
	 * @param precursorMolecule
	 * @throws Exception
	 */
	private Precursor(final IAtomContainer precursorMolecule) throws Exception {
		this.precursorMolecule = precursorMolecule;
		this.bondIndexToConnectedAtomIndeces = new short[this.getNonHydrogenBondCount()][2];
		this.ringBondToBelongingRingBondIndeces = new FastBitArray[this.precursorMolecule.getBondCount() + 1];
		this.aromaticBonds = new FastBitArray(this.getNonHydrogenBondCount());
		this.atomAdjacencyList = new short[getIndex(this.getNonHydrogenAtomCount() - 2,
				this.getNonHydrogenAtomCount() - 1) + 1];
		this.numberHydrogensConnectedToAtom = new Integer[this.getNonHydrogenAtomCount()];
		this.massesOfAtoms = new double[this.getNonHydrogenAtomCount()];
		this.initiliseAtomIndexToConnectedAtomIndeces();
		this.initialiseNumberHydrogens();
		this.initiliseBondIndexToConnectedAtomIndeces();
		this.initialiseRingBondsFastBitArray();
		this.initialiseAtomAdjacencyList();
		this.initialiseAtomMasses();
	}

	/**
	 * returns bond index + 1
	 * 
	 * @param x
	 * @param y
	 * @return
	 */
	public short getBondIndexFromAtomAdjacencyList(short x, short y) {
		return this.atomAdjacencyList[this.getIndex(x, y)];
	}

	/**
	 * returns atom indeces that are connected by bond with bondIndex
	 * 
	 * @param bondIndex
	 * @return
	 */
	public short[] getConnectedAtomIndecesOfAtomIndex(short atomIndex) {
		return this.atomIndexToConnectedAtomIndeces.get(atomIndex);
	}

	public short[] getConnectedAtomIndecesOfBondIndex(short bondIndex) {
		return this.bondIndexToConnectedAtomIndeces[bondIndex];
	}

	/**
	 * convert 2D matrix coordinates to 1D adjacency list coordinate
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	private int getIndex(int a, int b) {
		int row = a;
		int col = b;
		if (a > b) {
			row = b;
			col = a;
		}
		return row * this.getNonHydrogenAtomCount() + col - ((row + 1) * (row + 2)) / 2;
	}

	public double getMassOfAtom(int index) {
		return this.massesOfAtoms[index] + this.getNumberHydrogensConnectedToAtomIndex(index)
				* Constants.MONOISOTOPIC_MASSES.get(Constants.H_INDEX).doubleValue();
	}

	public int getNonHydrogenAtomCount() {
		return this.precursorMolecule.getAtomCount();
	}

	public int getNonHydrogenBondCount() {
		return this.precursorMolecule.getBondCount();
	}

	public int getNumberHydrogensConnectedToAtomIndex(int atomIndex) {
		return this.numberHydrogensConnectedToAtom[atomIndex].intValue();
	}

	public IAtomContainer getStructureAsIAtomContainer() {
		return this.precursorMolecule;
	}

	/**
	 * initialise 1D atom adjacency list
	 */
	private void initialiseAtomAdjacencyList() {
		for (int i = 0; i < this.getNonHydrogenBondCount(); i++) {
			java.util.Iterator<IAtom> atoms = this.precursorMolecule.getBond(i).atoms().iterator();
			this.atomAdjacencyList[getIndex(this.precursorMolecule.indexOf(atoms.next()),
					this.precursorMolecule.indexOf(atoms.next()))] = (short) (i + 1);
		}
	}

	private void initialiseAtomMasses() {

		for (int i = 0; i < this.getNonHydrogenAtomCount(); i++) {
			this.massesOfAtoms[i] = Constants.MONOISOTOPIC_MASSES
					.get(Constants.ELEMENTS.indexOf(this.precursorMolecule.getAtom(i).getSymbol())).doubleValue();
		}
	}

	private void initialiseNumberHydrogens() {

		for (int i = 0; i < this.getNonHydrogenAtomCount(); i++) {
			this.numberHydrogensConnectedToAtom[i] = this.precursorMolecule.getAtom(i).getImplicitHydrogenCount();
		}
	}

	/**
	 * initialise indeces belonging to a ring in the precursor molecule
	 */
	private void initialiseRingBondsFastBitArray() throws Exception {
		final AllRingsFinder allRingsFinder = new AllRingsFinder();
		final IRingSet ringSet = allRingsFinder.findAllRings(this.precursorMolecule);

		this.initialiseRingBondToBelongingRingBondIndecesFastBitArrays(ringSet);

		if (ringSet.getAtomContainerCount() != 0) {
			final Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
			java.util.Set<IBond> aromBonds = arom.findBonds(this.precursorMolecule);
			java.util.Iterator<IBond> it = aromBonds.iterator();
			while (it.hasNext()) {
				IBond currentBond = it.next();
				this.aromaticBonds.set(this.precursorMolecule.indexOf(currentBond), true);
			}
		}
	}

	/**
	 * initialises ringBondToBelongingRingBondIndeces FastBitArray array fast and
	 * easy way to retrieve all bond indeces belonging to a ring including the bond
	 * at specified index of that array
	 * 
	 * @param ringSet
	 */
	private void initialiseRingBondToBelongingRingBondIndecesFastBitArrays(IRingSet ringSet) {

		for (int i = 0; i < this.ringBondToBelongingRingBondIndeces.length; i++)
			this.ringBondToBelongingRingBondIndeces[i] = new FastBitArray(this.precursorMolecule.getBondCount() + 1);

		for (int i = 0; i < ringSet.getAtomContainerCount(); i++) {
			int[] indexes = new int[ringSet.getAtomContainer(i).getBondCount()];
			for (int j = 0; j < ringSet.getAtomContainer(i).getBondCount(); j++) {
				indexes[j] = this.precursorMolecule.indexOf(ringSet.getAtomContainer(i).getBond(j));
			}
			for (int j = 0; j < indexes.length; j++)
				this.ringBondToBelongingRingBondIndeces[indexes[j]].setBits(indexes);
		}
	}

	/**
	 * 
	 */
	private void initiliseAtomIndexToConnectedAtomIndeces() {

		for (int i = 0; i < this.getNonHydrogenAtomCount(); i++) {
			List<IAtom> connectedAtoms = this.precursorMolecule
					.getConnectedAtomsList(this.precursorMolecule.getAtom(i));
			short[] connectedAtomIndeces = new short[connectedAtoms.size()];
			for (int k = 0; k < connectedAtoms.size(); k++)
				connectedAtomIndeces[k] = (short) this.precursorMolecule.indexOf(connectedAtoms.get(k));
			this.atomIndexToConnectedAtomIndeces.add(i, connectedAtomIndeces);
		}
	}

	/**
	 * 
	 */
	private void initiliseBondIndexToConnectedAtomIndeces() {
		for (int i = 0; i < this.getNonHydrogenBondCount(); i++) {
			this.bondIndexToConnectedAtomIndeces[i][0] = (short) this.precursorMolecule
					.indexOf(this.precursorMolecule.getBond(i).getAtom(0));
			this.bondIndexToConnectedAtomIndeces[i][1] = (short) this.precursorMolecule
					.indexOf(this.precursorMolecule.getBond(i).getAtom(1));
		}
	}
}