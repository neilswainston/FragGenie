package uk.ac.liverpool.metfraglib.precursor;

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
import de.ipbhalle.metfraglib.interfaces.IMolecularFormula;
import de.ipbhalle.metfraglib.parameter.Constants;
import uk.ac.liverpool.metfraglib.fragment.Fragment;

public class Precursor { // implements IMolecularStructure {

	private java.util.ArrayList<short[]> atomIndexToConnectedAtomIndeces;
	private short[][] bondIndexToConnectedAtomIndeces;
	private FastBitArray[] ringBondToBelongingRingBondIndeces;
	private FastBitArray aromaticBonds;
	private short[] atomAdjacencyList;
	private byte[] numberHydrogensConnectedToAtom;
	private double[] massesOfAtoms;
	private IAtomContainer precursorMolecule;
	private double neutralMonoisotopicMass;
	private IMolecularFormula molecularFormula;

	private Precursor(IAtomContainer precursorMolecule) throws Exception {
		this.precursorMolecule = precursorMolecule;
		this.initiliseAtomIndexToConnectedAtomIndeces();
		this.initialiseNumberHydrogens();
		this.initiliseBondIndexToConnectedAtomIndeces();
		this.initialiseRingBondsFastBitArray();
		this.initialiseAtomAdjacencyList();
		this.initialiseAtomMasses();
	}

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

	public double getMeanNodeDegree() {
		double meanNodeDegree = 0.0;
		for (int i = 0; i < this.atomIndexToConnectedAtomIndeces.size(); i++) {
			meanNodeDegree += this.atomIndexToConnectedAtomIndeces.get(i).length;
		}
		meanNodeDegree /= this.atomIndexToConnectedAtomIndeces.size();
		return meanNodeDegree;
	}

	public int getNumNodeDegreeOne() {
		int numDegreeOne = 0;
		for (int i = 0; i < this.atomIndexToConnectedAtomIndeces.size(); i++) {
			numDegreeOne += this.atomIndexToConnectedAtomIndeces.get(i).length == 1 ? 1 : 0;
		}
		return numDegreeOne;
	}

	public int getNumberHydrogensConnectedToAtomIndex(int atomIndex) {
		return this.numberHydrogensConnectedToAtom[atomIndex];
	}

	/**
	 * 
	 */
	private void initiliseAtomIndexToConnectedAtomIndeces() {
		this.atomIndexToConnectedAtomIndeces = new java.util.ArrayList<>();
		for (int i = 0; i < this.getNonHydrogenAtomCount(); i++) {
			java.util.List<IAtom> connectedAtoms = this.precursorMolecule
					.getConnectedAtomsList(this.precursorMolecule.getAtom(i));
			short[] connectedAtomIndeces = new short[connectedAtoms.size()];
			for (int k = 0; k < connectedAtoms.size(); k++)
				connectedAtomIndeces[k] = (short) this.precursorMolecule.indexOf(connectedAtoms.get(k));
			this.atomIndexToConnectedAtomIndeces.add(i, connectedAtomIndeces);
		}
	}

	private void initialiseNumberHydrogens() {
		this.numberHydrogensConnectedToAtom = new byte[this.getNonHydrogenAtomCount()];
		for (int i = 0; i < this.getNonHydrogenAtomCount(); i++) {
			this.numberHydrogensConnectedToAtom[i] = (byte) (int) this.precursorMolecule.getAtom(i)
					.getImplicitHydrogenCount();
		}
	}

	private void initialiseAtomMasses() {
		this.massesOfAtoms = new double[this.getNonHydrogenAtomCount()];
		for (int i = 0; i < this.getNonHydrogenAtomCount(); i++) {
			this.massesOfAtoms[i] = Constants.MONOISOTOPIC_MASSES
					.get(Constants.ELEMENTS.indexOf(this.precursorMolecule.getAtom(i).getSymbol()));
		}
	}

	public double getMassOfAtom(int index) {
		return this.massesOfAtoms[index] + this.getNumberHydrogensConnectedToAtomIndex(index)
				* Constants.MONOISOTOPIC_MASSES.get(Constants.H_INDEX);
	}

	/**
	 * 
	 */
	private void initiliseBondIndexToConnectedAtomIndeces() {
		this.bondIndexToConnectedAtomIndeces = new short[this.getNonHydrogenBondCount()][2];

		for (int i = 0; i < this.getNonHydrogenBondCount(); i++) {
			this.bondIndexToConnectedAtomIndeces[i][0] = (short) this.precursorMolecule
					.indexOf(this.precursorMolecule.getBond(i).getAtom(0));
			this.bondIndexToConnectedAtomIndeces[i][1] = (short) this.precursorMolecule
					.indexOf(this.precursorMolecule.getBond(i).getAtom(1));
		}
	}

	/**
	 * initialise indeces belonging to a ring in the precursor molecule
	 */
	private void initialiseRingBondsFastBitArray() throws Exception {
		this.aromaticBonds = new FastBitArray(this.getNonHydrogenBondCount());
		AllRingsFinder allRingsFinder = new AllRingsFinder();
		IRingSet ringSet = allRingsFinder.findAllRings(this.precursorMolecule);
		this.initialiseRingBondToBelongingRingBondIndecesFastBitArrays(ringSet);
		if (ringSet.getAtomContainerCount() != 0) {
			Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
			java.util.Set<IBond> aromaticBonds = arom.findBonds(this.precursorMolecule);
			java.util.Iterator<IBond> it = aromaticBonds.iterator();
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
		this.ringBondToBelongingRingBondIndeces = new FastBitArray[this.precursorMolecule.getBondCount() + 1];
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
	 * initialise 1D atom adjacency list
	 */
	/*
	 * private void initialiseAtomAdjacencyList() { this.atomAdjacencyList = new
	 * short[getIndex(this.getNonHydrogenAtomCount() - 2,
	 * this.getNonHydrogenAtomCount() - 1) + 1]; for(int i = 0; i <
	 * this.getNonHydrogenAtomCount(); i++) { java.util.List<IAtom> connectedAtoms =
	 * this.precursorMolecule.getConnectedAtomsList(this.precursorMolecule.getAtom(i
	 * )); for(int k = 0; k < connectedAtoms.size(); k++) { int atomNumber =
	 * this.precursorMolecule.indexOf(connectedAtoms.get(k)); int bondNumber =
	 * this.precursorMolecule.indexOf(this.precursorMolecule.getAtom(i),
	 * connectedAtoms.get(k)); this.atomAdjacencyList[getIndex(i, atomNumber)] =
	 * (short)(bondNumber + 1); } } }
	 */

	private void initialiseAtomAdjacencyList() {
		this.atomAdjacencyList = new short[getIndex(this.getNonHydrogenAtomCount() - 2,
				this.getNonHydrogenAtomCount() - 1) + 1];
		for (int i = 0; i < this.getNonHydrogenBondCount(); i++) {
			java.util.Iterator<IAtom> atoms = this.precursorMolecule.getBond(i).atoms().iterator();
			this.atomAdjacencyList[getIndex(this.precursorMolecule.indexOf(atoms.next()),
					this.precursorMolecule.indexOf(atoms.next()))] = (short) (i + 1);
		}
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

	public Fragment toFragment() {
		return new Fragment(this);
	}

	public short getNumberHydrogens() {
		return this.molecularFormula.getNumberHydrogens();
	}

	public double getNeutralMonoisotopicMass() {
		return this.neutralMonoisotopicMass;
	}

	public int getNonHydrogenAtomCount() {
		return this.precursorMolecule.getAtomCount();
	}

	public int getNonHydrogenBondCount() {
		return this.precursorMolecule.getBondCount();
	}

	public IMolecularFormula getMolecularFormula() {
		return this.molecularFormula;
	}

	public IAtomContainer getStructureAsIAtomContainer() {
		return this.precursorMolecule;
	}
}