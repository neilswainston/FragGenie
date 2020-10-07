package uk.ac.liverpool.metfrag;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

public class Precursor {

	/**
	 * 
	 */
	private static final Map<String, Float> MONOISOTOPIC_MASSES = new HashMap<>();
	
	/**
	 * 
	 */
	private final IAtomContainer prec;
	
	/**
	 * 
	 */
	private final int[][] bondAtomIdxs;
	
	/**
	 * 
	 */
	private final int[][] connectedAtomIdxs;
	
	/**
	 * 
	 * @param precursor
	 */
	public Precursor(final IAtomContainer precursor) {
		this.prec = precursor;
		
		final int numAtoms = this.prec.getAtomCount();
		this.bondAtomIdxs = new int[numAtoms][numAtoms];
		this.connectedAtomIdxs = new int[numAtoms][];
		
		for (int i = 0; i < this.prec.getBondCount(); i++) {
			final Iterator<IAtom> atomsIterator = this.prec.getBond(i).atoms().iterator();
			final int idx1 = this.prec.indexOf(atomsIterator.next());
			final int idx2 = this.prec.indexOf(atomsIterator.next());
			this.bondAtomIdxs[idx1][idx2] = i;
			this.bondAtomIdxs[idx2][idx1] = i;
		}
		
		for (int i = 0; i < numAtoms; i++) {
			final List<IAtom> connected = this.prec.getConnectedAtomsList(this.prec.getAtom(i));
			final int[] connectedIdxs = new int[connected.size()];

			for (int k = 0; k < connected.size(); k++) {
				connectedIdxs[k] = this.prec.indexOf(connected.get(k));
			}
			
			this.connectedAtomIdxs[i] = connectedIdxs;
		}
	}

	/**
	 * 
	 * @return int
	 */
	int getAtomCount() {
		return this.prec.getAtomCount();
	}
	
	/**
	 * 
	 * @return int
	 */
	int getBondCount() {
		return this.prec.getBondCount();
	}
	
	/**
	 * 
	 * @param atomIdx
	 * @return String
	 */
	String getAtomSymbol(final int atomIdx) {
		final IAtom atom = this.prec.getAtom(atomIdx);
		return atom.getSymbol();
	}
	
	/**
	 * 
	 * @param atomIdx
	 * @return int
	 */
	int getAtomImplicitHydrogenCount(final int atomIdx) {
		final IAtom atom = this.prec.getAtom(atomIdx);
		return atom.getImplicitHydrogenCount().intValue();
	}
	
	/**
	 * 
	 * @param atomIdx
	 * @return int[]
	 */
	int[] getConnectedAtomIdxs(final int atomIdx) {
		return this.connectedAtomIdxs[atomIdx];
	}
	
	/**
	 * 
	 * @param atomIdx1
	 * @param atomIdx2
	 * @return int
	 */
	int getBondIdx(final int atomIdx1, final int atomIdx2) {
		return this.bondAtomIdxs[atomIdx1][atomIdx2];
	}
	
	/**
	 * 
	 * @param bondIdx
	 * @return Bond
	 */
	Bond getBond(final int bondIdx) {
		final IBond iBond = this.prec.getBond(bondIdx);
		final Iterator<IAtom> atoms = iBond.atoms().iterator();
		return new Bond(atoms.next(), atoms.next(), iBond.getOrder(), iBond.isAromatic());
	}
	
	/**
	 * 
	 * @param bondIdx
	 * @return int[]
	 */
	int[] getBondConnectedAtoms(final int bondIdx) {
		final IBond bond = this.prec.getBond(bondIdx);
		return new int[] { this.prec.indexOf(bond.getAtom(0)), this.prec.indexOf(bond.getAtom(1)) };
	}
	
	/**
	 * 
	 * @param atomIdx
	 * @return float
	 */
	float getAtomMass(final int atomIdx) {
		final IAtom atom = this.prec.getAtom(atomIdx);
		final String symbol = atom.getSymbol();
		final int hCount = atom.getImplicitHydrogenCount().intValue();
		
		return MONOISOTOPIC_MASSES.get(symbol).floatValue()
				+ hCount * MONOISOTOPIC_MASSES.get("H").floatValue(); //$NON-NLS-1$
	}
	
	static {
		MONOISOTOPIC_MASSES.put("[13C]", Float.valueOf(13.00335f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("C", Float.valueOf(12.00000f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Al", Float.valueOf(26.98154f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Am", Float.valueOf(243.06140f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ar", Float.valueOf(39.96238f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("As", Float.valueOf(74.92160f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("At", Float.valueOf(209.98710f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Au", Float.valueOf(196.96660f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("B", Float.valueOf(11.00931f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ba", Float.valueOf(137.90520f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Bi", Float.valueOf(208.98040f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Br", Float.valueOf(78.91834f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ca", Float.valueOf(39.96259f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cd", Float.valueOf(113.90340f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ce", Float.valueOf(139.90540f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cl", Float.valueOf(34.96885f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Co", Float.valueOf(58.93320f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cr", Float.valueOf(51.94051f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cu", Float.valueOf(62.92960f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("D", Float.valueOf(2.01410f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("F", Float.valueOf(18.99840f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Fe", Float.valueOf(55.93494f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ga", Float.valueOf(68.92558f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Gd", Float.valueOf(157.92410f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ge", Float.valueOf(73.92118f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("H", Float.valueOf(1.00783f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("He", Float.valueOf(4.00260f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Hg", Float.valueOf(201.97060f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("I", Float.valueOf(126.90450f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("In", Float.valueOf(114.90390f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("K", Float.valueOf(38.96371f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Li", Float.valueOf(7.01600f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Mg", Float.valueOf(23.98504f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Mn", Float.valueOf(54.93805f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Mo", Float.valueOf(97.90541f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("[15N]", Float.valueOf(15.00011f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("N", Float.valueOf(14.00307f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Na", Float.valueOf(22.98977f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ne", Float.valueOf(19.99244f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ni", Float.valueOf(57.93535f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("[18O]", Float.valueOf(17.99916f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("O", Float.valueOf(15.99491f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("P", Float.valueOf(30.97376f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Pb", Float.valueOf(207.97660f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Po", Float.valueOf(208.98240f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Pt", Float.valueOf(194.96480f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ra", Float.valueOf(226.02540f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Rb", Float.valueOf(84.91179f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Re", Float.valueOf(186.95580f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ru", Float.valueOf(101.90430f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("S", Float.valueOf(31.97207f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Sb", Float.valueOf(120.90380f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Se", Float.valueOf(79.91652f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Si", Float.valueOf(27.97693f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Sn", Float.valueOf(119.90220f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Tc", Float.valueOf(97.90722f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Te", Float.valueOf(129.90620f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ti", Float.valueOf(47.94795f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Tl", Float.valueOf(204.97440f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("U", Float.valueOf(238.05080f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("V", Float.valueOf(50.94396f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("W", Float.valueOf(183.95090f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Y", Float.valueOf(88.90585f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Zn", Float.valueOf(63.92915f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Zr", Float.valueOf(89.90470f)); //$NON-NLS-1$
	}
}