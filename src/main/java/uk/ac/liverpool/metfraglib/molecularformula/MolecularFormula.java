package uk.ac.liverpool.metfraglib.molecularformula;

import org.openscience.cdk.interfaces.IIsotope;

import de.ipbhalle.metfraglib.FastBitArray;
import de.ipbhalle.metfraglib.parameter.Constants;
import uk.ac.liverpool.metfraglib.precursor.Precursor;

/**
 * 
 * efficient representation of a molecular formula
 * 
 * @author cruttkie
 *
 */
public class MolecularFormula {

	/**
	 * 
	 */
	private final Precursor precursorMolecule;
	
	private byte[] atomsAsIndeces;
	private short[] numberOfAtoms;
	
	private short numberHydrogens;
	
	/**
	 * 
	 * @param precursorMolecule
	 * @param atomsFastBitArray
	 * @throws Exception
	 */
	public MolecularFormula(Precursor precursorMolecule, FastBitArray atomsFastBitArray) throws Exception {
		initialise(precursorMolecule);
		this.precursorMolecule = precursorMolecule;

		this.numberHydrogens = 0;
		for(int i = 0; i < this.numberOfAtoms.length; i++)
			this.numberOfAtoms[i] = 0;
		for(int i = 0; i < this.precursorMolecule.getStructureAsIAtomContainer().getAtomCount(); i++) {
			String currentAtomSymbol = this.getAtomSymbol(this.precursorMolecule.getStructureAsIAtomContainer().getAtom(i));
			byte atomNumber = (byte)Constants.ELEMENTS.indexOf(currentAtomSymbol);

			if(atomsFastBitArray.get(i)) {
				for(int ii = 0; ii < this.atomsAsIndeces.length; ii++) {
					if(this.atomsAsIndeces[ii] == atomNumber) {
						this.numberOfAtoms[ii]++;
						break;
					}
				}
				this.numberHydrogens += this.precursorMolecule.getNumberHydrogensConnectedToAtomIndex(i);
			}
		}
	}
	
	/**
	 * initialise molecular formula by molecular structure
	 * 
	 * @param molecule
	 * @throws ExplicitHydrogenRepresentationException 
	 */
	private void initialise(Precursor precursorMolecule2) 
			throws 	de.ipbhalle.metfraglib.exceptions.ExplicitHydrogenRepresentationException, 
					de.ipbhalle.metfraglib.exceptions.AtomTypeNotKnownFromInputListException {
		this.numberHydrogens = 0;
		java.util.Map<Byte, Short> elementsToCount = new java.util.HashMap<>();
		int numberElementsPresentInPrecursorMolecule = 0;
		
		for(int i = 0; i < precursorMolecule2.getStructureAsIAtomContainer().getAtomCount(); i++) {
			String currentAtomSymbol = this.getAtomSymbol(precursorMolecule2.getStructureAsIAtomContainer().getAtom(i));
			byte byteToAtomSymbol = (byte)Constants.ELEMENTS.indexOf(currentAtomSymbol);

			if(byteToAtomSymbol == -1) {
				throw new de.ipbhalle.metfraglib.exceptions.AtomTypeNotKnownFromInputListException(currentAtomSymbol + " not found"); //$NON-NLS-1$
			}
			if(currentAtomSymbol.equals("H")) //$NON-NLS-1$
				throw new de.ipbhalle.metfraglib.exceptions.ExplicitHydrogenRepresentationException();
			if(elementsToCount.containsKey(byteToAtomSymbol))
				elementsToCount.put(byteToAtomSymbol, (short)(elementsToCount.get(byteToAtomSymbol) + 1));
			else {
				elementsToCount.put(byteToAtomSymbol, (short)1);
				numberElementsPresentInPrecursorMolecule++;
			}
			this.numberHydrogens += precursorMolecule2.getStructureAsIAtomContainer().getAtom(i).getImplicitHydrogenCount();
		}
		
		java.util.Iterator<Byte> keys = elementsToCount.keySet().iterator();
		int index = 0;
		this.atomsAsIndeces = new byte[numberElementsPresentInPrecursorMolecule];
		this.numberOfAtoms = new short[numberElementsPresentInPrecursorMolecule];
		while(keys.hasNext()) {
			byte atomSymbol = keys.next();
			this.atomsAsIndeces[index] = atomSymbol;
			this.numberOfAtoms[index] = elementsToCount.get(atomSymbol);
			index++;
		}
	}

	public short getNumberHydrogens() {
		return this.numberHydrogens;
	}
	
	public void setNumberHydrogens(short numberHydrogens) {
		if(numberHydrogens >= 0) this.numberHydrogens = numberHydrogens;
	}
	
	public String getAtomSymbol(IIsotope atom) {
		String symbol = atom.getSymbol();
		if(atom.getMassNumber() != null)
			symbol = "[" + atom.getMassNumber() + symbol + "]"; //$NON-NLS-1$ //$NON-NLS-2$
		return symbol;
	}
}
