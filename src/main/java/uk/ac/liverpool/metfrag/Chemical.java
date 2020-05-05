/**
 * 
 */
package uk.ac.liverpool.metfrag;

/**
 * 
 * @author neilswainston
 */
public class Chemical {
	/**
	 * 
	 */
	private final String identifier;
	
	/**
	 * 
	 * @param identifier
	 */
	public Chemical(final String identifier) {
		this.identifier = identifier;
	}
	// InChI,MonoisotopicMass,CompoundName,Identifier,MolecularFormula,SMILES
	
	@Override
	public String toString() {
		return this.identifier;
	}
}