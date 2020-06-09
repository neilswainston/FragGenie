package uk.ac.liverpool.metfraglib.fragment;

import de.ipbhalle.metfraglib.exceptions.AtomTypeNotKnownFromInputListException;
import de.ipbhalle.metfraglib.interfaces.IFragment;
import de.ipbhalle.metfraglib.interfaces.IMolecularStructure;

public abstract class AbstractFragment implements IFragment {

	protected int ID;
	protected byte treeDepth;
	// protected final IMolecularStructure precursorMolecule;
	protected boolean hasMatched;
	protected boolean isBestMatchedFragment;
	protected boolean isValidFragment;
	protected boolean discardedForFragmentation;

	/**
	 * constructor setting precursor molecule of fragment
	 * 
	 * @param precursor
	 */
	public AbstractFragment() {
		// this.precursorMolecule = precursorMolecule;
		this.treeDepth = 0;
		this.hasMatched = false;
		this.isValidFragment = false;
		this.discardedForFragmentation = false;
		this.isBestMatchedFragment = false;
	}

	/**
	 * 
	 * @return
	 * @throws AtomTypeNotKnownFromInputListException
	 */

	@Override
	public void setID(int ID) {
		this.ID = ID;
	}

	@Override
	public int getID() {
		return this.ID;
	}

	/**
	 * default zero
	 * 
	 * @return
	 */

	@Override
	public byte getTreeDepth() {
		return this.treeDepth;
	}

	@Override
	public boolean hasMatched() {
		return this.hasMatched;
	}

	@Override
	public void setHasMatched() {
		this.hasMatched = true;
	}

	@Override
	public void setIsBestMatchedFragment(boolean value) {
		this.isBestMatchedFragment = value;
	}

	public boolean isBestMatchedFragment() {
		return this.isBestMatchedFragment;
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

	/**
	 * 
	 * @return
	 */

	@Override
	public abstract IFragment clone(IMolecularStructure precursorMolecule);
}
