package de.ipbhalle.metfraglib.interfaces;

import de.ipbhalle.metfraglib.exceptions.AtomTypeNotKnownFromInputListException;
import de.ipbhalle.metfraglib.list.FragmentList;

public interface IFragmenter {

	public FragmentList generateFragments() throws AtomTypeNotKnownFromInputListException;

	public Byte getMaximumTreeDepth();

	public void setMaximumTreeDepth(Byte maximumTreeDepth);

	public IMolecularStructure getPrecursorMolecule();

	/**
	 * delete all objects
	 */
	public void nullify();

}
