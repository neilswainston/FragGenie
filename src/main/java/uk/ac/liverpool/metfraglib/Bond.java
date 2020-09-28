package uk.ac.liverpool.metfraglib;

// import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IElement;

/**
 * 
 * @author neilswainston
 */
public class Bond {
	
	/**
	 * 
	 */
	private final IElement atom1;
	
	/**
	 * 
	 */
	private final IElement atom2;
	
	/**
	 * 
	 */
	private final Order order;
	
	/**
	 * 
	 */
	private final boolean aromatic;
	
	/**
	 * 
	 * @param atom1
	 * @param atom2
	 * @param order
	 * @param aromatic
	 */
	public Bond(final IElement atom1, final IElement atom2, final Order order, final boolean aromatic) {
		this.atom1 = atom1;
		this.atom2 = atom2;
		this.order = order;
		this.aromatic = aromatic;
	}

	/**
	 * 
	 * @return int
	 */
	public int encode() {
		final int atom1Encoded = this.atom1.getAtomicNumber().intValue() << 11;
		final int atom2Encoded = this.atom2.getAtomicNumber().intValue() << 4;
		final int orderEncoded = this.order.ordinal() << 1;
		final int aromaticEncoded = this.aromatic ? 1 : 0;
		return atom1Encoded + atom2Encoded + orderEncoded + aromaticEncoded;
	}
	
	/**
	 * 
	 * @param encoded
	 * @return Bond
	 */
	public static Bond decode(final int encoded) {
		// final IsotopeFactory factory = null;
		
		// final int atomicNumber1 = (encoded & 2 ^ 17) >>> 11;
		// final int atomicNumber2 = (encoded & 2 ^ 10) >>> 4;
		final int orderOrdinal = (encoded & 2 ^ 3) >>> 1;
		final boolean aromatic = (encoded & 1) == 1;
		
		return new Bond(null, null, Order.values()[orderOrdinal], aromatic);
		// return new Bond(factory.getElement(atomicNumber1), factory.getElement(atomicNumber2), Order.values()[orderOrdinal], aromatic);
	}
	
	@Override
	public String toString() {
		return Integer.toString(this.encode());
	}
}