package uk.ac.liverpool.metfrag;

import java.io.IOException;
import java.util.List;

import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.config.Isotopes;
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
	private final IElement elem1;
	
	/**
	 * 
	 */
	private final IElement elem2;
	
	/**
	 * 
	 */
	private final Order ord;
	
	/**
	 * 
	 */
	private final boolean arom;
	
	/**
	 * 
	 * @param element1
	 * @param element2
	 * @param order
	 * @param aromatic
	 */
	public Bond(final IElement element1, final IElement element2, final Order order, final boolean aromatic) {
		this.elem1 = element1.getAtomicNumber().intValue() < element2.getAtomicNumber().intValue() ? element1 : element2;
		this.elem2 = element1.getAtomicNumber().intValue() < element2.getAtomicNumber().intValue() ? element2 : element1;
		this.ord = order;
		this.arom = aromatic;
	}

	/**
	 * 
	 * @return int
	 */
	public int encode() {
		final int element1Encoded = this.elem1.getAtomicNumber().intValue() << 11;
		final int element2Encoded = this.elem2.getAtomicNumber().intValue() << 4;
		final int orderEncoded = this.ord.ordinal() << 1;
		final int aromaticEncoded = this.arom ? 1 : 0;
		return element1Encoded + element2Encoded + orderEncoded + aromaticEncoded;
	}
	
	/**
	 * 
	 * @param encoded
	 * @return Bond
	 * @throws IOException 
	 */
	public static Bond decode(final int encoded) throws IOException {
		final IsotopeFactory factory = Isotopes.getInstance();
		
		final int elementicNumber1 = (encoded & (int)Math.pow(2, 18) - 1) >>> 11;
		final int elementicNumber2 = (encoded & (int)Math.pow(2, 11) - 1) >>> 4;
		final int orderOrdinal = (encoded & (int)Math.pow(2, 4) - 1) >>> 1;
		final boolean aromatic = (encoded & 1) == 1;
		
		return new Bond(factory.getElement(elementicNumber1), factory.getElement(elementicNumber2), Order.values()[orderOrdinal], aromatic);
	}
	
	/**
	 * 
	 * @param filter
	 * @return boolean
	 */
	public boolean passesFilter(final List<Object> filter) {
		if(filter.get(0) != null && filter.get(0) instanceof List<?>) {
			final List<?> symbols = (List<?>)filter.get(0);
			
			if(symbols.get(0) != this.elem1.getSymbol()) {
				return false;
			}
			
			if(symbols.get(1) != this.elem2.getSymbol()) {
				return false;
			}
		}
		
		if(filter.get(1) != null && filter.get(1) != this.ord) {
			return false;
		}
		
		if(filter.get(2) != null && ((Boolean)filter.get(2)).booleanValue() != this.arom) {
			return false;
		}
		
		return true;
	}
	
	@Override
	public String toString() {
		return Integer.toString(this.encode());
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (this.arom ? 1231 : 1237);
		result = prime * result + ((this.elem1 == null) ? 0 : this.elem1.hashCode());
		result = prime * result + ((this.elem2 == null) ? 0 : this.elem2.hashCode());
		result = prime * result + ((this.ord == null) ? 0 : this.ord.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
			
		if (obj == null) {
			return false;
		}
			
		if (getClass() != obj.getClass()) {
			return false;
		}
			
		final Bond other = (Bond) obj;
		
		if (this.arom != other.arom) {
			return false;
		}
			
		if (this.elem1 == null) {
			if (other.elem1 != null) {
				return false;
			}	
		}
		else if (!this.elem1.equals(other.elem1)) {
			return false;
		}
			
		if (this.elem2 == null) {
			if (other.elem2 != null) {
				return false;
			}
				
		} else if (!this.elem2.equals(other.elem2)) {
			return false;
		}
			
		if (this.ord != other.ord) {
			return false;
		}
			
		return true;
	}
}