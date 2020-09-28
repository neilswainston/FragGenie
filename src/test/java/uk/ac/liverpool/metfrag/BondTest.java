/**
 * 
 */
package uk.ac.liverpool.metfrag;

import static org.junit.Assert.assertEquals;

import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import org.junit.Test;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.interfaces.IBond.Order;

/**
 * 
 * @author neilswainston
 */
public class BondTest {

	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testEncodeDecode() throws Exception {
		final IsotopeFactory factory = Isotopes.getInstance();
		final int MIN_ATOMIC_NUM = 1;
		final int MAX_ATOMIC_NUM = 32;
		final int TESTS = 1024;
		
		for(int i = 0; i < TESTS; i++) {
			final int atomNum1 = ThreadLocalRandom.current().nextInt(MIN_ATOMIC_NUM, MAX_ATOMIC_NUM + 1);
			final int atomNum2 = ThreadLocalRandom.current().nextInt(MIN_ATOMIC_NUM, MAX_ATOMIC_NUM + 1);
			final Order order = Order.values()[ThreadLocalRandom.current().nextInt(0, Order.values().length)];
			final boolean aromatic = new Random().nextBoolean();
			
			try {
				final Bond bond = new Bond(factory.getElement(atomNum1), factory.getElement(atomNum2), order, aromatic);
				final Bond bond2 = Bond.decode(bond.encode());
				assertEquals(bond, bond2);
			}
			catch(@SuppressWarnings("unused") NullPointerException e) {
				// Take no action.
			}
			
		}
	}
}