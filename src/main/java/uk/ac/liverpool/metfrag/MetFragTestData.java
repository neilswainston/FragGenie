/**
 * 
 */
package uk.ac.liverpool.metfrag;

/**
 * @author neilswainston
 *
 */
public class MetFragTestData {

	/**
	 * 
	 */
	public final static String[] SMILES = new String[] { "C(C(=O)O)OC1=NC(=C(C(=C1Cl)N)Cl)F", //$NON-NLS-1$
			"C(C(=O)OC1=NC(=C(C(=C1Cl)N)Cl)F)O", //$NON-NLS-1$
			"C1(=C(C(=NC(=C1Cl)F)C(C(=O)O)O)Cl)N", //$NON-NLS-1$
			"CC(=O)OOC1=NC(=C(C(=C1Cl)N)Cl)F", //$NON-NLS-1$
			"C(C(CSCP(=O)([O-])[O-])C(=O)[O-])C(=O)[O-]" //$NON-NLS-1$
	};

	/**
	 * 
	 */
	public final static double[] MZ = new double[] { 90.97445, 106.94476, 110.0275, 115.98965, 117.9854, 124.93547,
			124.99015, 125.99793, 133.95592, 143.98846, 144.99625, 146.0041, 151.94641, 160.96668, 163.00682, 172.99055,
			178.95724, 178.97725, 180.97293, 196.96778, 208.9678, 236.96245, 254.97312 };

	public final static int[] INTEN = new int[] { 681, 274, 110, 95, 384, 613, 146, 207, 777, 478, 352, 999, 962, 387,
			782, 17, 678, 391, 999, 720, 999, 999, 999 };
}