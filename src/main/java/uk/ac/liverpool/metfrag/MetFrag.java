/**
 * 
 */
package uk.ac.liverpool.metfrag;

import de.ipbhalle.metfraglib.interfaces.ICandidate;
import de.ipbhalle.metfraglib.list.ScoredCandidateList;
import de.ipbhalle.metfraglib.parameter.VariableNames;
import de.ipbhalle.metfraglib.process.CombinedMetFragProcess;
import de.ipbhalle.metfraglib.settings.MetFragGlobalSettings;

/**
 * 
 * @author neilswainston
 */
public class MetFrag {
	
	/**
	 * 
	 * @return String
	 * @throws Exception
	 */
	public String match() throws Exception {
		final String peakListFilePath = ClassLoader.getSystemResource("peaklist_file_example_1.txt").getFile();
		final String candidateListFilePath = ClassLoader.getSystemResource("candidate_file_example_1.txt").getFile();
		
		final MetFragGlobalSettings settings = new MetFragGlobalSettings();
		//set peaklist path and candidate list path
		settings.set(VariableNames.PEAK_LIST_PATH_NAME, peakListFilePath);
		settings.set(VariableNames.LOCAL_DATABASE_PATH_NAME, candidateListFilePath);
		//set needed parameters
		settings.set(VariableNames.RELATIVE_MASS_DEVIATION_NAME, 5.0);
		settings.set(VariableNames.ABSOLUTE_MASS_DEVIATION_NAME, 0.001);
		settings.set(VariableNames.PRECURSOR_NEUTRAL_MASS_NAME, 253.966126);
		settings.set(VariableNames.METFRAG_DATABASE_TYPE_NAME, "LocalCSV");
		
		final CombinedMetFragProcess metfragProcess = new CombinedMetFragProcess(settings);
		
		metfragProcess.retrieveCompounds();
		metfragProcess.run();
	
		final ScoredCandidateList scoredCandidateList = (ScoredCandidateList)metfragProcess.getCandidateList();
		
		ICandidate correctCandidate = null;
		
		for(int i = 0; i < scoredCandidateList.getNumberElements(); i++) {
			String inchikey1 = (String)scoredCandidateList.getElement(i).getProperty(VariableNames.INCHI_KEY_1_NAME);
			
			if(inchikey1.equals("MEFQWPUMEMWTJP")) {
				correctCandidate = scoredCandidateList.getElement(i);
			}
		}
		
		int numberPeaksUsed = scoredCandidateList.getNumberPeaksUsed();
		
		if(correctCandidate != null) {
			final int numberPeaksExplained = correctCandidate.getMatchList().getNumberElements();
			final double fragmenterScore = (Double)correctCandidate.getProperty(VariableNames.METFRAG_FRAGMENTER_SCORE_NAME);
			final double score = (Double)correctCandidate.getProperty(VariableNames.FINAL_SCORE_COLUMN_NAME);
			
			int rank = 0;
			
			for(int i = 0; i < scoredCandidateList.getNumberElements(); i++) {
				if(((Double)scoredCandidateList.getElement(i).getProperty(VariableNames.FINAL_SCORE_COLUMN_NAME)).doubleValue() >= score)
					rank++;
			}
		}
		
		
		return "MetFrag match";
	}
}
