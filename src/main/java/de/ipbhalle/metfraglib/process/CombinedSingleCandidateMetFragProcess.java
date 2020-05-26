package de.ipbhalle.metfraglib.process;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import de.ipbhalle.metfraglib.fragmenterassignerscorer.AbstractFragmenterAssignerScorer;
import de.ipbhalle.metfraglib.fragmenterassignerscorer.TopDownFragmenterAssignerScorer;
import de.ipbhalle.metfraglib.interfaces.ICandidate;
import de.ipbhalle.metfraglib.list.SortedTandemMassPeakList;
import de.ipbhalle.metfraglib.parameter.VariableNames;
import de.ipbhalle.metfraglib.settings.Settings;

public class CombinedSingleCandidateMetFragProcess implements Runnable {
	
	//usually this array is of size 1
	//for special cases additional candidates may be generated out of one
	private ICandidate candidate;
	
	/**
	 * 
	 */
	private final SortedTandemMassPeakList peakList;
	
	//pre-processing candidate filter collection
	//used to check before the processing of a candidate whether it fulfills all criteria to be processed 
	// private PreProcessingCandidateFilterCollection preProcessingCandidateFilterCollection;
	//if true candidate was processed successfully if not an error might have occured
	//reference to settings object
	private Settings settings;
	
	private Logger logger = Logger.getLogger(CombinedSingleCandidateMetFragProcess.class);

	
	/**
	 * each candidate processing thread has its own settings object and the candidate to process
	 * 
	 * @param settings
	 * @param candidate
	 */
	public CombinedSingleCandidateMetFragProcess(final Settings settings, final ICandidate candidate, final SortedTandemMassPeakList peakList, final Level logLevel) {
		this.settings = settings;
		this.logger.setLevel(logLevel);
		this.candidate = candidate;
		this.peakList = peakList;
	}
	
	/**
	 * runs the single candidate metfrag process
	 * the actual processing is done with the AbstractFragmenterAssignerScorer object
	 */
	@Override
	public void run() {
		try {
			this.settings.set(VariableNames.CANDIDATE_NAME, this.candidate);
			
			//define the fragmenterAssignerScorer
			AbstractFragmenterAssignerScorer fas = new TopDownFragmenterAssignerScorer(this.settings, this.candidate, this.peakList);
			fas.setCandidate(this.candidate);
			fas.initialise();
			fas.calculate();
			fas.assignInterimScoresResults();
			this.candidate.resetPrecursorMolecule();
			fas.assignScores();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
}