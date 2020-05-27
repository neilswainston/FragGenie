package de.ipbhalle.metfraglib.process;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import de.ipbhalle.metfraglib.fragmenterassignerscorer.TopDownFragmenterAssignerScorer;
import de.ipbhalle.metfraglib.interfaces.ICandidate;
import de.ipbhalle.metfraglib.list.SortedTandemMassPeakList;
import de.ipbhalle.metfraglib.parameter.VariableNames;
import de.ipbhalle.metfraglib.settings.Settings;

public class CombinedSingleCandidateMetFragProcess implements Runnable {
	
	/**
	 * 
	 */
	private final ICandidate candidate;
	
	/**
	 * 
	 */
	private final SortedTandemMassPeakList peakList;
	
	/**
	 * 
	 */
	private final Settings settings;
	
	/**
	 * 
	 */
	private final Logger logger = Logger.getLogger(CombinedSingleCandidateMetFragProcess.class);

	
	/**
	 * each candidate processing thread has its own settings object and the candidate to process
	 * 
	 * @param settings
	 * @param candidate
	 */
	public CombinedSingleCandidateMetFragProcess(final Settings settings, final ICandidate candidate, final SortedTandemMassPeakList peakList, final Level logLevel) {
		this.settings = settings;
		this.settings.set(VariableNames.CANDIDATE_NAME, candidate);
		
		this.logger.setLevel(logLevel);
		this.candidate = candidate;
		this.peakList = peakList;
	}
	
	/*
	 * 
	 */
	@Override
	public void run() {
		try {
			final TopDownFragmenterAssignerScorer scorer = new TopDownFragmenterAssignerScorer(this.settings, this.candidate, this.peakList);
			scorer.calculate();
			scorer.assignInterimScoresResults();
			this.candidate.resetPrecursorMolecule();
			scorer.assignScores();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
}