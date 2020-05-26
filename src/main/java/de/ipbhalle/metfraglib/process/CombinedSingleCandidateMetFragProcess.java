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
	
	//fragments the candidate, assignes fragments to m/z peaks and scores
	private AbstractFragmenterAssignerScorer fas;
	
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
			this.fas = new TopDownFragmenterAssignerScorer(this.settings, this.candidate, this.peakList);
			//sets the candidate to be processed
			this.fas.setCandidate(candidate);
			//inits the candidate, fragmenter, scores objects
			this.fas.initialise();
			
			/*
			 * do the actual work
			 * fragment candidate, assign fragments and score
			 */
			this.fas.calculate();
			//removing score assignment and shifted to CombinedMetFragProcess after postCalculating scores
			this.fas.assignInterimScoresResults();
			//set the reference to the scored candidate(s)
			this.candidate = this.fas.getCandidate();

			// ((ProcessingStatus)this.settings.get(VariableNames.PROCESS_STATUS_OBJECT_NAME)).checkNumberFinishedCandidates();
			this.candidate.resetPrecursorMolecule();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public void assignScores() {
		this.fas.assignScores();
	}
	
	public ICandidate getScoredPrecursorCandidates() {
		return candidate;
	}
}