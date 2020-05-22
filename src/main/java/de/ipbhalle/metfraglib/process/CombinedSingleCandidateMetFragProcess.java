package de.ipbhalle.metfraglib.process;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import de.ipbhalle.metfraglib.fragmenterassignerscorer.AbstractFragmenterAssignerScorer;
import de.ipbhalle.metfraglib.interfaces.ICandidate;
import de.ipbhalle.metfraglib.parameter.VariableNames;
import de.ipbhalle.metfraglib.settings.Settings;

public class CombinedSingleCandidateMetFragProcess implements Runnable {
	
	//usually this array is of size 1
	//for special cases additional candidates may be generated out of one
	private ICandidate[] scoredPrecursorCandidates;
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
	public CombinedSingleCandidateMetFragProcess(Settings settings, ICandidate candidate) {
		this.settings = settings;
		this.logger.setLevel((Level)this.settings.get(VariableNames.LOG_LEVEL_NAME));
		this.scoredPrecursorCandidates = new ICandidate[] {candidate};
	}
	
	/**
	 * runs the single candidate metfrag process
	 * the actual processing is done with the AbstractFragmenterAssignerScorer object
	 */
	@Override
	public void run() {
		try {
			this.settings.set(VariableNames.CANDIDATE_NAME, this.scoredPrecursorCandidates[0]);
			
			//define the fragmenterAssignerScorer
			this.fas = (AbstractFragmenterAssignerScorer) Class.forName((String)this.settings.get(VariableNames.METFRAG_ASSIGNER_SCORER_NAME)).getConstructor(Settings.class, ICandidate.class).newInstance(this.settings, this.scoredPrecursorCandidates[0]);
			//sets the candidate to be processed
			this.fas.setCandidate(this.scoredPrecursorCandidates[0]);
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
			this.scoredPrecursorCandidates = this.fas.getCandidates();

			((ProcessingStatus)this.settings.get(VariableNames.PROCESS_STATUS_OBJECT_NAME)).checkNumberFinishedCandidates();
			this.scoredPrecursorCandidates[0].resetPrecursorMolecule();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public void singlePostCalculateScores() throws Exception {
		this.fas.singlePostCalculateScore();
	}
	
	public void assignScores() {
		this.fas.assignScores();
	}
	
	public ICandidate[] getScoredPrecursorCandidates() {
		return this.scoredPrecursorCandidates;
	}
}