package de.ipbhalle.metfraglib.fragmenterassignerscorer;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import de.ipbhalle.metfraglib.collection.ScoreCollection;
import de.ipbhalle.metfraglib.fragmenter.AbstractTopDownFragmenter;
import de.ipbhalle.metfraglib.fragmenter.TopDownNeutralLossFragmenter;
import de.ipbhalle.metfraglib.interfaces.ICandidate;
import de.ipbhalle.metfraglib.interfaces.IScore;
import de.ipbhalle.metfraglib.list.MatchList;
import de.ipbhalle.metfraglib.parameter.VariableNames;
import de.ipbhalle.metfraglib.score.NewFragmenterScore;
import de.ipbhalle.metfraglib.settings.Settings;

public abstract class AbstractFragmenterAssignerScorer {

	//matchlist contains all fragment peak matches 
	protected MatchList matchList;
	//settings object
	protected Settings settings;
	//candidate(s) to be processed
	protected ICandidate candidate;
	//fragmenter performing in silico fragmention
	protected AbstractTopDownFragmenter fragmenter;
	//score collection containg all scores
	protected ScoreCollection scoreCollection;
	//final scores is weighted sum of all scorecollection scores
	protected double finalScore;

	protected Logger logger = Logger.getLogger(AbstractFragmenterAssignerScorer.class);
	
	/**
	 * 
	 * @param settings
	 * @param candidate
	 */
	public AbstractFragmenterAssignerScorer(Settings settings, ICandidate candidate) throws Exception {
		this.settings = settings;
		this.finalScore = 0;
		this.logger.setLevel((Level)this.settings.get(VariableNames.LOG_LEVEL_NAME));
		this.candidate = candidate;

		this.candidate.initialisePrecursorCandidate();
		/*
		 * initialise fragmenter
		 */
		this.fragmenter = new TopDownNeutralLossFragmenter(this.candidate, this.settings);
		
		/*
		 * initialise score
		 */
		String[] score_types = (String[])this.settings.get(VariableNames.METFRAG_SCORE_TYPES_NAME);
		IScore[] scores = new IScore[score_types.length];
		for(int i = 0; i < score_types.length; i++) {
			logger.debug("\t\tinitialising " + score_types[i]);
			scores[i] = new NewFragmenterScore(this.settings);
		}
		this.scoreCollection = new ScoreCollection(scores);
	}
	
	public abstract void calculate();

	public void assignScores() {
		/*
		 * generate the result as scored candidate and set the scores as candidate property
		 */
		String[] score_types = (String[])this.settings.get(VariableNames.METFRAG_SCORE_TYPES_NAME);
		for(int i = 0; i < score_types.length; i++) {
			if(scoreCollection.getScore(i).getValue() != null) {
				if(!scoreCollection.getScore(i).isUserDefinedPropertyScore()) {
					this.candidate.setProperty(score_types[i], scoreCollection.getScore(i).getValue());
				}
			}
		}
	}
	
	public void assignInterimScoresResults() {
		/*
		 * generate the result as scored candidate and set the scores as candidate property
		 */
		String[] score_types = (String[])this.settings.get(VariableNames.METFRAG_SCORE_TYPES_NAME);
		for(int i = 0; i < score_types.length; i++) {
			if(scoreCollection.getScore(i).getValue() != null) {
				if(scoreCollection.getScore(i).hasInterimResults() && !scoreCollection.getScore(i).isCandidatePropertyScore()) {
					this.candidate.setProperty(score_types[i] + "_Values", scoreCollection.getScore(i).getOptimalValuesToString());
				}
			}
		}
	}
}