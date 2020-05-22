package de.ipbhalle.metfraglib.process;

import java.lang.reflect.InvocationTargetException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import de.ipbhalle.metfraglib.additionals.BondEnergies;
import de.ipbhalle.metfraglib.exceptions.ScorePropertyNotDefinedException;
import de.ipbhalle.metfraglib.interfaces.ICandidate;
import de.ipbhalle.metfraglib.interfaces.IDatabase;
import de.ipbhalle.metfraglib.interfaces.IPeakListReader;
import de.ipbhalle.metfraglib.list.AbstractPeakList;
import de.ipbhalle.metfraglib.list.CandidateList;
import de.ipbhalle.metfraglib.list.DefaultPeakList;
import de.ipbhalle.metfraglib.list.ScoredCandidateList;
import de.ipbhalle.metfraglib.parameter.ClassNames;
import de.ipbhalle.metfraglib.parameter.VariableNames;
import de.ipbhalle.metfraglib.settings.MetFragGlobalSettings;
import de.ipbhalle.metfraglib.settings.MetFragSingleProcessSettings;
import de.ipbhalle.metfraglib.settings.Settings;

public class CombinedMetFragProcess implements Runnable {

	//settings object containing all parameters
	private final MetFragGlobalSettings globalSettings;
	//database object for candidate retrieval
	private IDatabase database;
	//peaklist reader to generate the peaklist -> m/z intensity 
	private IPeakListReader peakListReader; 
	//candidate filters
	// private PreProcessingCandidateFilterCollection preProcessingCandidateFilterCollection;
	// private PostProcessingCandidateFilterCollection postProcessingCandidateFilterCollection;
	//candidate list -> later also containing the scored candidates
	private CandidateList sortedScoredCandidateList;

	private boolean threadStoppedExternally = false;
	//threads to process single candidates
	private CombinedSingleCandidateMetFragProcess[] processes;
	//process status object -> stores values about metfrag's processing status
	private ProcessingStatus processingStatus;
	
	private ExecutorService executer;
	
	private Logger logger = Logger.getLogger(CombinedMetFragProcess.class);
	
	/**
	 * constructore needs settings object
	 * 
	 * @param globalSettings
	 * @throws ClassNotFoundException 
	 * @throws SecurityException 
	 * @throws NoSuchMethodException 
	 * @throws InvocationTargetException 
	 * @throws IllegalArgumentException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 */
	public CombinedMetFragProcess(final MetFragGlobalSettings globalSettings, final Level logLevel) throws InstantiationException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException, ClassNotFoundException {
		this.globalSettings = globalSettings;
		//set log level
		this.logger.setLevel(logLevel);
		
		this.processingStatus = new ProcessingStatus(this.globalSettings);
		this.globalSettings.set(VariableNames.PROCESS_STATUS_OBJECT_NAME, this.processingStatus);
		
		//initialise database
		this.database = (IDatabase) Class.forName(ClassNames.getClassNameOfDatabase((String)this.globalSettings.get(VariableNames.METFRAG_DATABASE_TYPE_NAME))).getConstructor(Settings.class).newInstance(this.globalSettings);
		
		//init peaklist reader
		this.peakListReader = (IPeakListReader) Class.forName((String)this.globalSettings.get(VariableNames.METFRAG_PEAK_LIST_READER_NAME)).getConstructor(Settings.class).newInstance(this.globalSettings);
		
		//init bond energies
		this.globalSettings.set(VariableNames.BOND_ENERGY_OBJECT_NAME, new BondEnergies());
	}
	
	/*
	 * retrieve the candidates from the database 
	 */
	public boolean retrieveCompounds() throws Exception {
		this.processes = null;
		java.util.ArrayList<String> databaseCandidateIdentifiers = this.database.getCandidateIdentifiers();
		if(this.globalSettings.containsKey(VariableNames.MAXIMUM_CANDIDATE_LIMIT_TO_STOP_NAME) && this.globalSettings.get(VariableNames.MAXIMUM_CANDIDATE_LIMIT_TO_STOP_NAME) != null) {
			int limit = (Integer)this.globalSettings.get(VariableNames.MAXIMUM_CANDIDATE_LIMIT_TO_STOP_NAME);
			if(limit < databaseCandidateIdentifiers.size()) {
				this.logger.info(databaseCandidateIdentifiers.size() + " candidate(s) exceeds the defined limit (MaxCandidateLimitToStop = " + limit + ")");
				return false;
			}
		}
		this.sortedScoredCandidateList = this.database.getCandidateByIdentifier(databaseCandidateIdentifiers);
		return true;
	}
	
	/*
	 * starts global metfrag process that starts a single thread for each candidate
	 * 
	 * @see java.lang.Runnable#run()
	 */
	@Override
	public void run() {
		this.processes = null;
		this.threadStoppedExternally = false;
			
		/*
		 * read peak list and store in settings object
		 * store database object
		 */
		try {
			this.globalSettings.set(VariableNames.PEAK_LIST_NAME, this.peakListReader.read());
		} catch (Exception e) {
			this.logger.error("Error when reading peak list.");
			e.printStackTrace();
			return;
		}
		this.globalSettings.set(VariableNames.MINIMUM_FRAGMENT_MASS_LIMIT_NAME, ((DefaultPeakList)this.globalSettings.get(VariableNames.PEAK_LIST_NAME)).getMinimumMassValue());

		this.processes = new CombinedSingleCandidateMetFragProcess[this.sortedScoredCandidateList.getNumberElements()];

		//reset processing status
		this.processingStatus.setProcessStatusString("Processing Candidates");
		this.processingStatus.setNumberCandidates(this.sortedScoredCandidateList.getNumberElements());
		this.processingStatus.setNumberFinishedCandidates(0);
		this.processingStatus.setNextPercentageValue(1);
		//initialise all necessary score parameters
		
		/*
		 * prepare single MetFrag threads
		 */
		for(int i = 0; i < this.sortedScoredCandidateList.getNumberElements(); i++) 
		{
			/*
			 * local settings for each thread stores a reference to the global settings
			 */
			MetFragSingleProcessSettings singleProcessSettings = new MetFragSingleProcessSettings(this.globalSettings);
			/*
			 * necessary to define number of hydrogens and make the implicit
			 */
			this.sortedScoredCandidateList.getElement(i).setUseSmiles((Boolean)this.globalSettings.get(VariableNames.USE_SMILES_NAME));
			CombinedSingleCandidateMetFragProcess scmfp = new CombinedSingleCandidateMetFragProcess(singleProcessSettings, this.sortedScoredCandidateList.getElement(i));
			// scmfp.setPreProcessingCandidateFilterCollection(this.preProcessingCandidateFilterCollection);
			
			this.processes[i] = scmfp;
		}
		
		/*
		 * define executer thread to run MetFrag process
		 */
		this.executer = Executors.newFixedThreadPool((Byte)this.globalSettings.get(VariableNames.NUMBER_THREADS_NAME));
		/* 
		 * ###############
		 * 	run processes
		 * ###############
		 */
		for(CombinedSingleCandidateMetFragProcess scmfp : this.processes) {
			this.executer.execute(scmfp);
		}
		this.executer.shutdown(); 
	    while(!this.executer.isTerminated())
		{
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
	    if(this.threadStoppedExternally) {
	    	return;
	    }
	    /*
	     * retrieve the result
	     */
	    ScoredCandidateList scoredCandidateList = new ScoredCandidateList();
	    if(this.processes == null) return;
	    
	    /**
	     * perform post processing of scores
	     */
	    this.globalSettings.set(VariableNames.METFRAG_PROCESSES_NAME, this.processes);
	    
		for(CombinedSingleCandidateMetFragProcess scmfp : this.processes) {
			/*
			 * check whether the single run was successful
			 */
			try {
				scmfp.assignScores();
			} catch (Exception e) {
				e.printStackTrace();
				return;
			}

			ICandidate candidate = scmfp.getScoredPrecursorCandidates();
			scoredCandidateList.addElement(candidate);
		}
		/*
		 * normalise scores of the candidate list 
		 */
		try {
			this.sortedScoredCandidateList = scoredCandidateList.normaliseScores(
				(Double[])this.globalSettings.get(VariableNames.METFRAG_SCORE_WEIGHTS_NAME), 
				(String[])this.globalSettings.get(VariableNames.METFRAG_SCORE_TYPES_NAME),
				(String[])this.globalSettings.get(VariableNames.SCORE_NAMES_NOT_TO_SCALE)
			);
		} catch (ScorePropertyNotDefinedException e) {
			this.logger.error(e.getMessage());
		}
		
		/*
		 * set number of peaks used for processing
		 */
		((ScoredCandidateList)this.sortedScoredCandidateList).setNumberPeaksUsed(((AbstractPeakList)this.globalSettings.get(VariableNames.PEAK_LIST_NAME)).getNumberPeaksUsed());
		
		this.logger.info(this.processingStatus.getNumberPreFilteredCandidates().get() + " candidate(s) were discarded before processing due to pre-filtering");
		this.logger.info(this.processingStatus.getNumberErrorCandidates().get() + " candidate(s) discarded during processing due to errors");
		// this.logger.info(this.postProcessingCandidateFilterCollection.getNumberPostFilteredCandidates() + " candidate(s) discarded after processing due to post-filtering");
		this.logger.info("Stored " + this.sortedScoredCandidateList.getNumberElements() + " candidate(s)");
		
		this.processingStatus.setProcessStatusString("Processing Candidates");
	}
	
	public CandidateList getCandidateList() {
		return this.sortedScoredCandidateList;
	}
}