package de.ipbhalle.metfraglib.process;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import de.ipbhalle.metfraglib.additionals.BondEnergies;
import de.ipbhalle.metfraglib.database.LocalCSVDatabase;
import de.ipbhalle.metfraglib.interfaces.IDatabase;
import de.ipbhalle.metfraglib.interfaces.IPeakListReader;
import de.ipbhalle.metfraglib.list.CandidateList;
import de.ipbhalle.metfraglib.list.SortedTandemMassPeakList;
import de.ipbhalle.metfraglib.parameter.VariableNames;
import de.ipbhalle.metfraglib.peaklistreader.FilteredTandemMassPeakListReader;
import de.ipbhalle.metfraglib.settings.MetFragSingleProcessSettings;
import de.ipbhalle.metfraglib.settings.Settings;

public class CombinedMetFragProcess implements Runnable {

	/**
	 * Settings object containing all parameters.
	 */
	private final Settings settings;
	
	/**
	 * Candidate list -> later also containing the scored candidates.
	 */
	private CandidateList candidateList;
	
	/**
	 * 
	 */
	private final int numThreads = 1;

	/**
	 * 
	 */
	private Logger logger = Logger.getLogger(CombinedMetFragProcess.class);
	
	/**
	 * 
	 * @param settings
	 * @param logLevel
	 * @throws Exception
	 */
	public CombinedMetFragProcess(final Settings settings, final Level logLevel) throws Exception {
		this.settings = settings;
		this.settings.set(VariableNames.BOND_ENERGY_OBJECT_NAME, new BondEnergies());
		
		this.logger.setLevel(logLevel);
		
		final IDatabase database = new LocalCSVDatabase(this.settings);
		java.util.ArrayList<String> databaseCandidateIdentifiers = database.getCandidateIdentifiers();
		this.candidateList = database.getCandidateByIdentifier(databaseCandidateIdentifiers);
		
		
	}
	
	/*
	 * starts global metfrag process that starts a single thread for each candidate
	 * 
	 * @see java.lang.Runnable#run()
	 */
	@Override
	public void run() {
		try {
			final IPeakListReader peakListReader = new FilteredTandemMassPeakListReader(this.settings);
			final SortedTandemMassPeakList peakList = (SortedTandemMassPeakList)peakListReader.read();
			final CombinedSingleCandidateMetFragProcess[] processes = new CombinedSingleCandidateMetFragProcess[this.candidateList.getNumberElements()];
			
			for(int i = 0; i < this.candidateList.getNumberElements(); i++) 
			{
				this.candidateList.getElement(i).setUseSmiles(true);
				final MetFragSingleProcessSettings singleProcessSettings = new MetFragSingleProcessSettings(this.settings);
				processes[i] = new CombinedSingleCandidateMetFragProcess(singleProcessSettings, this.candidateList.getElement(i), peakList, this.logger.getLevel());
			}
			
			final ExecutorService executer = Executors.newFixedThreadPool(this.numThreads);

			for(Runnable process : processes) {
				executer.execute(process);
			}
			
			executer.shutdown();
			
		    while(!executer.isTerminated())
			{
				Thread.sleep(1000);
			}
			
			this.logger.info("Stored " + this.candidateList.getNumberElements() + " candidate(s)");
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * 
	 * @return CandidateList
	 */
	public CandidateList getCandidateList() {
		return this.candidateList;
	}
}