import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;


public class MainTests {
	
	public static int numBlocks =-1;
	public static int[] numAptPerBlock = null;
	public static ArrayList<ArrayList<Double>> blockLocations = new ArrayList<ArrayList<Double>>();
	public static ArrayList<int[]> capacities = new ArrayList<int[]>();
	public static int numTypes = -1;
	public static int[] numAgentPerType = null;
	public static int numAgents = 0;
	public static ArrayList<Agent> agents = null;
	public static double variance = 10;
	

	public static void generateRandomHDBData(int numBlocks, int numAptPerBlock, int numAgents, int[] numAgentsPerType, int numTypes, boolean locationBased){
		
		MainTests.numBlocks = numBlocks;
		MainTests.numAptPerBlock = new int[numBlocks];
		MainTests.numTypes = numTypes;
		
		MainTests.agents = new ArrayList<Agent>();					
		MainTests.numAgentPerType = new int[numTypes];
		MainTests.numAgents = numAgents;

		for (int l=0;l<numTypes;l++){
			MainTests.numAgentPerType[l] = numAgentsPerType[l];
		}
		MainTests.capacities = new ArrayList<int[]>();
		MainTests.blockLocations = new ArrayList<ArrayList<Double>>();
		for (int b=0;b<MainTests.numBlocks;b++){
			ArrayList<Double> loc = new ArrayList<Double>();
			loc.add(Math.random());
			loc.add(Math.random());
			MainTests.blockLocations.add(loc);
			MainTests.numAptPerBlock[b] = numAptPerBlock;
			int[] capacities_b = new int[MainTests.numTypes];
			for (int l=0;l<numTypes;l++){
				capacities_b[l]= (int)(numAptPerBlock*1.0/MainTests.numTypes);
			}
			MainTests.capacities.add(capacities_b);
		}
		if (locationBased){
			for (int l=0;l<MainTests.numTypes;l++){
				for (int i=0;i<MainTests.numAgentPerType[l];i++){
					Agent a = new Agent(l, new double[]{Math.random(),Math.random()});
					a.sortItems();
					MainTests.agents.add(a);
				}
			}
		}else{//ethnicity based
			for (int l=0;l<MainTests.numTypes;l++){
				double[] loc = new double[]{Math.random(),Math.random()};
				for (int i=0;i<MainTests.numAgentPerType[l];i++){
					Agent a = new Agent(l, loc);
					MainTests.agents.add(a);
					a.sortItems();
				}
			}
		}
	}

	public static double lotteryProcess(){
		
		double totalUtility = 0;
		
		//Generate one agent permutation:
		ArrayList<Integer> ordering = new ArrayList<Integer>();
		for (int i=0;i<numAgents;i++){
			ordering.add(i);
		}
		Collections.shuffle(ordering);
		
		//Create variables for current capacities and flat availabilities:
		ArrayList<int[]> currentCapacities = new ArrayList<int[]>();
		ArrayList<boolean[]> notAvailable = new ArrayList<boolean[]>();
		for (int b=0;b<MainTests.numBlocks;b++){
			int numApt_b = MainTests.numAptPerBlock[b];
			notAvailable.add(new boolean[numApt_b]);
			int[] currentCapacitiesb = new int[MainTests.numTypes];
			for (int l=0;l<MainTests.numTypes;l++){
				currentCapacitiesb[l] = MainTests.capacities.get(b)[l];
			}
			currentCapacities.add(currentCapacitiesb);
		}
		
		//Simulate flat selections:
		for (int i=0;i<MainTests.numAgents;i++){
			
			int agent = ordering.get(i);//this agent must choose a flat at this step
			int type = MainTests.agents.get(agent).getType();
			
			//Flats/items ordered by preference in each block:
			ArrayList<int[]> sortedItems = MainTests.agents.get(agent).getSortedItems();
			//Agent's utilities:
			ArrayList<double[]> utilities = MainTests.agents.get(agent).getUtilities();
			
			//Find best available item for the agent:
			double max = -1;
			int bestBlock = -1;
			int bestApt = -1;
			for (int b=0;b<MainTests.numBlocks;b++){
				int numApt_b = MainTests.numAptPerBlock[b];
				int pos = 0;
				if (currentCapacities.get(b)[type] > 0){//if the type capacity has not been reached in this block
					int sortedItemsb = -1;
					boolean continu = true; 
					while(pos < numApt_b && continu){//find the most preferred item among the flats that are still available in this block
						sortedItemsb = sortedItems.get(b)[pos];
						continu = notAvailable.get(b)[sortedItemsb];
						pos = pos+1;
					}
					if (!continu){//if an item can be allocated to the agent in this block
						double val = utilities.get(b)[sortedItemsb];
						if (max<val){
							max = val;
							bestBlock = b;
							bestApt = sortedItemsb;
						}
					}
				}
			}
			if (max != -1){//if the agent gets a flat
				currentCapacities.get(bestBlock)[type] = currentCapacities.get(bestBlock)[type]-1;
				notAvailable.get(bestBlock)[bestApt] = true;
				totalUtility = totalUtility + max;
			}
		}
		return totalUtility;	
	}
	
	
	public static void main(String[] args){

		String directory = "Results/";
		
		// set only one the following variables to true:
		boolean HDB_realData = false; // run experiments using the data collected from HDB website
		boolean varyHDB = false; // run experiments with real data but varying the number of agents, types and blocks
		boolean chicago = false; // run experiments using the data obtained from the Chicago public school allocation website
		
		if (HDB_realData){
			
			// set only one the following variables to true:
			boolean randomUtility = false; //random normalized utility modelies
			boolean locationBasedUtility = false; // for each agent, generate a random location LOC, and draw utilities from a normal distribution with mean = 1/dist(LOC,Block)
			boolean locationEthnicityBasedUtility = false;
			boolean incomePriceBasedUtility = false;
			
			// if not random utility, choose the variance for the utility model:
			MainTests.variance = 10;
			
			//Choose the number of runs:
			int numTests = 100;
				
			//Data obtained from the HDB website:
			double[] quotas_lambda_pq = {0.87,0.25,0.15};
			double[] actual_prop = {0.741,0.134,0.125};
			double[] income_per_agentType = {7326,4575,7664};
			
			//Load data on block locations:
			ProcessResultsFile fileLocations = new ProcessResultsFile();
			fileLocations.runParser("Data/housing/locations");
			MainTests.blockLocations = fileLocations.getVectors();
			
			//Load data on the numbers of flats:
			ProcessResultsFile fileNumApts = new ProcessResultsFile();
			fileNumApts.runParser("Data/housing/numApts");
			ArrayList<ArrayList<Double>> numApts = fileNumApts.getVectors();
			
			//Load data on the types of flats:
			ProcessResultsFile fileNumAptsPerType = new ProcessResultsFile();
			fileNumAptsPerType.runParser("Data/housing/numFlatsPerType");
			
			//Load data on the number of flats per type:
			ArrayList<ArrayList<Double>> numFlatsPerType = fileNumAptsPerType.getVectors();
			ProcessResultsFile filePricePerType = new ProcessResultsFile();
			
			//Load data on flat prices:
			filePricePerType.runParser("Data/housing/prices");
			ArrayList<ArrayList<Double>> pricePerType = filePricePerType.getVectors();
				
			//Initialization:
			MainTests.numTypes = 3; //Chinese, Malay, Indian/other
			MainTests.numBlocks = MainTests.blockLocations.size();
			MainTests.numAptPerBlock = new int[MainTests.numBlocks];
			MainTests.numAgentPerType = new int[numTypes];
			int tmp = 0;
			for (int b=0;b<MainTests.numBlocks;b++){
				MainTests.numAptPerBlock[b] = numApts.get(b).get(0).intValue();
				tmp = tmp + MainTests.numAptPerBlock[b];
				int[] capacities_b = new int[MainTests.numTypes];
				for (int l=0;l<numTypes;l++){
					capacities_b[l]= (int)(MainTests.numAptPerBlock[b]*quotas_lambda_pq[l]);
				}
				MainTests.capacities.add(capacities_b);
			}
			for (int l=0;l<MainTests.numTypes;l++){
				MainTests.numAgentPerType[l] = (int)(1350*actual_prop[l]);
				MainTests.numAgents = MainTests.numAgents + MainTests.numAgentPerType[l];
			}
			
			//Compute minAlpha for UB1:
			double minALPHA = Double.MAX_VALUE;
			for (int b=0;b<MainTests.numBlocks;b++){
				for (int l=0;l<MainTests.numTypes;l++){
					double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
					if (alpha < minALPHA){
						minALPHA = alpha;
					}
				}
			}
			
			//Compute min and max coordinates for agents' preferred locations:
			double[] minCoordinate = new double[]{Double.MAX_VALUE,Double.MAX_VALUE};
			double[] maxCoordinate = new double[]{-Double.MAX_VALUE,-Double.MAX_VALUE};
			for (int b=0;b<MainTests.numBlocks;b++){
				double locations_x =  MainTests.blockLocations.get(b).get(0).doubleValue();
				if (locations_x < minCoordinate[0]){
					minCoordinate[0] = locations_x;
				}
				if (locations_x > maxCoordinate[0]){
					maxCoordinate[0] = locations_x;
				}
				double locations_y =  MainTests.blockLocations.get(b).get(1).doubleValue();
				if (locations_y < minCoordinate[1]){
					minCoordinate[1] = locations_y;
				}
				if (locations_y > maxCoordinate[1]){
					maxCoordinate[1] = locations_y;
				}
			}
			double maxDist = Math.sqrt(Math.pow(minCoordinate[0]-maxCoordinate[0], 2) + Math.pow(minCoordinate[1]-maxCoordinate[1], 2));
			
			//Start running tests:
			for (int t=0;t<numTests;t++){

				if (randomUtility){
		
					//Generate result files:
					ProcessResultsFile fileResults = new ProcessResultsFile();
					fileResults.createResultsFile(directory+"randomUtility-var"+MainTests.variance);
					String fileError=directory+"error-randomUtility-var"+MainTests.variance+".txt";
					
					//Generate the agents (with random utilities):
					MainTests.agents = new ArrayList<Agent>();
					ArrayList<Double> results = new ArrayList<Double>();	
					for (int l=0;l<MainTests.numTypes;l++){
						for (int i=0;i<MainTests.numAgentPerType[l];i++){
							Agent a = new Agent(l); 
							a.sortItems();
							MainTests.agents.add(a);				
						}
					}	
					
					//Compute OPT:
					double optWithout = Solver.assignmentWithOrWithoutConstraintsHDB(false);
					results.add(optWithout);
					
					//Compute OPT_C:
					double optConstraints = Solver.assignmentWithOrWithoutConstraintsHDB(true);
					results.add(optWithout/optConstraints);
		
					//Compute Lot:
					double meanLot=0;
					int numLot = 50;
					for (int l=0;l<numLot;l++){
						meanLot = meanLot+ MainTests.lotteryProcess();
					}
					meanLot = meanLot*1.0/numLot;
					results.add(optWithout/meanLot);

					//Compute UB1:
					results.add(1.0/minALPHA);
							
					//Compute UB2:
					double betaStar = Solver.computeDeltaStar(optWithout);
					results.add(betaStar);
					double sum = 0;
					for (int l=0;l<MainTests.numTypes;l++){
						double minl = Double.MAX_VALUE;
						for (int b=0;b<MainTests.numBlocks;b++){
							double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
							if (alpha < minl){
								minl = alpha;
							}
						}
						sum = sum + minl * MainTests.numAgentPerType[l]*1.0/MainTests.numAgents;
					}
					results.add(1.0/(sum*betaStar));
					
					//Write results:
					fileResults.addVector(results, results.size(), 0);
					ProcessResultsFile fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"randomUtility-var"+MainTests.variance);
					fileOctave.writeOctaveMeans("totalUtilitiesBound1BetaBound2");
							
					//standard error:
					File FicError = new File(fileError);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileError);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> opts = fileResults.getVectors();
						double average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(1);
						}
						average = average*1.0/opts.size();
						double err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(1)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						String text = "constraint: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size())) +"\n" ;
								
						average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(2);
						}
						average = average*1.0/opts.size();
						err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(2)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						text = text + "lottery: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size())) +"\n";	
								
						average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(5);
						}
						average = average*1.0/opts.size();
						err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(5)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						text = text + "bound: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size()));
							
						FichierOctave.writeBytes(text);
							
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
				if (locationBasedUtility){
					
					//Generate result files:
					ProcessResultsFile fileResults = new ProcessResultsFile();
					fileResults.createResultsFile(directory+"locationBased-var"+MainTests.variance);
					String fileError=directory+"error-locationBased-var"+MainTests.variance+".txt";
				
					//Generate the agents (with location based utilities):
					MainTests.agents = new ArrayList<Agent>();					
					ArrayList<Double> results = new ArrayList<Double>();
					for (int l=0;l<MainTests.numTypes;l++){
						for (int i=0;i<MainTests.numAgentPerType[l];i++){
							double[] bestPosition = new double[]{Math.random()*(maxCoordinate[0]- minCoordinate[0]) + minCoordinate[0], Math.random()*(maxCoordinate[1]- minCoordinate[1])+minCoordinate[1]};
							Agent a = new Agent(l,bestPosition); 
							a.sortItems();
							MainTests.agents.add(a);
								
						}	
					}
					
					//Compute OPT:
					double optWithout = Solver.assignmentWithOrWithoutConstraintsHDB(false);
					results.add(optWithout);
						
					//Compute OPT_C:
					double optConstraints = Solver.assignmentWithOrWithoutConstraintsHDB(true);
					results.add(optWithout/optConstraints);
		
					//Compute Lot:
					double meanLot=0;
					int numLot = 50;
					for (int l=0;l<numLot;l++){
						meanLot = meanLot+ MainTests.lotteryProcess();
					}
					meanLot = meanLot*1.0/numLot;
					results.add(optWithout/meanLot);
						
					//Compute UB1:
					results.add(1.0/minALPHA);
						
					//Compute UB2:
					double betaStar = Solver.computeDeltaStar(optWithout);
					results.add(betaStar);
					double sum = 0;
					for (int l=0;l<MainTests.numTypes;l++){
						double minl = Double.MAX_VALUE;
						for (int b=0;b<MainTests.numBlocks;b++){
							double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
							if (alpha < minl){
								minl = alpha;
							}
						}
						sum = sum + minl * MainTests.numAgentPerType[l]*1.0/MainTests.numAgents;
					}
					results.add(1.0/(sum*betaStar));
						
					//Write results:
					fileResults.addVector(results, results.size(), 0);
					ProcessResultsFile fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"locationBased-var"+MainTests.variance);
					fileOctave.writeOctaveMeans("totalUtilitiesBound1BetaBound2");
							
					//standard error:
					File FicError = new File(fileError);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileError);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> opts = fileResults.getVectors();
						double average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(1);
						}
						average = average*1.0/opts.size();
						double err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(1)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						String text = "constraint: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size())) +"\n" ;
								
						average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(2);
						}
						average = average*1.0/opts.size();
						err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(2)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						text = text + "lottery: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size())) +"\n";
								
								
						average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(5);
						}
						average = average*1.0/opts.size();
						err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(5)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						text = text + "bound: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size()));
							
						FichierOctave.writeBytes(text);
							
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
				if (locationEthnicityBasedUtility){
					
					//Generate result files:
					ProcessResultsFile fileResults = new ProcessResultsFile();
					fileResults.createResultsFile(directory+"locationEthnicityBased-var"+MainTests.variance);
					String fileError=directory+"error-locationEthnicityBased-var"+MainTests.variance+".txt";
					
					//Generate the agents (with ethnicity based utilities):
					MainTests.agents = new ArrayList<Agent>();
					ArrayList<Double> results = new ArrayList<Double>();
					for (int l=0;l<MainTests.numTypes;l++){
						double[] bestPosition = new double[]{Math.random()*(maxCoordinate[0]- minCoordinate[0]) + minCoordinate[0], Math.random()*(maxCoordinate[1]- minCoordinate[1])+minCoordinate[1]};
						for (int i=0;i<MainTests.numAgentPerType[l];i++){
							Agent a = new Agent(l,bestPosition); 
							a.sortItems();
							MainTests.agents.add(a);
						}
					}
					
					//Compute OPT:
					double optWithout = Solver.assignmentWithOrWithoutConstraintsHDB(false);
					results.add(optWithout);
						
					//Compute OPT_C:
					double optConstraints = Solver.assignmentWithOrWithoutConstraintsHDB(true);
					results.add(optWithout/optConstraints);

					//Compute Lot:
					double meanLot=0;
					int numLot = 50;
					for (int l=0;l<numLot;l++){
						meanLot = meanLot+ MainTests.lotteryProcess();
					}
					meanLot = meanLot*1.0/numLot;
					results.add(optWithout/meanLot);
							
					//Compute UB1:
					results.add(1.0/minALPHA);
						
					//Compute UB2:
					double betaStar = Solver.computeDeltaStar(optWithout);
					results.add(betaStar);
					double sum = 0;
					for (int l=0;l<MainTests.numTypes;l++){
						double minl = Double.MAX_VALUE;
						for (int b=0;b<MainTests.numBlocks;b++){
							double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
							if (alpha < minl){
								minl = alpha;
							}
						}
						sum = sum + minl * MainTests.numAgentPerType[l]*1.0/MainTests.numAgents;
					}
					results.add(1.0/(sum*betaStar));
							
					//Write results:
					fileResults.addVector(results, results.size(), 0);
					ProcessResultsFile fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"locationEthnicityBased-var"+MainTests.variance);
					fileOctave.writeOctaveMeans("totalUtilitiesBound1BetaBound2");
							
					//standard error:
					File FicError = new File(fileError);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileError);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> opts = fileResults.getVectors();
						double average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(1);
						}
						average = average*1.0/opts.size();
						double err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(1)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						String text = "constraint: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size())) +"\n" ;
									
						average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(2);
						}
						average = average*1.0/opts.size();
						err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(2)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						text = text + "lottery: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size())) +"\n";
												
						average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(5);
						}
						average = average*1.0/opts.size();
						err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(5)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						text = text + "bound: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size()));
								
						FichierOctave.writeBytes(text);
								
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
				if (incomePriceBasedUtility){
					
					//Generate the result files:
					ProcessResultsFile fileResults = new ProcessResultsFile();
					fileResults.createResultsFile(directory+"incomePriceBasedUtility-var"+MainTests.variance);
					String fileError=directory+"error-incomePriceBasedUtility-var"+MainTests.variance+".txt";
					
					//Generate the agents (with income based utilities):
					MainTests.agents = new ArrayList<Agent>();
					ArrayList<Double> results = new ArrayList<Double>();
					ArrayList<ArrayList<Double>> prices = new ArrayList<ArrayList<Double>>();
					for (int b=0;b<MainTests.numBlocks;b++){
						ArrayList<Double> price_b = new ArrayList<Double>();
						int numTypes = numFlatsPerType.get(b).size();
						ArrayList<Double> true_price_b = pricePerType.get(b);
						for (int type=0;type<numTypes;type++){
							double numFlats_t = numFlatsPerType.get(b).get(type);
							for (int flat=0;flat<numFlats_t;flat++){
								price_b.add((Math.random()*(true_price_b.get(2*type+1)-true_price_b.get(2*type) + true_price_b.get(2*type)))/120);//Price per month during 10 years
							}
						}
						prices.add(price_b);
					}
					Random r = new Random();
					for (int l=0;l<MainTests.numTypes;l++){
						for (int i=0;i<MainTests.numAgentPerType[l];i++){
							double income_i = r.nextGaussian()*Math.sqrt(MainTests.variance)+income_per_agentType[l];
							Agent a = new Agent(l,income_i,prices); 
							a.sortItems();
							MainTests.agents.add(a);
						}
					}
					
					//Compute OPT:
					double optWithout = Solver.assignmentWithOrWithoutConstraintsHDB(false);
					results.add(optWithout);
							
					//Compute OPT_C:
					double optConstraints = Solver.assignmentWithOrWithoutConstraintsHDB(true);
					results.add(optWithout/optConstraints);
		
					//Compute Lot:
					double meanLot=0;
					int numLot = 50;
					for (int l=0;l<numLot;l++){
						meanLot = meanLot+ MainTests.lotteryProcess();
					}
					meanLot = meanLot*1.0/numLot;
					results.add(optWithout/meanLot);
						
					//Compute UB1:
					results.add(1.0/minALPHA);
						
					//Compute UB2:
					double betaStar = Solver.computeDeltaStar(optWithout);
					results.add(betaStar);
					double sum = 0;
					for (int l=0;l<MainTests.numTypes;l++){
						double minl = Double.MAX_VALUE;
						for (int b=0;b<MainTests.numBlocks;b++){
							double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
							if (alpha < minl){
								minl = alpha;
							}
						}
						sum = sum + minl * MainTests.numAgentPerType[l]*1.0/MainTests.numAgents;
					}
					results.add(1.0/(sum*betaStar));
							
					//Write results:
					fileResults.addVector(results, results.size(), 0);
					ProcessResultsFile fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"incomePriceBasedUtility-var"+MainTests.variance);
					fileOctave.writeOctaveMeans("totalUtilitiesBound1BetaBound2");
							
					//standard error:
					File FicError = new File(fileError);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileError);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> opts = fileResults.getVectors();
						double average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(1);
						}
						average = average*1.0/opts.size();
						double err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(1)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						String text = "constraint: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size())) +"\n" ;
								
						average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(2);
						}
						average = average*1.0/opts.size();
						err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(2)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						text = text + "lottery: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size())) +"\n";
									
						average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(5);
						}
						average = average*1.0/opts.size();
						err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(5)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						text = text + "bound: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size()));
							
						FichierOctave.writeBytes(text);
							
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}	
		
		if (varyHDB){
			
			//set only one of the following variables to true:
			boolean varyNumberAgents = false;// run experiments varying the number of agents
			boolean varyNumberBlocks = false;// run experiments varying the number of blocks
			boolean varyNumberTypes = false; // run experiments varying the number of agent types
			
			//Define the utility model:
			boolean locationBased = true;//if false, then ethnicityLocationBased
			MainTests.variance = 10;
			
			//Choose the number of runs:
			int numTests = 100;
		
			if (varyNumberAgents){
				
				//Define how to vary the number of agents:
				int[] currentNumAgents = {500,1000,2000,3000};
				//Choose the number of blocks, the number of flat per block, and the number of types:
				int numBlocks = 10;
				int numAptPerBlock = 100;
				int numTypes = 5;
	
				//Generate the result files:
				ProcessResultsFile fileResultsOpt = new ProcessResultsFile();
				String fileErrorOpt="";
				ProcessResultsFile fileResultsOptConst = new ProcessResultsFile();
				String fileErrorOptConst="";
				ProcessResultsFile fileResultsLottery = new ProcessResultsFile();
				String fileErrorLottery="";
				ProcessResultsFile fileResultsUpBound1 = new ProcessResultsFile();
				String fileErrorUpBound1="";
				ProcessResultsFile fileResultsBetaStar = new ProcessResultsFile();
				ProcessResultsFile fileResultsUpBound2 = new ProcessResultsFile();
				String fileErrorUpBound2="";
				
				fileResultsOpt.createResultsFile(directory+"Opt-randomData-varyNumberAgents-var"+MainTests.variance);
				fileErrorOpt=directory+"error-Opt-randomData-varyNumberAgents-var"+MainTests.variance+".txt";
				
				fileResultsOptConst.createResultsFile(directory+"OptConst-randomData-varyNumberAgents-var"+MainTests.variance);
				fileErrorOptConst=directory+"error-OptConst-randomData-varyNumberAgents-var"+MainTests.variance+".txt";
				
				fileResultsLottery.createResultsFile(directory+"Lottery-randomData-varyNumberAgents-var"+MainTests.variance);
				fileErrorLottery=directory+"error-Lottery-randomData-varyNumberAgents-var-"+MainTests.variance+".txt";
					
				fileResultsBetaStar.createResultsFile(directory+"BetaStar-randomData-varyNumberAgents-var"+MainTests.variance);
					
				fileResultsUpBound1.createResultsFile(directory+"UpBound1-randomData-varyNumberAgents-var"+MainTests.variance);
				fileErrorUpBound1=directory+"error-UpBound1-randomData-varyNumberAgents-var-"+MainTests.variance+".txt";
					
				fileResultsUpBound2.createResultsFile(directory+"UpBound2-randomData-varyNumberAgents-var"+MainTests.variance);
				fileErrorUpBound2=directory+"error-UpBound2-randomData-varyNumberAgents-var"+MainTests.variance+".txt";
	
				//Start running tests:
				for (int t=0;t<numTests;t++){
					ArrayList<Double> Opt = new ArrayList<Double> ();
					ArrayList<Double> OptConst = new ArrayList<Double> ();
					ArrayList<Double> Lottery = new ArrayList<Double>();
					ArrayList<Double> upBound1 = new ArrayList<Double> ();
					ArrayList<Double> BetaStar = new ArrayList<Double> ();
					ArrayList<Double> UpBound2 = new ArrayList<Double> ();
	
					//Vary the number of agents:
					for (int i=0;i<currentNumAgents.length;i++){
					
						//Generate random data: 
						int numAgents = currentNumAgents[i];
						int[] numAgentsPerType = new int[numTypes];
						int sum = 0;
						for (int j=0;j<numTypes-1;j++){
							numAgentsPerType[j] = numAgents/numTypes;
							sum = sum + numAgentsPerType[j];
						}
						numAgentsPerType[numTypes-1] = numAgentsPerType[numTypes-1]+numAgents-sum;	
						MainTests.generateRandomHDBData(numBlocks,numAptPerBlock,numAgents,numAgentsPerType,numTypes, locationBased);
						
						//Compute OPT
						double optWithout = Solver.assignmentWithOrWithoutConstraintsHDB(false);
						Opt.add(optWithout);
						
						//Compute OPT_C
						double optConstraints = Solver.assignmentWithOrWithoutConstraintsHDB(true);
						OptConst.add(optConstraints);
	
						//Compute Lot
						double meanLot=0;
						int numLot = 100;
						for (int l=0;l<numLot;l++){
							meanLot = meanLot+ MainTests.lotteryProcess();
						}
						meanLot = meanLot*1.0/numLot;
						Lottery.add(meanLot);
						
						//Compute UB1:
						double minALPHA = Double.MAX_VALUE;
						for (int b=0;b<MainTests.numBlocks;b++){
							for (int l=0;l<MainTests.numTypes;l++){
								double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
								if (alpha < minALPHA){
									minALPHA = alpha;
								}
							}
						}	
						upBound1.add(1.0/minALPHA);
						
						//Compute UB2:
						double betaStar = Solver.computeDeltaStar(optWithout);
						BetaStar.add(betaStar);
						double summ = 0;
						for (int l=0;l<MainTests.numTypes;l++){
							double minl = Double.MAX_VALUE;
							for (int b=0;b<MainTests.numBlocks;b++){
								double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
								if (alpha < minl){
									minl = alpha;
								}
							}
							summ = summ + minl * MainTests.numAgentPerType[l]*1.0/MainTests.numAgents;
						}
						UpBound2.add(1.0/(summ*betaStar));
					}
	
					//Write the results:
					fileResultsOpt.addVector(Opt, Opt.size(), 0);
					ProcessResultsFile fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"Opt-randomData-varyNumberAgents-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("OptWithoutConstraint");
						
					fileResultsOptConst.addVector(OptConst, OptConst.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"OptConst-randomData-varyNumberAgents-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("OptConstraint");
						
					fileResultsLottery.addVector(Lottery, Lottery.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"Lottery-randomData-varyNumberAgents-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("Lottery");
						
					fileResultsUpBound1.addVector(upBound1, upBound1.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"UpBound1-randomData-varyNumberAgents-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("UpBound1");
						
					fileResultsBetaStar.addVector(BetaStar, BetaStar.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"BetaStar-randomData-varyNumberAgents-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("BetaStar");
						
					fileResultsUpBound2.addVector(UpBound2, UpBound2.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"UpBound2-randomData-varyNumberAgents-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("UpBound2");
					
					//standard error:
					File FicError = new File(fileErrorOpt);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> opts = fileResultsOpt.getVectors();
						String text = "";
						for (int k=0;k<currentNumAgents.length;k++){
							double average = 0;
							for (int l=0;l<opts.size();l++){
								average = average +opts.get(l).get(k);
							}
							average = average*1.0/opts.size();
							double err=0;
							for (int l=0;l<opts.size();l++){
								err = err +Math.pow(opts.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/opts.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					FicError = new File(fileErrorOptConst);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> optsConst = fileResultsOptConst.getVectors();
						String text = "";
						for (int k=0;k<currentNumAgents.length;k++){
							double average = 0;
							for (int l=0;l<optsConst.size();l++){
								average = average +optsConst.get(l).get(k);
							}
							average = average*1.0/optsConst.size();
							double err=0;
							for (int l=0;l<optsConst.size();l++){
								err = err +Math.pow(optsConst.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/optsConst.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					FicError = new File(fileErrorLottery);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> lottery = fileResultsLottery.getVectors();
						String text = "";
						for (int k=0;k<currentNumAgents.length;k++){
							double average = 0;
							for (int l=0;l<lottery.size();l++){
								average = average +lottery.get(l).get(k);
							}
							average = average*1.0/lottery.size();
							double err=0;
							for (int l=0;l<lottery.size();l++){
								err = err +Math.pow(lottery.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/lottery.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}	
						
					FicError = new File(fileErrorUpBound1);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> UpBound1 = fileResultsUpBound1.getVectors();
						String text = "";
						for (int k=0;k<currentNumAgents.length;k++){
							double average = 0;
							for (int l=0;l<UpBound1.size();l++){
								average = average +UpBound1.get(l).get(k);
							}
							average = average*1.0/UpBound1.size();
							double err=0;
							for (int l=0;l<UpBound1.size();l++){
								err = err +Math.pow(UpBound1.get(l).get(k)-average,2);
							}
							err= Math.sqrt(err/UpBound1.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}	
					
					FicError = new File(fileErrorUpBound2);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> upBound2 = fileResultsUpBound2.getVectors();
						String text = "";
						for (int k=0;k<currentNumAgents.length;k++){
							double average = 0;
							for (int l=0;l<upBound2.size();l++){
								average = average +upBound2.get(l).get(k);
							}
							average = average*1.0/upBound2.size();
							double err=0;
							for (int l=0;l<upBound2.size();l++){
								err = err +Math.pow(upBound2.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/upBound2.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}	
				}
			}	
			
			if (varyNumberBlocks){
				
				//Define how to vary the number of blocks:
				int[] currentNumBlocks = {5,10,15,20};
				//Choose the number of flats per block, the number of agents and the number of type
				int numAptPerBlock = 100;
				int numAgents = 2000;
				int numTypes = 5;
					
				//Generate result files:
				ProcessResultsFile fileResultsOpt = new ProcessResultsFile();
				String fileErrorOpt="";
				ProcessResultsFile fileResultsOptConst = new ProcessResultsFile();
				String fileErrorOptConst="";
				ProcessResultsFile fileResultsLottery = new ProcessResultsFile();
				String fileErrorLottery="";
				ProcessResultsFile fileResultsUpBound1 = new ProcessResultsFile();
				String fileErrorUpBound1="";
				ProcessResultsFile fileResultsBetaStar = new ProcessResultsFile();
				ProcessResultsFile fileResultsUpBound2 = new ProcessResultsFile();
				String fileErrorUpBound2="";
				
				fileResultsOpt.createResultsFile(directory+"Opt-randomData-varyNumberBlocks-var"+MainTests.variance);
				fileErrorOpt=directory+"error-Opt-randomData-varyNumberBlocks-var"+MainTests.variance+".txt";
					
				fileResultsOptConst.createResultsFile(directory+"OptConst-randomData-varyNumberBlocks-var"+MainTests.variance);
				fileErrorOptConst=directory+"error-OptConst-randomData-varyNumberBlocks-var"+MainTests.variance+".txt";
					
				fileResultsLottery.createResultsFile(directory+"Lottery-randomData-varyNumberBlocks-var"+MainTests.variance);
				fileErrorLottery=directory+"error-Lottery-randomData-varyNumberBlocks-var-"+MainTests.variance+".txt";
					
				fileResultsBetaStar.createResultsFile(directory+"BetaStar-randomData-varyNumberBlocks-var"+MainTests.variance);
					
				fileResultsUpBound1.createResultsFile(directory+"UpBound1-randomData-varyNumberBlocks-var"+MainTests.variance);
				fileErrorUpBound1=directory+"error-UpBound1-randomData-varyNumberBlocks-var-"+MainTests.variance+".txt";
					
				fileResultsUpBound2.createResultsFile(directory+"UpBound2-randomData-varyNumberBlocks-var"+MainTests.variance);
				fileErrorUpBound2=directory+"error-UpBound2-randomData-varyNumberBlocks-var"+MainTests.variance+".txt";
	
				//Start running tests:
				for (int t=0;t<numTests;t++){
					ArrayList<Double> Opt = new ArrayList<Double> ();
					ArrayList<Double> OptConst = new ArrayList<Double> ();
					ArrayList<Double> Lottery = new ArrayList<Double>();
					ArrayList<Double> upBound1 = new ArrayList<Double> ();
					ArrayList<Double> BetaStar = new ArrayList<Double> ();
					ArrayList<Double> UpBound2 = new ArrayList<Double> ();
	
					//Vary the number of blocks:
					for (int i=0;i<currentNumBlocks.length;i++){
						
						//Generate random data:
						int numBlocks = currentNumBlocks[i];
						int[] numAgentsPerType = new int[numTypes];
						int sum = 0;
						for (int j=0;j<numTypes-1;j++){
							numAgentsPerType[j] = numAgents/numTypes;
							sum = sum + numAgentsPerType[j];
						}
						numAgentsPerType[numTypes-1] = numAgentsPerType[numTypes-1]+numAgents-sum;
						MainTests.generateRandomHDBData(numBlocks,numAptPerBlock,numAgents,numAgentsPerType,numTypes, locationBased);
						
						//Compute OPT:
						double optWithout = Solver.assignmentWithOrWithoutConstraintsHDB(false);
						Opt.add(optWithout);
						
						//Compute OPT_C:
						double optConstraints = Solver.assignmentWithOrWithoutConstraintsHDB(true);
						OptConst.add(optConstraints);
	
						//Compute Lot:
						double meanLot=0;
						int numLot = 50;
						for (int l=0;l<numLot;l++){
							meanLot = meanLot+ MainTests.lotteryProcess();
						}
						meanLot = meanLot*1.0/numLot;
						Lottery.add(meanLot);
						
						//Compute UB1:
						double minALPHA = Double.MAX_VALUE;
						for (int b=0;b<MainTests.numBlocks;b++){
							for (int l=0;l<MainTests.numTypes;l++){
								double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
								if (alpha < minALPHA){
									minALPHA = alpha;
								}
							}
						}	
						upBound1.add(1.0/minALPHA);
						
						//Compute UB2:
						double betaStar = Solver.computeDeltaStar(optWithout);
						BetaStar.add(betaStar);
						double summ = 0;
						for (int l=0;l<MainTests.numTypes;l++){
							double minl = Double.MAX_VALUE;
							for (int b=0;b<MainTests.numBlocks;b++){
								double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
								if (alpha < minl){
									minl = alpha;
								}
							}
							summ = summ + minl * MainTests.numAgentPerType[l]*1.0/MainTests.numAgents;
						}
						UpBound2.add(1.0/(summ*betaStar));
					}
					
					//Write the results:
					fileResultsOpt.addVector(Opt, Opt.size(), 0);
					ProcessResultsFile fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"Opt-randomData-varyNumberBlocks-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("OptWithoutConstraint");
						
					fileResultsOptConst.addVector(OptConst, OptConst.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"OptConst-randomData-varyNumberBlocks-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("OptConstraint");
						
					fileResultsLottery.addVector(Lottery, Lottery.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"Lottery-randomData-varyNumberBlocks-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("Lottery");
						
					fileResultsUpBound1.addVector(upBound1, upBound1.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"UpBound1-randomData-varyNumberBlocks-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("UpBound1");
						
					fileResultsBetaStar.addVector(BetaStar, BetaStar.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"BetaStar-randomData-varyNumberBlocks-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("BetaStar");
						
					fileResultsUpBound2.addVector(UpBound2, UpBound2.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"UpBound2-randomData-varyNumberBlocks-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("UpBound2");
					
					//standard error
					File FicError = new File(fileErrorOpt);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> opts = fileResultsOpt.getVectors();
						String text = "";
						for (int k=0;k<currentNumBlocks.length;k++){
							double average = 0;
							for (int l=0;l<opts.size();l++){
								average = average +opts.get(l).get(k);
							}
							average = average*1.0/opts.size();
							double err=0;
							for (int l=0;l<opts.size();l++){
								err = err +Math.pow(opts.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/opts.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					FicError = new File(fileErrorOptConst);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> optsConst = fileResultsOptConst.getVectors();
						String text = "";
						for (int k=0;k<currentNumBlocks.length;k++){
							double average = 0;
							for (int l=0;l<optsConst.size();l++){
								average = average +optsConst.get(l).get(k);
							}
							average = average*1.0/optsConst.size();
							double err=0;
							for (int l=0;l<optsConst.size();l++){
								err = err +Math.pow(optsConst.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/optsConst.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					FicError = new File(fileErrorLottery);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> lottery = fileResultsLottery.getVectors();
						String text = "";
						for (int k=0;k<currentNumBlocks.length;k++){
							double average = 0;
							for (int l=0;l<lottery.size();l++){
								average = average +lottery.get(l).get(k);
							}
							average = average*1.0/lottery.size();
							double err=0;
							for (int l=0;l<lottery.size();l++){
								err = err +Math.pow(lottery.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/lottery.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
							
					FicError = new File(fileErrorUpBound1);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> UpBound1 = fileResultsUpBound1.getVectors();
						String text = "";
						for (int k=0;k<currentNumBlocks.length;k++){
							double average = 0;
							for (int l=0;l<UpBound1.size();l++){
								average = average +UpBound1.get(l).get(k);
							}
							average = average*1.0/UpBound1.size();
							double err=0;
							for (int l=0;l<UpBound1.size();l++){
								err = err +Math.pow(UpBound1.get(l).get(k)-average,2);
							}
							err= Math.sqrt(err/UpBound1.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}	
					
					FicError = new File(fileErrorUpBound2);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> upBound2 = fileResultsUpBound2.getVectors();
						String text = "";
						for (int k=0;k<currentNumBlocks.length;k++){
							double average = 0;
							for (int l=0;l<upBound2.size();l++){
								average = average +upBound2.get(l).get(k);
							}
							average = average*1.0/upBound2.size();
							double err=0;
							for (int l=0;l<upBound2.size();l++){
								err = err +Math.pow(upBound2.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/upBound2.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}	
				}
			}
			
			if (varyNumberTypes){
				
				//Define how to vary the number of types:
				int[] currentNumTypes = {5,10,15,20};
				//Choose the number of blocks, the number of flats per block and the number of agents:
				int numBlocks = 10;
				int numAptPerBlock = 100;
				int numAgents = 2000;
	
				//Generate the result files:
				ProcessResultsFile fileResultsOpt = new ProcessResultsFile();
				String fileErrorOpt="";
				ProcessResultsFile fileResultsOptConst = new ProcessResultsFile();
				String fileErrorOptConst="";
				ProcessResultsFile fileResultsLottery = new ProcessResultsFile();
				String fileErrorLottery="";
				ProcessResultsFile fileResultsUpBound1 = new ProcessResultsFile();
				String fileErrorUpBound1="";
				ProcessResultsFile fileResultsBetaStar = new ProcessResultsFile();
				ProcessResultsFile fileResultsUpBound2 = new ProcessResultsFile();
				String fileErrorUpBound2="";
				
				fileResultsOpt.createResultsFile(directory+"Opt-randomData-varyNumberTypes-var"+MainTests.variance);
				fileErrorOpt=directory+"error-Opt-randomData-varyNumberTypes-var"+MainTests.variance+".txt";
					
				fileResultsOptConst.createResultsFile(directory+"OptConst-randomData-varyNumberTypes-var"+MainTests.variance);
				fileErrorOptConst=directory+"error-OptConst-randomData-varyNumberTypes-var"+MainTests.variance+".txt";
					
				fileResultsLottery.createResultsFile(directory+"Lottery-randomData-varyNumberTypes-var"+MainTests.variance);
				fileErrorLottery=directory+"error-Lottery-randomData-varyNumberTypes-var-"+MainTests.variance+".txt";
					
				fileResultsBetaStar.createResultsFile(directory+"BetaStar-randomData-varyNumberTypes-var"+MainTests.variance);
					
				fileResultsUpBound1.createResultsFile(directory+"UpBound1-randomData-varyNumberTypes-var"+MainTests.variance);
				fileErrorUpBound1=directory+"error-UpBound1-randomData-varyNumberTypes-var-"+MainTests.variance+".txt";
					
				fileResultsUpBound2.createResultsFile(directory+"UpBound2-randomData-varyNumberTypes-var"+MainTests.variance);
				fileErrorUpBound2=directory+"error-UpBound2-randomData-varyNumberTypes-var"+MainTests.variance+".txt";
				
				//Start running tests:
				for (int t=0;t<numTests;t++){
					ArrayList<Double> Opt = new ArrayList<Double> ();
					ArrayList<Double> OptConst = new ArrayList<Double> ();
					ArrayList<Double> Lottery = new ArrayList<Double>();
					ArrayList<Double> upBound1 = new ArrayList<Double> ();
					ArrayList<Double> BetaStar = new ArrayList<Double> ();
					ArrayList<Double> UpBound2 = new ArrayList<Double> ();
	
					//Vary the number of agent types:
					for (int i=0;i<currentNumTypes.length;i++){
						
						//Generate the random data:
						int numTypes = currentNumTypes[i];
						int[] numAgentsPerType = new int[numTypes];
						int sum = 0;
						for (int j=0;j<numTypes-1;j++){
							numAgentsPerType[j] = numAgents/numTypes;
							sum = sum + numAgentsPerType[j];
						}
						numAgentsPerType[numTypes-1] = numAgentsPerType[numTypes-1]+numAgents-sum;	
						MainTests.generateRandomHDBData(numBlocks,numAptPerBlock,numAgents,numAgentsPerType,numTypes, locationBased);	
						
						//Compute OPT:
						double optWithout = Solver.assignmentWithOrWithoutConstraintsHDB(false);
						Opt.add(optWithout);
						
						//Compute OPT_C:
						double optConstraints = Solver.assignmentWithOrWithoutConstraintsHDB(true);
						OptConst.add(optConstraints);
	
						//Compute Lot:
						double meanLot=0;
						int numLot = 50;
						for (int l=0;l<numLot;l++){
							meanLot = meanLot+ MainTests.lotteryProcess();
						}
						meanLot = meanLot*1.0/numLot;
						Lottery.add(meanLot);
						
						//Compute UB1:
						double minALPHA = Double.MAX_VALUE;
						for (int b=0;b<MainTests.numBlocks;b++){
							for (int l=0;l<MainTests.numTypes;l++){
								double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
								if (alpha < minALPHA){
									minALPHA = alpha;
								}
							}
						}	
						upBound1.add(1.0/minALPHA);
						
						//Compute UB2:
						double betaStar = Solver.computeDeltaStar(optWithout);
						BetaStar.add(betaStar);
						
						double summ = 0;
						for (int l=0;l<MainTests.numTypes;l++){
							double minl = Double.MAX_VALUE;
							for (int b=0;b<MainTests.numBlocks;b++){
								double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
								if (alpha < minl){
									minl = alpha;
								}
							}
							summ = summ + minl * MainTests.numAgentPerType[l]*1.0/MainTests.numAgents;
						}
						UpBound2.add(1.0/(summ*betaStar));
					}
					
					//Write the results:
					fileResultsOpt.addVector(Opt, Opt.size(), 0);
					ProcessResultsFile fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"Opt-randomData-varyNumberTypes-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("OptWithoutConstraint");
						
					fileResultsOptConst.addVector(OptConst, OptConst.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"OptConst-randomData-varyNumberTypes-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("OptConstraint");
						
					fileResultsLottery.addVector(Lottery, Lottery.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"Lottery-randomData-varyNumberTypes-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("Lottery");
						
					fileResultsUpBound1.addVector(upBound1, upBound1.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"UpBound1-randomData-varyNumberTypes-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("UpBound1");
						
					fileResultsBetaStar.addVector(BetaStar, BetaStar.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"BetaStar-randomData-varyNumberTypes-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("BetaStar");
						
					fileResultsUpBound2.addVector(UpBound2, UpBound2.size(), 0);
					fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"UpBound2-randomData-varyNumberTypes-var"+MainTests.variance);	
					fileOctave.writeOctaveMeans("UpBound2");
					
					//standard error
					File FicError = new File(fileErrorOpt);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> opts = fileResultsOpt.getVectors();
						String text = "";
						for (int k=0;k<currentNumTypes.length;k++){
							double average = 0;
							for (int l=0;l<opts.size();l++){
								average = average +opts.get(l).get(k);
							}
							average = average*1.0/opts.size();
							double err=0;
							for (int l=0;l<opts.size();l++){
								err = err +Math.pow(opts.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/opts.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					FicError = new File(fileErrorOptConst);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> optsConst = fileResultsOptConst.getVectors();
						String text = "";
						for (int k=0;k<currentNumTypes.length;k++){
							double average = 0;
							for (int l=0;l<optsConst.size();l++){
								average = average +optsConst.get(l).get(k);
							}
							average = average*1.0/optsConst.size();
							double err=0;
							for (int l=0;l<optsConst.size();l++){
								err = err +Math.pow(optsConst.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/optsConst.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					FicError = new File(fileErrorLottery);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> lottery = fileResultsLottery.getVectors();
						String text = "";
						for (int k=0;k<currentNumTypes.length;k++){
							double average = 0;
							for (int l=0;l<lottery.size();l++){
								average = average +lottery.get(l).get(k);
							}
							average = average*1.0/lottery.size();
							double err=0;
							for (int l=0;l<lottery.size();l++){
								err = err +Math.pow(lottery.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/lottery.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
						
					FicError = new File(fileErrorUpBound1);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> UpBound1 = fileResultsUpBound1.getVectors();
						String text = "";
						for (int k=0;k<currentNumTypes.length;k++){
							double average = 0;
							for (int l=0;l<UpBound1.size();l++){
								average = average +UpBound1.get(l).get(k);
							}
							average = average*1.0/UpBound1.size();
							double err=0;
							for (int l=0;l<UpBound1.size();l++){
								err = err +Math.pow(UpBound1.get(l).get(k)-average,2);
							}
							err= Math.sqrt(err/UpBound1.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}	
					
					FicError = new File(fileErrorUpBound2);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileErrorOpt);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> upBound2 = fileResultsUpBound2.getVectors();
						String text = "";
						for (int k=0;k<currentNumTypes.length;k++){
							double average = 0;
							for (int l=0;l<upBound2.size();l++){
								average = average +upBound2.get(l).get(k);
							}
							average = average*1.0/upBound2.size();
							double err=0;
							for (int l=0;l<upBound2.size();l++){
								err = err +Math.pow(upBound2.get(l).get(k)-average,2);;
							}
							err= Math.sqrt(err/upBound2.size());
							text = text+err+" ";
						}
						FichierOctave.writeBytes(text);
					}catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}	
				}
			}
		}
		
		
		
		
		if (chicago){
		
			//Choose the number of children and the number of runs:
			int numberOfChildren = 5000;
			int numTests = 100;
			
			//Choose how to vary the variance for the utility model (locationBased):
			double[] vars = {0,10,50,100};
			
			//data obtained from the Chicago school allocation website:
			MainTests.numTypes = 4;
			double[] quotas = {0.25,0.25,0.25,0.25};
			double[] actual = {0.271,0.275,0.236, 0.218};
			
			//read the school locations
			ProcessResultsFile fileSchoolsLocations = new ProcessResultsFile();
			fileSchoolsLocations.runParser("Data/chicago/schools_locations");
			MainTests.blockLocations = new ArrayList<ArrayList<Double>>();
			for (int i=0;i<37;i++){
				MainTests.blockLocations.add(fileSchoolsLocations.getVectors().get(i));
			}
			
			//read the number of students for each school
			ProcessResultsFile fileSchoolsStudents = new ProcessResultsFile();
			fileSchoolsStudents.runParser("Data/chicago/schools_students");
			ArrayList<ArrayList<Double>> numStudentsPerSchool = new ArrayList<ArrayList<Double>>();
			for (int i=0;i<37;i++){
				numStudentsPerSchool.add(fileSchoolsStudents.getVectors().get(i));
			}
			
			//Initialization:
			MainTests.numBlocks = MainTests.blockLocations.size();
			MainTests.numAptPerBlock = new int[MainTests.numBlocks];
			MainTests.numAgentPerType = new int[numTypes];
			int tmp = 0;
			for (int b=0;b<MainTests.numBlocks;b++){
				MainTests.numAptPerBlock[b] = numStudentsPerSchool.get(b).get(0).intValue()/9; //the average number of students (first year)
				tmp = tmp + MainTests.numAptPerBlock[b];
				int[] capacities_b = new int[MainTests.numTypes];
				for (int l=0;l<numTypes;l++){
					capacities_b[l]= (int)(MainTests.numAptPerBlock[b]*quotas[l]);
				}
				MainTests.capacities.add(capacities_b);
			}
			for (int l=0;l<MainTests.numTypes;l++){
				MainTests.numAgentPerType[l] = (int)(numberOfChildren*actual[l]); //
				MainTests.numAgents = MainTests.numAgents + MainTests.numAgentPerType[l];
			}	
			double minALPHA = 0.25;
			
			//collect tract locations:
			ArrayList<ArrayList<Double>> tractsLoc1 = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> tractsLoc2 = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> tractsLoc3 = new ArrayList<ArrayList<Double>>();
			ArrayList<ArrayList<Double>> tractsLoc4 = new ArrayList<ArrayList<Double>>();		
			for (int k=1;k<14;k++){
				ProcessResultsFile tractsLocations = new ProcessResultsFile();
				tractsLocations.runParser("Data/chicago/tracts_locations_"+k);// for each tract, we have a vector that gives the coordinates
				ArrayList<ArrayList<Double>> tractsLoc = tractsLocations.getVectors();
				
				ProcessResultsFile tractsTiers = new ProcessResultsFile();
				tractsTiers.runParser("Data/chicago/tracts_tier_student_"+k); //for each tract, we have a vector with two components (the corresponding tier, the number of students) 
				ArrayList<ArrayList<Double>> tractsTier = tractsTiers.getVectors();

				//Sort tracts with respect to tiers:
				for (int j=0;j<tractsLoc.size();j++){
					if (tractsTier.get(j).get(0).intValue()==1){
						tractsLoc1.add(tractsLoc.get(j));
					}
					if (tractsTier.get(j).get(0).intValue()==2){
						tractsLoc2.add(tractsLoc.get(j));
					}
					if (tractsTier.get(j).get(0).intValue()==3){
						tractsLoc3.add(tractsLoc.get(j));
					}
					if (tractsTier.get(j).get(0).intValue()==4){
						tractsLoc4.add(tractsLoc.get(j));
					}	
				}
			}
			
			//start running tests:
			for (int t=0;t<numTests;t++){				
				
				for (int v = 0;v< vars.length;v++){
				
					MainTests.variance = vars[v];
					
					//Create result files
					ProcessResultsFile fileResults = new ProcessResultsFile();
					fileResults.createResultsFile(directory+"distanceBased-var"+MainTests.variance);
					String fileError=directory+"error-distanceBased-var"+MainTests.variance+".txt";
					
					//Create the children
					MainTests.agents = new ArrayList<Agent>();					
					ArrayList<Double> results = new ArrayList<Double>();			
					for (int l=0;l<MainTests.numTypes;l++){
						for (int i=0;i<MainTests.numAgentPerType[l];i++){
							double[] bestPosition = null; //for every child of type l, the preferred location is randomly selected among the tracts of tier l
							if (l == 0){
								int pos = (int)(Math.random()*tractsLoc1.size());
								bestPosition = new double[]{tractsLoc1.get(pos).get(0).doubleValue(),tractsLoc1.get(pos).get(1).doubleValue()};
							}
							if (l == 1){
								int pos = (int)(Math.random()*tractsLoc2.size());
								bestPosition = new double[]{tractsLoc2.get(pos).get(0).doubleValue(),tractsLoc2.get(pos).get(1).doubleValue()};
							}
							if (l == 2){
								int pos = (int)(Math.random()*tractsLoc3.size());
								bestPosition = new double[]{tractsLoc3.get(pos).get(0).doubleValue(),tractsLoc3.get(pos).get(1).doubleValue()};
							}
							if (l == 3){
								int pos = (int)(Math.random()*tractsLoc4.size());
								bestPosition = new double[]{tractsLoc4.get(pos).get(0).doubleValue(),tractsLoc4.get(pos).get(1).doubleValue()};
							}
							Agent a = new Agent(l,bestPosition,MainTests.variance); // the utilities are defined using the normal distribution
							a.chooseTwentySchools(); //only 20 schools will have a positive utility because children are allowed to apply to at most 20 schools
							a.sortItems(); //sort schools according to utility
							MainTests.agents.add(a);
								
						}	
					}
					
					//OPT without diversity constraints
					double optWithout = Solver.assignmentWithOrWithoutConstraintsChicago(false);
					results.add(optWithout);
					
					//OPT with diversity constraints:
					double optConstraints = Solver.assignmentWithOrWithoutConstraintsChicago(true);
					results.add(optWithout/optConstraints);
			
					//average total utility for the lottery mechanism:
					double meanLot=0;
					int numLot = 100;
					for (int l=0;l<numLot;l++){
						meanLot = meanLot+ MainTests.lotteryProcess();
					}
					meanLot = meanLot*1.0/numLot;
					results.add(optWithout/meanLot);
						
					//UB1:
					results.add(1.0/minALPHA);
								
					double betaStar = Solver.computeDeltaStar(optWithout);	
					results.add(betaStar);
					
					//UB2:
					double sum = 0;
					for (int l=0;l<MainTests.numTypes;l++){
						double minl = Double.MAX_VALUE;
						for (int b=0;b<MainTests.numBlocks;b++){
							double alpha = MainTests.capacities.get(b)[l]*1.0/MainTests.numAptPerBlock[b];
							if (alpha < minl){
								minl = alpha;
							}
						}
						sum = sum + minl * MainTests.numAgentPerType[l]*1.0/MainTests.numAgents;
					}
					results.add(1.0/(sum*betaStar));
						
					//Write results: 
					fileResults.addVector(results, results.size(), 0);
					ProcessResultsFile fileOctave = new ProcessResultsFile();
					fileOctave.runParser(directory+"distanceBased-var"+MainTests.variance);
					fileOctave.writeOctaveMeans("totalUtilitiesBound1BetaBound2");
								
					//standard error:
					File FicError = new File(fileError);
					try {
						FicError.createNewFile();
						RandomAccessFile FichierOctave = new RandomAccessFile(FicError,"rw");
						if (FichierOctave.length()!=0){
							FicError.delete();
							FicError = new File(fileError);
							FichierOctave = new RandomAccessFile(FicError,"rw");
						}
						ArrayList<ArrayList<Double>> opts = fileResults.getVectors();
						double average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(1);
						}
						average = average*1.0/opts.size();
						double err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(1)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						String text = "constraint: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size())) +"\n" ;
									
						average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(2);
						}
						average = average*1.0/opts.size();
						err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(2)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						text = text + "lottery: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size())) +"\n";
											
						average = 0;
						for (int i=0;i<opts.size();i++){
							average = average +opts.get(i).get(5);
						}
						average = average*1.0/opts.size();
						err=0;
						for (int i=0;i<opts.size();i++){
							err = err +Math.pow(opts.get(i).get(5)-average,2);;
						}
						err= Math.sqrt(err/opts.size());
						text = text + "bound: sd = "+ err + ", se = "+(err/Math.sqrt(opts.size()));
							
						FichierOctave.writeBytes(text);
								
					} catch (IOException e) {
						// TODO Auto-generated catch block						
						e.printStackTrace();
					}
					
						
					
				}
			}
		}
			
			
	}
		
		

}