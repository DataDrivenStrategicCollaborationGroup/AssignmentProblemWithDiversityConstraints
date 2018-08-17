import java.util.ArrayList;
import java.util.Random;

public class Agent {

	private ArrayList<double[]> utilities;
	private ArrayList<int[]> sortedItems = null;
	private int type;
	
	public Agent(int type){//random normalized utility
		
		this.type = type;
		this.utilities = new ArrayList<double[]>();
		double sum = 0;
		for (int b=0;b<MainTests.numBlocks;b++){
			int numApt_b = MainTests.numAptPerBlock[b];
			double[] utilities_b = new double[numApt_b];
			for (int j=0;j<numApt_b;j++){
				utilities_b[j] = Math.random(); 
				sum = sum + utilities_b[j];
			}
			utilities.add(utilities_b);
		}
		//normalization:
		for (int b=0;b<MainTests.numBlocks;b++){
			int numApt_b = MainTests.numAptPerBlock[b];
			double[] utilities_b = utilities.get(b);
			for (int j=0;j<numApt_b;j++){
				utilities_b[j] = utilities_b[j]*1.0/sum; 
			}
		}
	}
	
	public Agent(int type, double[] preferredLocation){//location based utility
		this.type = type;
		this.utilities = new ArrayList<double[]>();
		double sum = 0;
		for (int b=0;b<MainTests.numBlocks;b++){
			ArrayList<Double> blockLocation = MainTests.blockLocations.get(b);
			Random r = new Random();
			double mean = 1.0/Math.sqrt(Math.pow(preferredLocation[0]-blockLocation.get(0) , 2) + Math.pow(preferredLocation[1]-blockLocation.get(1) , 2));
			double var = MainTests.variance;
			int numApt_b = MainTests.numAptPerBlock[b];
			double[] utilities_b = new double[numApt_b];
			for (int j=0;j<numApt_b;j++){
				utilities_b[j] = mean + Math.max(0, r.nextGaussian()*var); 
				sum = sum + utilities_b[j];
			}
			this.utilities.add(utilities_b);
		}
		//normalization:
		for (int b=0;b<MainTests.numBlocks;b++){
			int numApt_b = MainTests.numAptPerBlock[b];
			double[] utilities_b = utilities.get(b);
			for (int j=0;j<numApt_b;j++){
				utilities_b[j] = utilities_b[j]*1.0/sum; 
			}
		}
	}
	
	public Agent(int type, double[] preferredLocation, double var){//distance based utility - chicago
		this.type = type;
		this.utilities = new ArrayList<double[]>();
		double sum = 0;
		for (int b=0;b<MainTests.numBlocks;b++){
			ArrayList<Double> blockLocation = MainTests.blockLocations.get(b);
			Random r = new Random();
			double val = 1.0/Math.sqrt(Math.pow(preferredLocation[0]-blockLocation.get(0) , 2) + Math.pow(preferredLocation[1]-blockLocation.get(1) , 2)) +Math.max(0, r.nextGaussian()*var);  
			int numApt_b = MainTests.numAptPerBlock[b];
			double[] utilities_b = new double[numApt_b];
			for (int j=0;j<numApt_b;j++){
				utilities_b[j] = val; 
				sum = sum + utilities_b[j];
			}
			this.utilities.add(utilities_b);
		}
		
	}
	
	public Agent(int type, double income, ArrayList<ArrayList<Double>> prices){//price based utilities
		
		this.type = type;
		this.utilities = new ArrayList<double[]>();
		double sum = 0;
		for (int b=0;b<MainTests.numBlocks;b++){
			//System.out.println("block : " + b);
			double[] utilities_b = new double[MainTests.numAptPerBlock[b]];
			ArrayList<Double> price_b = prices.get(b);
			for (int t=0;t<MainTests.numAptPerBlock[b];t++){
				//System.out.println("flat : " + t);
				utilities_b[t] = 1.0/Math.pow(income*1.0/3-price_b.get(t),2);
				sum = sum + utilities_b[t];
			}
			utilities.add(utilities_b);
		}
		//normalization:
		for (int b=0;b<MainTests.numBlocks;b++){
			int numApt_b = MainTests.numAptPerBlock[b];
			double[] utilities_b = utilities.get(b);
			for (int j=0;j<numApt_b;j++){
				utilities_b[j] = utilities_b[j]*1.0/sum; 
			}
		}
	}	

	public void sortItems(){//sort item by preference in each block
		this.sortedItems = new ArrayList<int[]>();
		for (int b=0;b<MainTests.numBlocks;b++){
			int numApt_b = MainTests.numAptPerBlock[b];
			int[] sortedItemsb = new int[numApt_b];
			double[] sortedValues = new double[numApt_b];
			double[] utilities_b = utilities.get(b);
			sortedItemsb[0] = 0;
			for (int j=0;j<numApt_b;j++){
				sortedValues[j] = utilities_b[j];
			}
			for (int j=1;j<numApt_b;j++){
				double val = sortedValues[j];
				int pos = j;
				while (pos>0){
					double valPosMoins1 = sortedValues[pos-1];
					if (val > valPosMoins1){
						sortedValues[pos] = valPosMoins1; 
						sortedItemsb[pos] = sortedItemsb[pos-1];
						pos = pos-1;
					}else{
						break;
					}
				}
				sortedValues[pos] = val;
				sortedItemsb[pos] = j;
			}
			this.sortedItems.add(sortedItemsb);
		}
	}
	
	public void chooseTwentySchools(){// for chicago public school assignment problem
		double[] utilityPerSchool = new double[MainTests.numBlocks];
		int[] indexSortedSchool = new int[MainTests.numBlocks];
		for (int b=0;b<MainTests.numBlocks;b++){
			utilityPerSchool[b] = this.utilities.get(0)[0];
			indexSortedSchool[b] = b;
		}
		for (int b=1;b<MainTests.numBlocks;b++){
			double val = utilityPerSchool[b];
			int pos = b;
			while (pos>0){
				double valPosMoins1 = utilityPerSchool[pos-1];
				if (val > valPosMoins1){
					utilityPerSchool[pos] = valPosMoins1; 
					indexSortedSchool[pos] = indexSortedSchool[pos-1];
					pos = pos-1;
				}else{
					break;
				}
			}
		}
		double sum = 0;
		for (int b=0;b<20;b++){
			sum = sum + this.utilities.get(indexSortedSchool[b])[0]*MainTests.numAptPerBlock[indexSortedSchool[b]];
		}
		for (int b=0;b<20;b++){
			double[] utilities_b = this.utilities.get(indexSortedSchool[b]);
			for (int i=0;i<utilities_b.length;i++){
				utilities_b[i]= utilities_b[i]*1.0/sum;
			}
		}
		for (int b=20;b<MainTests.numBlocks;b++){
			double[] utilities_b = this.utilities.get(indexSortedSchool[b]);
			for (int i=0;i<utilities_b.length;i++){
				utilities_b[i] = -1;
			}
		}
	}
	
	public ArrayList<double[]> getUtilities(){
		return this.utilities;
	}
	
	public ArrayList<int[]> getSortedItems(){
		return this.sortedItems;
	}
	
	public int getType(){
		return this.type;
	}
}
