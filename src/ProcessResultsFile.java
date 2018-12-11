import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;


import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

import java.util.ArrayList;
import java.util.StringTokenizer;
import java.lang.Integer;

public class ProcessResultsFile extends DefaultHandler {
	
	private String ficXML;
	private ArrayList<ArrayList<Double>> vectors = new ArrayList<ArrayList<Double>>();;
	private double[] tmpVector;
	private int numVectors = 0;
	private int indexValue;
	private String tmpVal;
	private ArrayList<Double> meanResults = null;

	public void createResultsFile(String ficXML){
		this.ficXML = ficXML+".txt";
		File FicXML = new File(this.ficXML);
		String text = "";
		
		try {
			FicXML.createNewFile();
			RandomAccessFile Fichier = new RandomAccessFile(FicXML,"rw");
			if (Fichier.length()== 0){
				text = text + "<?xml version=\"1.0\" encoding=\"utf-8\" standalone=\"no\"?>\n<!DOCTYPE results SYSTEM \"../results.dtd\">\n\n<results>\n</results>";
				Fichier.writeBytes(text);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
	}
	
	public void deleteResultsFile(){
		File FicXML = new File(this.ficXML);
		FicXML.delete();
	}
	
	public void deleteResultsFile(String fic){
		File FicXML = new File(fic);
		FicXML.delete();
	}

	
	public void runParser(String ficXML) {
		this.ficXML = ficXML+".txt";
		parseDocument(); 
	}

	public void parseDocument() {
		
		SAXParserFactory factory = SAXParserFactory.newInstance();
		try {
		
			//On obtient une instance de parser :
			SAXParser parser = factory.newSAXParser();
			factory.setValidating(true);
            factory.setNamespaceAware(true);
			//lancement de la lecture :
            System.out.println(this.ficXML);
			parser.parse(new File(this.ficXML), this);
			
		}catch(SAXException se) {
			System.out.println("probleme1");
			se.printStackTrace();
		}catch(ParserConfigurationException pce) {
			System.out.println("probleme2");
			pce.printStackTrace();
		}catch (IOException ie) {
			System.out.println("probleme3");
			ie.printStackTrace();
		}
		
	}
	
	public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
	//Réinitialisation
		this.tmpVal = "";
		if(qName.equalsIgnoreCase("vector")) {
			this.vectors.add(new ArrayList<Double>());
		}
	}
	
	public void characters(char[] ch, int start, int length) throws SAXException {
		this.tmpVal = new String(ch,start,length);
		//System.out.println(this.tmpVal);
		//System.out.println(this.numVectors);
	}
	
	public void endElement(String uri, String localName, String qName) throws SAXException {

		if(qName.equalsIgnoreCase("vector")) {
			//On ajoute la matrice à la liste "matrices"
			this.numVectors = this.numVectors+1;
			
		}else if (qName.equalsIgnoreCase("value")){
			this.vectors.get(this.numVectors).add(Double.parseDouble(this.tmpVal));
			this.indexValue = this.indexValue+1;
		}
	}
	
	public void addVector(ArrayList<Double> vector, int size, double numStepsToFind){
		this.numVectors = this.numVectors+1;
		this.vectors.add(vector);
		String chaine = "\n\t<vector SizeVector=\""+ size+"\" realCost=\""+numStepsToFind+"\">";
		FileWriter myWriter = null;
		int offset;
		System.out.println(this.ficXML);
		try{
			RandomAccessFile Fichier = new RandomAccessFile(this.ficXML,"rw");
			offset = (int) Fichier.length() - "</results>".length() - 1;
			Fichier.seek(offset);
			for(int i=0;i<size;i++)
				chaine = chaine + "\n\t\t<value>" + vector.get(i) + "</value>";
			chaine = chaine + "\n\t</vector>\n</results>";
			Fichier.writeBytes(chaine);

		}catch(IOException ex){
		    ex.printStackTrace();
		}finally{
		  if(myWriter != null){
		     try {
				myWriter.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		  }
		}
	}
	
	
	public void meanResults(double maxTaille){
		this.meanResults = new ArrayList<Double>();
		for (int i=0;i<maxTaille;i++){
			double mean = 0;
			for (int j=0;j<this.numVectors;j++){
				if (i<this.vectors.get(j).size()){
					mean = mean + this.vectors.get(j).get(i);
				}
			}
			this.meanResults.add(mean*1.0/this.numVectors);	
		}
	}
	
	
	
	public void writeOctaveMeans(String name){
		String ficOctave = this.ficXML+"-octave.txt";
		File FicOctave = new File(ficOctave);
		
		try {
			FicOctave.createNewFile();
			RandomAccessFile FichierOctave = new RandomAccessFile(FicOctave,"rw");
			if (FichierOctave.length()!=0){
				FicOctave.delete();
				FicOctave = new File(ficOctave);
				FichierOctave = new RandomAccessFile(FicOctave,"rw");
			}
			double maxTaille = this.vectors.get(0).size();
			for (int i=1;i<this.numVectors;i++){
				if (maxTaille<this.vectors.get(i).size())
					maxTaille = this.vectors.get(i).size();
			}
			String text = "# name: "+name+"\n# type : matrix\n# rows : 1\n# columns : "+ maxTaille+"\n";
			this.meanResults(maxTaille);
			for(int i=0;i<this.meanResults.size();i++){
				float t = this.meanResults.get(i).floatValue();
				text = text + t +" ";
			}
			FichierOctave.writeBytes(text);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
	}
	
	
	

	
	
	public ArrayList<Double> getMeanResults(){
		return this.meanResults;
	}
	
	public ArrayList<ArrayList<Double>> getVectors(){
		return this.vectors;
	}
}