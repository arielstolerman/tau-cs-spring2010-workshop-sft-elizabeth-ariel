/*
 * 
 * Workshop Learning & Coding Theory Project - The SFT Algorithm
 * TAU, Spring Semester 2010
 * Elizabeth Firman and Ariel Stolerman
 * 
 * Filename: XMLFunction.java
 * Description: Code for class XMLFunction, implementation of the abstract class Function for XML inputs
 * 
 */

package Function;

import java.util.*;
import java.io.File;
import java.io.IOException;
import javax.xml.parsers.*;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import SFT.*;

/**
 * @author Elizabeth Firman and Ariel Stolerman
 * This class is an implementation of the abstract class Function used by class SFT, which defines a Function
 * object by a given XML file that contains description of functions over Z_N -> C by elements in Z_N and their
 * coefficients.
 * TODO: add XML input specification.
 */
public class XMLFourierPolynomial extends Function {
	
	private Map<String,FourierPolynomial> polynomials = null;
	private XMLParser parser;
	private boolean isRandom;
	private long[] maxAlpha;
	
	/**
	 * Constructs a XMLFunction object from an input XML file as described in the class documentation.
	 * @param XMLInputFile
	 * 				The input XML file.
	 * @param G
	 * 				The vector of values describing G, i.e. Cartesian product of Z_Ni.
	 * @throws FunctionException
	 * 				If an XML parsing error occurred, I/O exception or invalid input.
	 */
	public XMLFourierPolynomial(File XMLInputFile, long[] G) throws FunctionException{
		super(G);
		maxAlpha = new long[G.length];
		for (int i=0; i<G.length; i++) maxAlpha[i] = -1;
		parser = new XMLParser(XMLInputFile, polynomials, G); // will set terms according to the XML input
		polynomials = parser.polynomials;
		if (polynomials == null || polynomials.isEmpty())
			throw new FunctionException("The XML must contain at least one polynomial.");
		isRandom = parser.runId.equalsIgnoreCase("random");
		maxAlpha = parser.maxAlpha;
	}

	/**
	 * Returns the value of the function for the input element in G.
	 * If the Function is defined as random, a polynomial will randomly be selected from the list
	 * of Fourier polynomials and its value for the input element will be returned. Otherwise, the sole
	 * polynomial defined in the input XML file's value for the input element will be returned.
	 * @param elem
	 * 			The element whose this function's value is calculated.
	 * @return
	 * 			The value of the function for the input element in G.
	 */
	public Complex getValue(long[] elem){
		String randIndex = (int)Math.ceil(Math.random()*polynomials.size())+"";
		return polynomials.get(randIndex).getValue(elem);
	}
	
	/*
	 * getters
	 */
	
	/**
	 * Specifies weather this XMLFunction object is "random" as described in the class documentation.
	 * That is, is value calculation performed randomly over all the polynomials given in the input XML file
	 * or only a specific polynomial is specified.
	 * @return
	 * 			True if and only if the value calculation mode is random over all the polynomials in the input
	 * 			XML file.
	 */
	public boolean isRandom(){
		return isRandom;
	}
	
	/*
	 * setters
	 */
	
	@Override
	/**
	 * Sets G to a new vector of values.
	 * @param G
	 * 			The vector of values describing G, i.e. Cartesian product of Z_Ni, the domain of the function
	 * @throws FunctionException
	 * 			If one of the given values is less than or equals to 0 or less than or equals to the largest
	 * 			element whose term is specified in one or more of the input polynomials.
	 */
	public void setG(long[] G) throws FunctionException{
		if (G.length != this.G.length)
			throw new FunctionException("the given G is of wrong length.");
		
		// check the maxAlpha condition
		int k = G.length;
		for (int i=0; i<k; i++){
			if (maxAlpha[i] >= G[i])
				throw new
				FunctionException("The input polynomials contain elements greater than or equal to one of the given values. " +
				"Cannot change G.");
		}
		// call super that will check the greater than 0 condition
		super.setG(G);
	}
	
	/* *****************************
	 * private functions and classes
	 *******************************/
	
	/**
	 * Enum to indicate the current scope of the XML file
	 */
	private enum Tag{
		FUNCTIONS,FUNCTION,TERM,ALPHA,COORD,RECOEFF,IMCOEFF,END;
	}
	
	/**
	 * The private XMLParser class is used for parsing an XML input file as described in the XMLFourierPolynomial
	 * class documentation, and constructing a list of polynomials to be used on the base XMLFourierPolynomial
	 * 'getValue()' method.
	 */
	private class XMLParser extends DefaultHandler{
		/*
		 * members
		 */
		
		// parsing state indicators
		private Tag currTag;
		private String runId;
		private String funcId = null;
		private double recoeff, imcoeff;
		private long[] alpha = null;
		private long coord = -1;
		private int coordIndex = -1;
		
		// XMLFourierFunction inputs
		private long[] maxAlpha;
		private long[] G;
		private Map<String,FourierPolynomial> polynomials;
		
		public XMLParser(File XMLInputFile, Map<String,FourierPolynomial> polynomials, long[] G) throws FunctionException{
			this.polynomials = polynomials;
			this.G = G;
			int k = G.length;
			maxAlpha = new long[k];
			for (int i=0; i<k; i++){
				maxAlpha[i] = -1;
			}
			parseDocument(XMLInputFile);
		}
		
		/**
		 * Main parsing procedure
		 * parses the input XML file
		 * @throws FunctionException
		 * 			If an error occurred during parsing, including XML parsing error or I/O exception 
		 */
		public void parseDocument(File XMLInputFile) throws FunctionException {
			
			Debug.log("XMLParser -> parseDocument started");
			
			try {
				//get a factory
				SAXParserFactory spf = SAXParserFactory.newInstance();
				//get a new instance of parser
				SAXParser sp = spf.newSAXParser();
				//parse the file and also register this class for call backs
				sp.parse(XMLInputFile, this);
			}catch(SAXException se) {
				throw new FunctionException("SAXException had occurred.\n"+se.getMessage());
			}catch(ParserConfigurationException pce) {
				throw new FunctionException("ParserConfigurationException had occurred.\n"+pce.getMessage());
			}catch (IOException ie) {
				throw new FunctionException("IOException had occurred.\n"+ie.getMessage());
			}
			
			Debug.log("XMLParser -> parseDocument completed");
		}
		
		/* ***************
		 *  Event handlers
		 *****************/
		
		/**
		 * handlers for start tags
		 */
		public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
			
			if(qName.equalsIgnoreCase("functions")) {
				// set runid
				runId = attributes.getValue(0);
				// check if valid, else throw exception
				try{
					Integer.parseInt(runId);
				} catch (NumberFormatException nfe){
					if (!runId.equalsIgnoreCase("random"))
						throw new SAXException("functions opening tag has invalid runid (not \"random\" or an integer): "+runId);
				}
				
				// initialize polynomials map
				polynomials = new HashMap<String,FourierPolynomial>();
				// set current tag
				currTag = Tag.FUNCTIONS;
				
				Debug.log("<functions runid=\""+runId+"\">");
			}
			
			else if(qName.equalsIgnoreCase("function")){
				String funcId = attributes.getValue(0);
				// check if valid, else throw exception
				try{
					Integer.parseInt(funcId);
				} catch (NumberFormatException nfe){
					throw new SAXException("function opening tag has invalid id (not an integer): "+funcId);
				}
				// if this function is the one we're looking for, or we're in random mode, get it
				if (runId.equalsIgnoreCase("random") || runId.equals(funcId)){
					try{
						FourierPolynomial poly = new FourierPolynomial(G,funcId);
						// add to the polynomials map
						polynomials.put(funcId, poly);
						// set the current function that is parsed to this one
						this.funcId = funcId;
					}catch(FunctionException fe){
						throw new SAXException("a function exception had occurred during parsing for function with id "+funcId);
					}
				}
				// set current tag
				currTag = Tag.FUNCTION;
				
				Debug.log("\t<function id=\""+funcId+"\">");
			}
			else if(qName.equalsIgnoreCase("term")){
				// do nothing
				// set current tag
				currTag = Tag.TERM;
				
				Debug.log("\t\t<term>");
			}
			else if(qName.equalsIgnoreCase("alpha")){
				// set current tag
				currTag = Tag.ALPHA;
				// initialize alpha
				alpha = new long[G.length];
				
				Debug.log("\t\t\t<alpha>");
			}
			else if(qName.equalsIgnoreCase("coord")){
				// set current tag
				currTag = Tag.COORD;
				
				// set coordIndex
				String ind = attributes.getValue(0);
				// check if valid, else throw exception
				try{
					coordIndex = Integer.parseInt(ind);
					if (coordIndex < 0 || coordIndex >= G.length)
						throw new SAXException("coord opening tag has invalid index, must be between 0 and "+G.length);
				} catch (NumberFormatException nfe){
					throw new SAXException("coord opening tag has invalid index (not an integer): "+ind);
				}
				
				Debug.log("\t\t\t\t<coord index=\""+coordIndex+"\">");
			}
			else if(qName.equalsIgnoreCase("reCoeff")){
				// set current tag
				currTag = Tag.RECOEFF;
				
				Debug.log("\t\t\t<reCoeff>");
			}
			else if(qName.equalsIgnoreCase("imCoeff")){
				// set current tag
				currTag = Tag.IMCOEFF;
				
				Debug.log("\t\t\t<imCoeff>");
			}
			// otherwise
			else {
				throw new SAXException("unrecognized opening tag: "+qName);
			}
		}
		
		/**
		 * handlers for end tags
		 */
		public void endElement(String uri, String localName,String qName) throws SAXException {
			int k = G.length;
			if(qName.equalsIgnoreCase("functions")) {
				// set current tag
				currTag = Tag.END;
				
				Debug.log("</functions>");
			}
			else if(qName.equalsIgnoreCase("function")){
				funcId = null;
				// set current tag
				currTag = Tag.FUNCTIONS;
				
				Debug.log("\t</function>");
			}
			else if(qName.equalsIgnoreCase("term")){
				// add a new term to the current polynomial (only if needed)
				FourierPolynomial p = polynomials.get(funcId);
				if (p != null) p.addUpdateTerm(alpha, recoeff, imcoeff);
				
				// update the maximum elements for each N in G seen so far
				for (int i=0; i<k; i++){
					if (alpha[i] > maxAlpha[i]) maxAlpha[i] = alpha[i];
				}
				
				// set current tag
				currTag = Tag.FUNCTION;
				
				Debug.log("\t\t</term>");
			}
			else if(qName.equalsIgnoreCase("alpha")){
				// set current tag
				currTag = Tag.TERM;
				
				// check that all coordinates were entered
				for (int i=0; i<k; i++){
					if (alpha[i] == -1) throw new SAXException("one of the elements for function "+funcId+" is missing value "+
							"for index "+i);
				}
				
				Debug.log("\t\t\t</alpha>");
			}
			else if(qName.equalsIgnoreCase("coord")){
				// set current tag
				currTag = Tag.ALPHA;
				// insert the coordinate into alpha (may override value)
				alpha[coordIndex] = coord;
				
				Debug.log("\t\t\t\t</coord>");
			}
			else if(qName.equalsIgnoreCase("reCoeff")){
				// set current tag
				currTag = Tag.TERM;
				
				Debug.log("\t\t\t</reCoeff>");
			}
			else if(qName.equalsIgnoreCase("imCoeff")){
				// set current tag
				currTag = Tag.TERM;
				
				Debug.log("\t\t\t</imCoeff>");
			}
			// otherwise
			else {
				throw new SAXException("unrecognized closing tag: "+qName);
			}
		}
		
		/**
		 * handler for text between open and close tags
		 */
		public void characters(char ch[], int start, int length) throws SAXException {
			String str = new String(ch, start, length);
			
			switch (currTag){
			case FUNCTIONS:
				// do nothing
				break;
			case FUNCTION:
				// do nothing
				break;
			case TERM:
				// do nothing
				break;
			case ALPHA:
				// do nothing
				break;
			case COORD:
				try{
					Long value = Long.parseLong(str);
					coord = value;
					Debug.log("\t\t\t\t\t"+str);
					if (coord < 0 || coord >= G[coordIndex])
						throw new SAXException("coordinate in index "+coordIndex+" not in range [0,1,...,G["+coordIndex+"]-1]: "+coord);
				} catch (NumberFormatException nfe){
					throw new SAXException("coordinate must be a number in range [0,1,...,G["+coordIndex+"]-1]: "+coord);
				}
				break;
			case RECOEFF:
				try{
					recoeff = Double.parseDouble(str);
					Debug.log("\t\t\t\t"+str);
				} catch (NumberFormatException nfe){
					throw new SAXException("reCoeff not a double");
				}
				break;
			case IMCOEFF:
				try{
					imcoeff = Double.parseDouble(str);
					Debug.log("\t\t\t\t"+str);
				} catch (NumberFormatException nfe){
					throw new SAXException("imCoeff not a double");
				}
				break;
			case END:
				throw new SAXException("XML parsing error");
			default:
				// do nothing
			}
		}
	}

	public Map<String, FourierPolynomial> getPolynomials() {
		return polynomials;
	}
	
}
