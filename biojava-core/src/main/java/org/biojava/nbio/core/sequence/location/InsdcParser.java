/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 01-21-2010
 */
package org.biojava.nbio.core.sequence.location;

import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DataSource;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.location.template.AbstractLocation;
import org.biojava.nbio.core.sequence.location.template.Location;
import org.biojava.nbio.core.sequence.location.template.Point;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Parser for working with INSDC style locations. This class supports the
 * full range of location types generated by Genbank, INSDC and ENA.
 *
 * @author ayates
 * @author jgrzebyta
 * @author Paolo Pavan
 */
public class InsdcParser {

	private boolean isSequenceCircular;
	private long sequenceLength;

	private final DataSource dataSource;

		/**
	 * parse a location. if group(1) is null than the feature is on the positive
	 * strand, group(2) start position, group(3) end position.
	 */
	protected static final Pattern singleLocationPattern = Pattern.compile("(?:([A-Za-z\\.0-9_]*?):)?(<?)(\\d+)(\\.{2}|\\^)?(>?)(\\d+)?(>?)?");
	/**
	 * Decodes a split pattern. Split patterns are a composition of multiple
	 * locationsString qualified by actions: join(location,location, ...
	 * location): The indicated elements should be joined (placed end-to-end) to
	 * form one contiguous sequence. order(location,location, ... location): The
	 * elements can be found in the specified order (5' to 3' direction),
	 * nothing is implied about their reasonableness
	 * bond(location,location...location): Found in protein files. These
	 * generally are used to describe disulfide bonds.
	 * complement(location,location...location): consider locations in their
	 * complement versus
	 *
	 * takes in input a comma splitted location string. The split must be done
	 * for outer level commas group(1) is the qualifier group(2) is the location
	 * string to getFeatures. In case of complex splits it will contain the
	 * nested expression
	 *
	 * Not really sure that they are not declared obsolete but they are still in
	 * several files.
	 */
	protected static final Pattern genbankSplitPattern = Pattern.compile("^\\s?(join|order|bond|complement|)\\(?(.+)\\)?");
	/**
	 * designed to recursively split a location string in tokens. Valid tokens
	 * are those divided by coma that are not inside a bracket. I. e. split on
	 * the comma only if that comma has zero, or an even number of quotes ahead
	 * of it.
	 */
	protected static final String locationSplitPattern = ",(?=([^\\(|\\)]+\\([^\\(|\\)]+\\))[^\\(|\\)]+)";
	/**
	 * these variables are used to compute the global start and end of complex
	 * features
	 */
	protected Integer featureGlobalStart, featureGlobalEnd;

	enum complexFeaturesAppendEnum {

		FLATTEN, HIERARCHICAL;
	}
	/**
	 * define the mode in which complex features should be appended in FLATTEN
	 * mode their will be appended to a single feature in HIERARCHICAL mode, the
	 * single mother feature will have a tree of features that will reflect the
	 * construction in genbank file
	 */
	private complexFeaturesAppendEnum complexFeaturesAppendMode = complexFeaturesAppendEnum.HIERARCHICAL;

	public void setComplexFeaturesAppendMode(complexFeaturesAppendEnum complexFeaturesAppendMode) {
		this.complexFeaturesAppendMode = complexFeaturesAppendMode;
	}

	public InsdcParser() {
		this(DataSource.ENA);
	}

	public InsdcParser(DataSource dataSource) {
		this.dataSource = dataSource;
	}

	public DataSource getDataSource() {
		return dataSource;
	}

	public void setSequenceCircular(boolean sequenceCircular) {
		isSequenceCircular = sequenceCircular;
	}

	public void setSequenceLength(long sequenceLength) {
		this.sequenceLength = sequenceLength;
	}

	/**
	 * Main method for parsing a location from a String instance
	 *
	 * @param locationString Represents a logical location
	 * @return The parsed location
	 * @throws ParserException thrown in the event of any error during parsing
	 */
	public Location parse(String locationString) throws ParserException {
		featureGlobalStart = Integer.MAX_VALUE;
		featureGlobalEnd = 1;

		Location l;
		List<Location> ll = parseLocationString(locationString, 1);

		if (ll.size() == 1) {
			l = ll.get(0);
		} else {
			l = new SimpleLocation(
					new SimplePoint(featureGlobalStart),
					new SimplePoint(featureGlobalEnd),
					Strand.UNDEFINED,
					isSequenceCircular,
					ll);
		}
		return l;
	}

	private List<Location> parseLocationString(String string, int versus) throws ParserException {
		Matcher m;
		List<Location> boundedLocationsCollection = new ArrayList<Location>();

		List<String> tokens = splitString(string);
		for (String t : tokens) {
			m = genbankSplitPattern.matcher(t);
			if (!m.find()) {
				throw new ParserException("Cannot interpret split pattern " + t
						+ "\nin location string:" + string);
			}
			String splitQualifier = m.group(1);
			String splitString = m.group(2);

			if (!splitQualifier.isEmpty()) {
				//recursive case
				int localVersus = splitQualifier.equalsIgnoreCase("complement") ? -1 : 1;
				List<Location> subLocations = parseLocationString(
						splitString, versus * localVersus);

				switch (complexFeaturesAppendMode) {
					case FLATTEN:
						boundedLocationsCollection.addAll(subLocations);
						break;
					case HIERARCHICAL:
						if (subLocations.size() == 1) {
							boundedLocationsCollection.addAll(subLocations);
						} else {
							Point min = Location.Tools.getMin(subLocations).getStart();
							Point max = Location.Tools.getMax(subLocations).getEnd();
							AbstractLocation motherLocation
									= new SimpleLocation(
											min,
											max
									);

							if (splitQualifier.equalsIgnoreCase("join")) {
								motherLocation = new InsdcLocations.GroupLocation(subLocations);
							}
							if (splitQualifier.equalsIgnoreCase("order")) {
								motherLocation = new InsdcLocations.OrderLocation(subLocations);
							}
							if (splitQualifier.equalsIgnoreCase("bond")) {
								motherLocation = new InsdcLocations.BondLocation(subLocations);
							}
							motherLocation.setStrand(getGroupLocationStrand(subLocations));
							boundedLocationsCollection.add(motherLocation);
						}
					break;
				}
			} else {
				//base case
				m = singleLocationPattern.matcher(splitString);
				if (!m.find()) {
					throw new ParserException("Cannot interpret location pattern " + splitString
							+ "\nin location string:" + string);
				}

				String accession = m.group(1);
				Strand s = versus == 1 ? Strand.POSITIVE : Strand.NEGATIVE;
				int start = Integer.valueOf(m.group(3));
				int end = m.group(6) == null ? start : Integer.valueOf(m.group(6));

				if (featureGlobalStart > start) {
					featureGlobalStart = start;
				}
				if (featureGlobalEnd < end) {
					featureGlobalEnd = end;
				}

				AbstractLocation l;
				if (start <= end) {
					l = new SimpleLocation(
							start,
							end,
							s
					);
				} else {
					// in case of location spanning the end point, Location contract wants sublocations
					AbstractLocation l5prime = new SimpleLocation(
							1,
							end,
							Strand.UNDEFINED
							);
					AbstractLocation l3prime = new SimpleLocation(
							start,
							(int) sequenceLength,
							Strand.UNDEFINED
							);

					l = new InsdcLocations.GroupLocation(
							new SimplePoint(start),
							new SimplePoint(end),
							s,
							isSequenceCircular,
							l5prime, l3prime
					);

				}

				if(m.group(4) != null && m.group(4).equals("^")) l.setBetweenCompounds(true);

				if (m.group(2).equals("<")) {
					l.setPartialOn5prime(true);
				}
				if (m.group(5) != null && (m.group(5).equals(">") || m.group(7).equals(">"))) {
					l.setPartialOn3prime(true);
				}

				if (!(accession == null || "".equals(accession))) l.setAccession(new AccessionID(accession));

				boundedLocationsCollection.add(l);

			}
		}

		return boundedLocationsCollection;
	}


	private List<String> splitString(String input) {
		List<String> result = new ArrayList<String>();
		int start = 0;
		int openedParenthesis = 0;
		for (int current = 0; current < input.length(); current++) {
			if (input.charAt(current) == '(') {
				openedParenthesis++;
			}
			if (input.charAt(current) == ')') {
				openedParenthesis--;
			}
			boolean atLastChar = (current == input.length() - 1);
			if (atLastChar) {
				result.add(input.substring(start));
			} else if (input.charAt(current) == ',' && openedParenthesis == 0) {
				result.add(input.substring(start, current));
				start = current + 1;
			}
		}
		return result;
	}

	private Strand getGroupLocationStrand(List<Location> ll){
		Strand returnStrand = null;

		for (Location l: ll) {
			if (returnStrand == null) returnStrand = l.getStrand();
			if (returnStrand != l.getStrand()) return Strand.UNDEFINED;
		}
		return returnStrand;
	}

	public static void main(String[] args){
		String[] testStrings = {
			"J00194.1:100..202",
			"A00001.5:34..45",
			"43..129",
			"bond(55,110)",
			"bond(34,35),join(56..80),complement(45,73)",
			"order(complement(30,40),70..80),bond(34,35),join(56,80),complement(45..56)",
			"join(join(complement(30,40),complement(70..80)),bond(34,35),join(56,80),complement(45..56))",
			"complement(join(complement(2000..4000),complement(70..80)),bond(34,35),join(56,80),complement(45..56))",

		};
		InsdcParser p = new InsdcParser();
		p.setComplexFeaturesAppendMode(complexFeaturesAppendEnum.HIERARCHICAL);

		for (String s: testStrings){
			Location l = p.parse(s);
			System.out.println(l.toString());
		}

	}

}
