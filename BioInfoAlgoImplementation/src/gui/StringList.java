package gui;

/**
 * Dummy class to report a multiple line String parameter to the user interface.
 * 
 * @author Martin Mann - http://www.bioinf.uni-freiburg.de/~mmann/
 *
 */
public class StringList {
	
	 /**
	  * Holds the content of the multiple line String (including the newline
	  * characters).
	  */
	private String value = null;
	
	
	/**
	 * Default construction of a dummy object.
	 */
	public StringList() {
		this.value = null;
	}

	 /**
	  * Constructs a dummy object.
	  * @param value the initial value to hold
	  */
	public StringList(String value) {
		this.value = value;
	}
	
	 /**
	  * Access to the String value of this object.
	  * @return the multiple line String object.
	  */
	public String toString() {
		return value;
	}

}
