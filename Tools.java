
package Berlekamp_Welch_Algorithm.test;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.util.Scanner;
import java.io.IOException;
import java.util.ArrayList;
import java.util.NoSuchElementException;
import java.math.BigInteger;
import java.math.BigDecimal;

public class Tools {
	
	public Tools() {
		
	}//constructor
	
	//Can only do maximum 11 bits long
	public int ConvertBinaryStringToPositiveInteger(String binaryString) {
		if (binaryString.length() > 11)
			return Integer.MIN_VALUE;
		
		int[] rev_2_pow = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024 };
		int   answer = 0;
		int   from_the_end = binaryString.length() - 1;
		for (int ii=from_the_end; ii >= 0; ii--) {
			char mychar = binaryString.charAt(ii);
			if (mychar == '1')
				answer += rev_2_pow[from_the_end - ii];
		}
        return answer;		
	}//ConvertBinaryStringToPositiveInteger
	

	/************************************************************************************************
	QR Code Log Antilog Table for Galois Field 256
	The following table contains the log and antilog values that are used in GF(256) arithmetic, 
	which is used for generating the error correction codes required for QR codes.
	*************************************************************************************************/
	
	private int[] Exponent_Of_Alpha_To_Integer = { 1,     2,   4,   8,  16,  32,  64, 128,  29,  58, 116, 232, 205, 135,  19,  38,
	                                             76,  152,  45,  90, 180, 117, 234, 201, 143,   3,   6,  12,  24,  48,  96, 192,
	                                             157,  39,  78, 156,  37,  74, 148,  53, 106, 212, 181, 119, 238, 193, 159,  35,
	                                             70,  140,   5,  10,  20,  40,  80, 160,  93, 186, 105, 210, 185, 111, 222, 161,
	                                             95,  190,  97, 194, 153,  47,  94, 188, 101, 202, 137,  15,  30,  60, 120, 240,
	                                             253, 231, 211, 187, 107, 214, 177, 127, 254, 225, 223, 163,  91, 182, 113, 226,
	                                             217, 175,  67, 134,  17,  34,  68, 136,  13,  26,  52, 104, 208, 189, 103, 206,
	                                             129,  31,  62, 124, 248, 237, 199, 147,  59, 118, 236, 197, 151,  51, 102, 204,
	                                             133,  23,  46,  92, 184, 109, 218, 169,  79, 158,  33,  66, 132,  21,  42,  84,
	                                             168,  77, 154,  41,  82, 164,  85, 170,  73, 146,  57, 114, 228, 213, 183, 115,
	                                             230, 209, 191,  99, 198, 145,  63, 126, 252, 229, 215, 179, 123, 246, 241, 255,
	                                             227, 219, 171,  75, 150,  49,  98, 196, 149,  55, 110, 220, 165,  87, 174,  65,
	                                             130,  25,  50, 100, 200, 141,   7,  14,  28,  56, 112, 224, 221, 167,  83, 166,
	                                             81,  162,  89, 178, 121, 242, 249, 239, 195, 155,  43,  86, 172,  69, 138,   9,
	                                             18,   36,  72, 144,  61, 122, 244, 245, 247, 243, 251, 235, 203, 139,  11,  22,
	                                             44,   88, 176, 125, 250, 233, 207, 131,  27,  54, 108, 216, 173,  71, 142,   1  };
	
	public int[] Table_Exponent_Of_Alpha_To_Integer()
	{
		return Exponent_Of_Alpha_To_Integer;
	}
	
	private int[] Integer_To_Exponent_Of_Alpha = {   Integer.MAX_VALUE, 0, 1, 25, 2, 50, 26, 198, 3, 223, 51, 238, 27, 104, 199, 75,
	                                               4, 100, 224, 14, 52, 141, 239, 129, 28, 193, 105, 248, 200, 8, 76, 113,
	                                               5, 138, 101, 47, 225, 36, 15, 33, 53, 147, 142, 218, 240, 18, 130, 69,
	                                              29, 181, 194, 125, 106, 39, 249, 185, 201, 154, 9, 120, 77, 228, 114, 166,
	                                               6, 191, 139, 98, 102, 221, 48, 253, 226, 152, 37, 179, 16, 145, 34, 136,
	                                              54, 208, 148, 206, 143, 150, 219, 189, 241, 210, 19, 92, 131, 56, 70, 64,
	                                              30, 66, 182, 163, 195, 72, 126, 110, 107, 58, 40, 84, 250, 133, 186, 61,
	                                             202, 94, 155, 159, 10, 21, 121, 43, 78, 212, 229, 172, 115, 243, 167, 87,
	                                               7, 112, 192, 247, 140, 128, 99, 13, 103, 74, 222, 237, 49, 197, 254, 24,
	                                             227, 165, 153, 119, 38, 184, 180, 124, 17, 68, 146, 217, 35, 32, 137, 46,
	                                              55, 63, 209, 91, 149, 188, 207, 205, 144, 135, 151, 178, 220, 252, 190, 97,
	                                             242, 86, 211, 171, 20, 42, 93, 158, 132, 60, 57, 83, 71, 109, 65, 162,
	                                              31, 45, 67, 216, 183, 123, 164, 118, 196, 23, 73, 236, 127, 12, 111, 246,
	                                             108, 161, 59, 82, 41, 157, 85, 170, 251, 96, 134, 177, 187, 204, 62, 90,
	                                             203, 89, 95, 176, 156, 169, 160, 81, 11, 245, 22, 235, 122, 117, 44, 215,
	                                              79, 174, 213, 233, 230, 231, 173, 232, 116, 214, 244, 234, 168, 80, 88, 175   };
	
	public int[] Table_Integer_To_Exponent_Of_Alpha()
	{
            return Integer_To_Exponent_Of_Alpha;
	}
        
        public ArrayList<Integer> IntegerTextFileToArrayList(String filepath)
        {
            ArrayList<Integer> myarray = new ArrayList<Integer>();
            Scanner inputReader = null;
            try {
                    File inputFile = new File(filepath);
                    inputReader = new Scanner(inputFile);
                    while (true)
                    {
                        String inputData = inputReader.nextLine();
                        Integer myint = Integer.parseInt(inputData);
                        myarray.add(myint);
                    }
            }
            catch (NoSuchElementException ex) {
                //EOF reached
                inputReader.close();
                return myarray;
            } 
            catch (Exception ex) {
                System.out.println("IntegerTextFileToArrayList: An exception occurred!");
                ex.printStackTrace();
                return null;
            }
        }//IntegerTextFileToArrayList
        
        public void ArrayListToIntegerTextFile(String outputfilepath, ArrayList<Integer> myarray)
        {
            try {
                FileWriter myWriter = new FileWriter(outputfilepath);
                for (int ii=0; ii < myarray.size(); ii++) {
                    myWriter.write(myarray.get(ii).toString());
                    myWriter.write('\n');
                }	
                myWriter.close();
            }
            catch (Exception ex) {
		System.out.println("ArrayListToIntegerTextFile: An exception occurred!");
		ex.printStackTrace();
            }
        }//ArrayListToIntegerTextFile
        
        public void ArrayListTo_Integer_ASCII_TextFile(String outputfilepath, ArrayList<Integer> myarray)
        {
            try {
                FileWriter myWriter = new FileWriter(outputfilepath);
                for (int ii=0; ii < myarray.size(); ii++) {
                    myWriter.write(myarray.get(ii).toString());
                    myWriter.write(" = ");
                    myWriter.write((char)((int)myarray.get(ii)));
                    myWriter.write('\n');
                }	
                myWriter.close();
            }
            catch (Exception ex) {
		System.out.println("ArrayListTo_Integer_ASCII_TextFile: An exception occurred!");
		ex.printStackTrace();
            }
        }//ArrayListTo_Integer_ASCII_TextFile
        
        
        public ArrayList<Integer> Add_ByteEncIndAndLength(ArrayList<Integer> block1,
                                                          ArrayList<Integer> block2,
                                                          ArrayList<Integer> block3,
                                                          ArrayList<Integer> block4,
                                                          int max_length)
        {
            //encoding == ENCODING_INDICATOR.BYTE == 0100 (4-bits)
            //length (8-bits)
            //d.c (each 8-bits, so not on byte-boundary)
            //add at end
            //encoding == ENCODING_INDICATOR.TERMINATOR == 0000 (4-bits)
            //pad to 8-bits-byte
            //pad 236 and 17 as bytes till max length
            ArrayList<Integer> good_data_only_in_correct_sequence = new ArrayList<Integer>();
            good_data_only_in_correct_sequence.addAll(block1);
            good_data_only_in_correct_sequence.addAll(block2);
            good_data_only_in_correct_sequence.addAll(block3);
            good_data_only_in_correct_sequence.addAll(block4);
            
            ArrayList<Integer> results = new ArrayList<Integer>();
            String byte_enc_str = "0100";
            Integer length = good_data_only_in_correct_sequence.size();
            String len_str = Integer.toBinaryString(length);
            String padded8str = "00000000";
            if (len_str.length() < 8)
                len_str = padded8str.substring(0, 8 - len_str.length()) + len_str;
            String ind_str = byte_enc_str + len_str;
            results.add(ConvertBinaryStringToPositiveInteger(ind_str.substring(0,8)));
            String startOfByte = ind_str.substring(8);
            for (int ii=0; ii < good_data_only_in_correct_sequence.size(); ii++) {
                padded8str = "00000000";
                String tempstr = Integer.toBinaryString(good_data_only_in_correct_sequence.get(ii));
                if (tempstr.length() < 8)
                    tempstr = padded8str.substring(0, 8 - tempstr.length()) + tempstr;
                
                tempstr = startOfByte + tempstr;
                results.add(ConvertBinaryStringToPositiveInteger(tempstr.substring(0,8)));
                startOfByte = tempstr.substring(8);
            }
            
            //startOfByte.length() == 4
            startOfByte += "0000"; //encoding == ENCODING_INDICATOR.TERMINATOR == 0000 (4-bits)
            results.add(ConvertBinaryStringToPositiveInteger(startOfByte));
            
            int end_padding_len = max_length - results.size();
            for (int ii=0; ii < end_padding_len; ii++)
            {
                if ((ii % 2) == 0)
                    results.add(236);
                else
                    results.add(17);
            }
       
                
            return results;
        }//Add_ByteEncIndAndLength
        
        
        public void Print_GF7_GenPolynomial(Integer alpha, Integer StartWithAlphaExpCoefficient)
        {
            System.out.println();
            System.out.println("***** Print_GF7_GenPolynomial_Start_With_AlphaToThe0  alpha==" + (int)alpha + "*****");
            System.out.println("***** Remember this is GF7-world not normal base10-decimal-world  alpha==" + (int)alpha + "*****");
            
            //GF-world not normal base10-decimal-world
            //(x + alpha^0)(x + alpha^1)
            //x^2 + (alpha^1 + alpha^0)x + alpha^1
            //(alpha^0)(x^2) + (alpha^1 + alpha^0)x^1 + alpha^1(x^0)
            //(alpha^0)(x^2) + (alpha^1 + 1)x^1 + alpha^1
            
            ArrayList<ArrayList<Integer>> big_gen_poly_alphadecimalcoeff_xtothexexp = new ArrayList<ArrayList<Integer>>();
            //minimum amount of error-correction-codewords == 2
            //maximum amount of error-correction-codewords == 7
            big_gen_poly_alphadecimalcoeff_xtothexexp.add(null); //0 error_correction_codewords
            big_gen_poly_alphadecimalcoeff_xtothexexp.add(null); //1 error_correction_codewords

            //term1 will be the big term that keeps getting bigger
            ArrayList<Integer> bigterm_alphadecimal = new ArrayList<Integer>();
            ArrayList<Integer> bigterm_x_exp = new ArrayList<Integer>();
            //in GF(whatever such as 256 or 7) world 
            //- is same as + in multiplication
            //n == number of E.C.C
            //(x - alpha^0).....(x - alpha^(n-1))
            //
            //so for 2 E.C.C
            //(x - alpha^0)(x - alpha^1)
            //((alpha^0)x^1 + (alpha^0)x^0) * ((alpha^0)x^1 + (alpha^1)x^0)
            //
            //term2 will be (x + alpha^(n-1)) 
            //it will always be total size 4
            Integer[] smallterm_alphadecimal = {0, 0};
            Integer[] smallterm_x_exp = {0, 0};

            ArrayList<Integer> alphadecimalcoeff_xtothexexp = new ArrayList<Integer>();
            //index0 is for biggest x to the xexp

            //answer will have 1 more than highest-term-count
            //since we start with (x + 1)(x + alpha^1) for n==2
            //that means there will be 3 terms
            //with highest x exponent of 2 as in x^2
            alphadecimalcoeff_xtothexexp.add(0); //coefficient for x^2
            alphadecimalcoeff_xtothexexp.add(0); //coefficient for x^1
            alphadecimalcoeff_xtothexexp.add(0); //coefficient for x^0
            Integer highest_x_exp = 2;

            for (Integer nn = 2; nn <= 7; nn++) {
                if (nn == 2) {
                    if (StartWithAlphaExpCoefficient == 0)
                    {
                         bigterm_alphadecimal.add(1);
                         bigterm_x_exp.add(1);
                         bigterm_alphadecimal.add(1);
                         bigterm_x_exp.add(0);
                    
                         smallterm_alphadecimal[0] = 1;
                         smallterm_x_exp[0] = 1;
                         smallterm_alphadecimal[1] = alpha;
                         smallterm_x_exp[1] = 0;
                    }
                    else if (StartWithAlphaExpCoefficient == 1)
                    {
                         bigterm_alphadecimal.add(1);
                         bigterm_x_exp.add(1);
                         bigterm_alphadecimal.add(alpha);
                         bigterm_x_exp.add(0);
                    
                         smallterm_alphadecimal[0] = 1;
                         smallterm_x_exp[0] = 1;
                         smallterm_alphadecimal[1] = (int)Math.pow(alpha, 2);
                         smallterm_x_exp[1] = 0;
                    }
                } else {
                    bigterm_alphadecimal.clear();
                    bigterm_x_exp.clear();
                    
                    //java has changed
                    //can not use
                    //coeff_x_to_the_xexp.clear
                    //must use
                    //new ArrayList<Integer>()
                    //because of gen_poly below
                    alphadecimalcoeff_xtothexexp = new ArrayList<Integer>();

                    //(x - a^0).....(x - a^(n-2))
                    
                    ArrayList<Integer> gen_poly = big_gen_poly_alphadecimalcoeff_xtothexexp.get(nn - 1);
                    for (int ii = 0; ii < gen_poly.size(); ii++) {
                        bigterm_alphadecimal.add(gen_poly.get(ii));
                        bigterm_x_exp.add(gen_poly.size() - ii - 1);
                        alphadecimalcoeff_xtothexexp.add(0);
                    }
                    alphadecimalcoeff_xtothexexp.add(0);

                    //smallterm will be 
                    //(a^0x^1 - a^(n-1)x^0)
                    //or
                    //(a^0x^1 - a^nx^0)
                    //it will always be total size 4
                    if (StartWithAlphaExpCoefficient == 0)
                    {
                        smallterm_alphadecimal[0] = 1;
                        smallterm_x_exp[0] = 1;
                        smallterm_alphadecimal[1] = (int)Math.pow(alpha, nn - 1);
                        smallterm_x_exp[1] = 0;
                    }
                    else if (StartWithAlphaExpCoefficient == 1)
                    {
                        smallterm_alphadecimal[0] = 1;
                        smallterm_x_exp[0] = 1;
                        smallterm_alphadecimal[1] = (int)Math.pow(alpha, nn);
                        smallterm_x_exp[1] = 0;
                    }

                    highest_x_exp = nn;
                }

                //  (x^(n-1) + ..... blah blah)* (x + alpha^(n-1))
                for (int bigterm_ii = 0; bigterm_ii < bigterm_alphadecimal.size(); bigterm_ii++) {
                    Integer temp1_alphadecimal = bigterm_alphadecimal.get(bigterm_ii);
                    int temp1_x_exp = bigterm_x_exp.get(bigterm_ii);
                    for (int smallterm_ii = 0; smallterm_ii < 2; smallterm_ii++) {
                        Integer temp2_alphadecimal = smallterm_alphadecimal[smallterm_ii];
                        Integer temp2_x_exp = smallterm_x_exp[smallterm_ii];

                        int x_exp = temp1_x_exp + temp2_x_exp;
                        Integer alphadecimal = temp1_alphadecimal * temp2_alphadecimal;
                        
                        //alpha is normally 2
                        //but in GF(7)-example
                        //alpha is 3

                        if (x_exp == highest_x_exp) {
                            //alphadecimal should always be 1 but keep code for human-readability-understanding
                            alphadecimalcoeff_xtothexexp.set(0, alphadecimal);
                        } else {
                            int index = highest_x_exp - x_exp;
                            Integer temp = alphadecimalcoeff_xtothexexp.get(index);
                            temp += alphadecimal;
                            temp %= 7;
                            alphadecimalcoeff_xtothexexp.set(index, temp);
                        }
                    }//term2_ii
                }//term1_ii

                big_gen_poly_alphadecimalcoeff_xtothexexp.add(alphadecimalcoeff_xtothexexp);
                System.out.println();
                System.out.print("E.C.C count == " + nn + " : ");
                
                for (int ii=0; ii < alphadecimalcoeff_xtothexexp.size(); ii++)
                {
                    System.out.print(alphadecimalcoeff_xtothexexp.get(ii) + "x^" + (alphadecimalcoeff_xtothexexp.size() - ii - 1));
                    if (ii != (alphadecimalcoeff_xtothexexp.size() - 1))
                        System.out.print(" + ");
                }
                
            }//for nn
        }//Print_GF7_GenPolynomial
        
        
    public int GF256multiply(int A, int B)
    {
        //GF(256) world
        if ((A==0) || (B==0))
            return 0;
        else
        {
            int A_as_alpha_exp = Table_Integer_To_Exponent_Of_Alpha()[A];
            int B_as_alpha_exp = Table_Integer_To_Exponent_Of_Alpha()[B];
            int total_exp = A_as_alpha_exp + B_as_alpha_exp;
            if (total_exp >= 256)
                total_exp %= 255;
            int gf256num = Table_Exponent_Of_Alpha_To_Integer()[total_exp];
            return gf256num;
        }
    }//GF256multiply
    
    
    public int GF256divide(int top, int bottom) throws Exception
    {
        //GF(256) world
        if (bottom==0)
            throw new Exception("GF256divide: bottom==0");
        
        if (top==0)
            return 0;
        
        int top_alpha_exp = Table_Integer_To_Exponent_Of_Alpha()[top];
        int bottom_alpha_exp = Table_Integer_To_Exponent_Of_Alpha()[bottom];
        int total_exp = top_alpha_exp - bottom_alpha_exp;
        if (total_exp < 0)
        {
            /***********************************************
            Negative exponents: As with regular arithmetic, 
            a^{-n}=(a^{-1})^{n}). 
            Therefore, calculating (alpha ^{-n}) involves 
            finding the multiplicative inverse of alpha(i.e., (alpha ^{-1}) 
            and then raising that result to the power of n.
            A multiplicative inverse is a number that, when multiplied by a given number, equals 1
            ************************************************/
            //in GF(256)
            //2 == 2^1
            //2^1 * 2^254 == 2^(1+254) == 2^255 == 1
            int pos_exp = total_exp * -1;
            total_exp = 254 + pos_exp;
        }
        if (total_exp >= 256)
            total_exp %= 255;

        int gf256num = Table_Exponent_Of_Alpha_To_Integer()[total_exp];
        return gf256num;
    }//GF256divide
                
        
}//class
