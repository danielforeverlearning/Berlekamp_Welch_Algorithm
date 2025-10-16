
package Berlekamp_Welch_Algorithm.test;

import java.util.ArrayList;
import java.math.BigInteger;

public class Main {

    
    public static void main(String[] args) {
        try
        {
            //ok try to code 
            //5.2.6 The Algorithm
            //pages 62-66
            
            /******************************************************
            page 46 Example 4.2.1
            ~m = (233, 211, 0, 7, 18)T
            Encoding systematically yields
            ~c = (233, 211, 0, 7, 18, 166, 14, 135)T (correct-bytes)
            ~c is then transmitted and an error occurs resulting in a received word
            ~r = (233, 117, 0, 7, 18, 166, 14, 135)T
            *******************************************************/
            //page64 Example 5.2.7
            //l == 8 (8-bits per byte)
            //received-bytes
            //int[] r_vector = { 233, 117, 0, 7, 18, 166, 14, 135 };
            //int k = 5; //number of message-bytes
            //Pages62_66_Modified_Welch_Berlekamp_Algorithm algo = new Pages62_66_Modified_Welch_Berlekamp_Algorithm(r_vector,k);
            //boolean success = algo.go();
            
            //Ok i think i figured it out i think .....
            //Berlekamp-Welch and Modified-Berlekamp-Welch
            //algorithms work on systematic-encoding where
            //you use a different polynomial u(x) called the  "message-polynomial" and you use ai=0,1,2,3,4,5 .....i-1
            //as inputs to x in the message-polynomial to get 
            //the first k bytes which are the message-bytes followed by the error-correction-bytes which
            //you put after the message-bytes.
            //BUT QR-codes first use polynomial-long-division on a polynomial by
            //a generator-polynomial and the error-correction-bytes which you
            //put after the message-bytes are the remainder-bytes from the polynomial-long-division.
            //This message-bytes and error-correction-bytes are not u(x)
            //you have to calculate u(x) first.
            
            /**********************************************************************
            2.3 Systematic Encoding
            It would be nice if the message appeared as a subset of the codeword, i.e.
            ~c = (m0, . . . mk−1, ck, . . . , cn−1)T   .
            In this case, recovering a message from the codeword is trivial; 
            one simply reads the message
            from the first k characters of the codeword. 
            This type of encoding is termed systematic,
            and may be achieved by altering the encoding algorithm.
            Previously we had encoded using the rule ~c = m(~α). 
            We will now introduce a new
            polynomial µ(x) that will generate the encoding systematically via ~c = µ(~α).
            µ(x) is the
            unique degree k − 1 or less polynomial such that µ(αi) = mi for 0 ≤ i < k. 
            We can find µ(x) explicitly using Lagrange interpolation :
             
            µ(x) = SERIESSUM( i=0 to k-1 (mi) * SERIESPRODUCT( 0 <= j < k, j!=i (x-aj)/(ai-aj)) )
            
            Given this definition, it is easy to see that ~c is indeed systematic.
            Systematic message encoding is implemented in Algorithm 2.3.1 
            and costs 2nk + 4k^2 − 2n + k operations. 
            Recovering the message from the codeword is trivial when using a
            systematic encoding as it is simply read from the first k characters.
            
            
            ok ..... holy cow,
            so here is example:
            page 24 arizona-paper:
            
            Example 2.1.1. Consider the first message from Example 1.4.1, 
            ~m = (233, 211, 0, 7, 18)T
            with l = 8 and k = 5. 
            We have the message polynomial
            m(x) = 233 + 211x + 0x^2 + 7x^3 + 18x^4
            To create a length n = 8 codeword we chose ~α = (0, 1, 2, 3, 4, 5, 6, 7)T
            Thus the codeword
            ~c is given by
                c0 = m(0) = 233
                c1 = m(1) = 47
                c2 = m(2) = 87
                c3 = m(3) = 131
                c4 = m(4) = 168
                c5 = m(5) = 2
                c6 = m(6) = 134
                c7 = m(7) = 62
                
            i.e. ~c = (233, 47, 87, 131, 168, 2, 134, 62)T
            
            
            Example 2.2.1. Continuing with Example 2.1.1, we may chose I = {0, 1, 2, 3, 4}. Then
            we must interpolate the points (0, 233), (1, 47), (2, 87), (3, 131), and (4, 168). 
            This means that
            
            m(x) = 233 * ((x − 1)/(0 − 1))((x − 2)/(0 − 2))((x − 3)/(0 − 3))((x − 4)/(0 − 4))
                   + 
                   47 * ((x − 0)/(1 − 0))((x − 2)/(1 − 2))((x − 3)/(1 − 3))((x − 4)/(1 − 4))
                   +
                   87 * ((x − 0)/(2 − 0))((x − 1)/(2 - 1))((x - 3)/(2 - 3))((x - 4)/(2 - 4))
                   + 
                   131 * ((x - 0)(3 - 0))((x - 1)/(3 - 1))((x - 2)(3 - 2))((x - 4)/(3 - 4))
                   +
                   168 * ((x - 0)(4 - 0))((x - 1)/(4 - 1))((x - 2)(4 - 2))((x - 3)/(4 - 3))
                 
                 = 18x^4 + 7x^3 + 0x^2 + 211x + 233
            
            
            ************************************************************************/
            
            Tools tool = new Tools();
            int[] u = { 65, 112, 10, 1, 233 }; //65x^4 + 112x^3 + 10x^2 + 1x^1 + 233x^0 
            int k = 5;
            ArrayList<Integer> c = new ArrayList<Integer>();
            
            for (int ai=0; ai <= 7; ai++)  //n == 8
            {
                if (ai==0)
                    c.add(u[4]);
                else
                {
                    int sum=0;
                    int x = ai;
                    int x_as_alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[x];
                    for (int uu=0; uu < u.length; uu++)
                    {
                        int coeff = u[uu];
                        if (coeff != 0)
                        {
                            int coeff_as_alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[coeff];
                            int x_exp = u.length - uu - 1;
                            int total_exp = coeff_as_alpha_exp + (x_as_alpha_exp * x_exp);
                            if (total_exp >= 256)
                                total_exp %= 255;
                            int gf256num = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                            sum ^= gf256num;
                        }
                    }
                    c.add(sum);
                }
            }
            //n  == c.length == 8
            //t^ == floor((n-k)/2)) == 1 == maximum number of error-bytes that can be detected and corrected
            System.out.println();
            System.out.print("m(x) = { ");
            for (int ii=0; ii < 5; ii++)
            {
                System.out.print(c.get(ii));
                if (ii != 4)
                    System.out.print(", ");
            }
            System.out.println("}");
            
            
            Integer[] r = new Integer[c.size()]; //received-bytes
            c.toArray(r);
            r[1] = 117; //error-index==1   correct-byte==211 changing it to 117 to simulate 1 error
            System.out.println();
            System.out.print("r(x) = { ");
            for (int ii=0; ii < r.length; ii++)
            {
                System.out.print(r[ii]);
                if (ii != (r.length-1))
                    System.out.print(", ");
            }
            System.out.println("}");
            Pages62_66_Modified_Welch_Berlekamp_Algorithm algo = new Pages62_66_Modified_Welch_Berlekamp_Algorithm(r, k);
            boolean success = algo.go();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
        
        //Compute the coefficients of
        //p(x) = PRODUCT( αi ∈ ~α (x - αi))
        //page121
        //Algorithm A.3.4 Constructing a polynomial from its roots
        
        PolynomialFunctions poly = new PolynomialFunctions();
        int[] a_vec = { 0, 1, 2, 3, 4};
        //GF(256) x + 0
        //GF(256) (x + 0)(x + 1) == x^2 + x
        //GF(256) (x^2 + x)(x + 2) == x^3 + 3x^2 + 2x
        //GF(256) (x^3 + 3x^2 + 2x)(x + 3) == x^4 + (2^25 ^ 2^25)x^3 + (2^50 ^ 2^1)x^2 + (2^26)x
        //                                 == x^4 + 0x^3 + (5 ^ 2)x^2 + 6x
        //                                 == x^4 + 0x^3 + 7x^2 + 6x
        //GF(256) (x^4 + 0x^3 + 7x^2 + 6x)(x + 4) == x^5 + (4 ^ 0)x^4 + (0 ^ 7)x^3 + (7*4 + 6)x^2 + (6*4)x
        //                                        == x^5 + 4x^4 + 7x^3 + (2^(198+2) ^ 6)x^2 + (2^(26+2))x
        //                                        == x^5 + 4x^4 + 7x^3 + (28 ^ 6)x^2 + 24x
        //                                        == x^5 + 4x^4 + 7x^3 + 26x^2 + 24x
        int[] pdeg5 = poly.POLYROOTS(a_vec);
        poly.Debug_Print(pdeg5); 
        
    }//main
    
}
