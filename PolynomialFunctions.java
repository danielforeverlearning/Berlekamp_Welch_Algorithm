
package Berlekamp_Welch_Algorithm.test;

import java.util.ArrayList;

public class PolynomialFunctions {
    
    public int   POLYDIVLIN_Remainder;
    public int[] POLYDIVLIN_qvec;
    
    public int[] poly_long_division_answer;
    public int[] poly_long_division_remainder;
    
    public PolynomialFunctions()
    {
        
    }//constructor
    
    
    public void Debug_Print(int[] poly)
    {
        System.out.println();
        int len = poly.length;
        for (int exp = len - 1; exp >= 0; exp--)
        {
            System.out.print(poly[exp]);
            if (exp==0)
                System.out.print(" ");
            else if (exp==1)
                System.out.print("x + ");
            else
                System.out.print("x^" + exp + " + ");
        }
        System.out.println();
    }//Debug_Print
    
    
    //page119
    //Algorithm A.3.2 Polynomial multiplication by x − α
    //(1) Compute the coefficients of f(x) = (x − α)p(x)
    //(2) Determine the degree of p(x) based on the length of ~p
    public int[] POLYMULTLIN(int[] p_vec, int a_ii)
    {
        int n = p_vec.length - 1;
        int[] f = new int[n + 2];
        
        //GF(256) world
        Tools tool = new Tools();
        //f[0] = a_ii * p_vec[0]
        int gf256num = tool.GF256multiply(a_ii, p_vec[0]);
        f[0] = gf256num;
        
        for (int ii=1; ii <= n; ii++)
        {
            int p_ii = p_vec[ii];
            gf256num = tool.GF256multiply(a_ii, p_ii);
            int p_iiminus1 = p_vec[ii-1];
            gf256num = p_iiminus1 ^ gf256num;
            f[ii] = gf256num;
        }
        f[n + 1] = p_vec[n];
        return f;
        
    }//POLYMULTLIN
    
    
    //Compute the coefficients of
    //p(x) = PRODUCT( αi ∈ ~α (x - αi))
    //page121
    //Algorithm A.3.4 Constructing a polynomial from its roots
    public int[] POLYROOTS(int[] a_vec)
    {
        int n = a_vec.length;
        int[] p_vec = new int[1];
        p_vec[0] = 1;
        for (int ii=0; ii < n; ii++)
            p_vec = POLYMULTLIN(p_vec, a_vec[ii]);
        return p_vec;
    }//POLYROOTS
    
    
    //Algorithm A.3.1 Polynomial evaluation (Horner’s Algorithm)
    //(1) Evaluate the polynomial ~p at the points ~a
    //(2) Determine the degree of p(x) based on the length of ~p
    //(3) m is the number of points to evaluate
    //(4) Loop across all of the points
    public int[] POLYEVAL(int[] p_vec, int[] a_vec)
    {
        Tools tool = new Tools();
        
        int n = p_vec.length - 1;
        int m = a_vec.length;
        int[] y = new int[m];
        for (int j=0; j < m; j++)
        {
            y[j] = p_vec[n];
            for (int i=(n-1); i >= 0; i--)
            {
                int gf256num = tool.GF256multiply(a_vec[j], y[j]);
                y[j] = p_vec[i] ^ gf256num;
            }
        }
        return y;
    }//POLYEVAL
    
    
    /*********************************************************************************************
        page25
        Example 2.2.1. Continuing with Example 2.1.1, we may chose I = {0, 1, 2, 3, 4}. 
        Then we must interpolate the points (0, 233), (1, 47), (2, 87), (3, 131), and (4, 168). 
        This means that  m(x) == 18x^4 + 7x^3 + 0x^2 + 211x + 233
        
        Example 2.3.1. Consider the first message from Example 1.4.1, ~m = (233, 211, 0, 7, 18)T
        with l = 8 and k = 5. We have the message polynomial
        
            µ(x) = 65x^4 + 112x^3 + 10x^2 + x + 233
            
        and systematic codeword ~c = (233, 211, 0, 7, 18, 166, 14, 135)T
    ************************************************************************************************/
    
    public int[] my_LAGRANGEINTERPOLATE(int[] a_vec, int[] B_vec) throws Exception
    {
        if (a_vec.length != B_vec.length)
            throw new Exception("my_LAGRANGEINTERPOLATE: a_vec.length != B_vec.length");
        
        Tools tool = new Tools();  
        int[] all_roots_poly = this.POLYROOTS(a_vec);
        
        //Initialize ~p as a length n vector of zeros
        int[] p = new int[a_vec.length];
        for (int ii=0; ii < a_vec.length; ii++)
            p[ii] = 0;
        
        for (int ii=0; ii < a_vec.length; ii++)
        {
            int[] bottom_poly = { a_vec[ii], 1 };
            //System.out.println("ii=" + ii);
            this.poly_long_division(all_roots_poly, bottom_poly);
            //this.Debug_Print(this.poly_long_division_answer);
            //this.Debug_Print(this.poly_long_division_remainder);
            //System.out.println();
            
            if (this.poly_long_division_remainder.length==1 && this.poly_long_division_remainder[0]==0)
            {
                int B_ii = B_vec[ii];
                
                //56x^4 + 2x^3 + 238x^2 + 238x + 233 
                //check_c = 233,211,0,54,18
                //
                //98x^7 + 51x^6 + 153x^5 + 169x^4 + 171x^3 + 49x^2 + 193x + 233 
                //check_c = 233,211,0,7,18,166,14,253
                //int[] temparray = { a_vec[ii] };
                //int[] bottomarray = this.POLYEVAL(this.poly_long_division_answer, temparray);
                //if (bottomarray.length != 1)
                //    throw new Exception("my_LAGRANGEINTERPOLATE: bottomarray.length != 1");
                //int mult = tool.GF256divide(B_ii, bottomarray[0]);
                
                
                //56x^4 + 2x^3 + 238x^2 + 238x + 233 
                //check_c = 233,211,0,54,18
                //
                //98x^7 + 51x^6 + 153x^5 + 169x^4 + 171x^3 + 49x^2 + 193x + 233 
                //check_c = 233,211,0,7,18,166,14,253
                int bottom = calc_product_aii_minus_ajj_where_ii_cannotequal_jj(a_vec, ii);
                int mult = tool.GF256divide(B_ii, bottom);
                
                 
                //56x^4 + 2x^3 + 238x^2 + 238x + 233 
                //check_c = 233,211,0,54,18
                //int temp2 = this.polycalc(this.poly_long_division_answer, a_vec[ii]);
                //int mult = tool.GF256divide(B_ii, temp2);
                
                
                //166x^4 + 222x^3 + 98x^2 + 231x + 74 
                //check_c = 74,183,0,3,233
                //int temp2 = this.polycalc(this.poly_long_division_answer, a_vec[ii]);
                //int top = tool.GF256multiply(B_ii, temp2);
                //int bottom = calc_product_aii_minus_ajj_where_ii_cannotequal_jj(a_vec, ii);
                //int mult = tool.GF256divide(top, bottom);
                
                //58x^4 + 50x^3 + 8x^2 + 214x + 155 
                //check_c = 155,77,0,122,224
                //int bottom = calc_product_aii_minus_ajj_where_ii_cannotequal_jj(a_vec, ii);
                //bottom = tool.GF256multiply(bottom, this.polycalc(this.poly_long_division_answer, a_vec[ii]));
                //int mult = tool.GF256divide(B_ii, bottom);
                
                int[] a_di_poly = poly_multiply_by_gf256num(this.poly_long_division_answer, mult);
                p = poly_add(p, a_di_poly);
            }
            else
                throw new Exception("my_LAGRANGEINTERPOLATE: poly_long_division_remainder != 0");
        }
        
        return p;
    }//my_LAGRANGEINTERPOLATE
    
    public int[] poly_multiply_by_gf256num(int[] poly, int gf256num)
    {
        Tools tool = new Tools();
        int[] answer_poly = new int[poly.length];
        for (int ii=0; ii < poly.length; ii++)
            answer_poly[ii] = tool.GF256multiply(poly[ii], gf256num);
        
        return answer_poly;
    }
    
    public int[] poly_add(int[] polyA, int[] polyB) throws Exception
    {
        if (polyA.length != polyB.length)
            throw new Exception("poly_add: polyA.length != polyB.length");
        
        int[] answer_poly = new int[polyA.length];
        for (int ii=0; ii < answer_poly.length; ii++)
        {
            answer_poly[ii] = polyA[ii] ^ polyB[ii];
        }
        return answer_poly;
    }
    
    public int[] roots_ii_cannotequal_jj(int[] a_vec, int ii)
    {
        //(x+1)(x+2) ==
        //x^2 + 2x + 1x + 2 ==
        //x^2 + (2 ^ 1)x + 2 ==
        //x^2 + 3x + 2
        
        //(x^2 + 3x + 2)(x+3) ==
        //x^3 + 3x^2 + 2x + 3x^2 + (2^25)(2^25)x + (2^1)(2^25) ==
        //x^3 + (3 ^ 3)x^2 + (2 ^ 2^50)x + 2^26 ==
        //x^3 + 0x^2 + (2 ^ 5)x + 6 ==
        //x^3 + 0x^2 + 7x + 6
        
        //(x^3 + 0x^2 + 7x + 6)(x+4) ==
        //x^4 + 0x^3 + 7x^2 + 6x + 4x^3 + GF(256)(7*4)x + GF(256)(6*4) ==
        //x^4 + 4x^3 + 7x^2 + (6 ^ 2^200)x + (2^28) ==
        //x^4 + 4x^3 + 7x^2 + (6 ^ 28)x + 24 ==
        //x^4 + 4x^3 + 7x^2 + 26x + 24
        
        //(0^1)(0^2)(0^3)(0^4) ==
        //GF256(1 * 2 * 3 * 4) ==
        //2^0 * 2^1 * 2^25 * 2^2 ==
        //2^(0+1+25+2)==
        //2^28 ==
        //24
        Tools tool = new Tools();
        int[] deg = null;
        
        for (int jj=0; jj < a_vec.length; jj++)
        {
            if (ii != jj)
            {
                if (deg == null) //empty
                {
                    deg = new int[2];
                    deg[0] = a_vec[jj];
                    deg[1] = 1; //always1 because 1x^1
                }
                else //it is long not empty
                {
                    deg = this.POLYMULTLIN(deg, jj);
                }
            }
        }
        return deg;
    }//roots_ii_cannotequal_jj
    
    
    private int calc_product_aii_minus_ajj_where_ii_cannotequal_jj(int[] a_vec, int ii)
    {
        Tools tool = new Tools();
        int mult = 1;
        int a_ii = a_vec[ii];
        for (int jj=0; jj < a_vec.length; jj++)
        {
            if (ii != jj)
            {
                int a_jj = a_vec[jj];
                int sum = a_ii ^ a_jj;
                mult = tool.GF256multiply(sum, mult);
            }
        }
        return mult;
    }//calc_product_aii_minus_ajj_where_ii_cannotequal_jj
    
    
    //Construct di(x) = PRODUCT(j!=i, from j=0 to n-1 (x − αj)) = d(x)/(x−αi)
    
    
    public int polycalc(int[] poly, int x)
    {
        //poly[0] == coefficient for x^0
        //poly[1] == coefficient for x^1
        //poly[2] == coefficient for x^2
        //.....
        //poly[poly.length - 1] == coefficient for x^(poly.length - 1)
        
        Tools tool = new Tools();
        
        int sum = poly[0];
        for (int x_exp=1; x_exp < poly.length; x_exp++)
        {
            if (x!=0) //0 raised to anything positive still 0, sum remains same
            {
                int x_as_alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[x];
                int total_exp = x_as_alpha_exp * x_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256_x_raised_to_x_exp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                int coeff = poly[x_exp];
                int gf256num = tool.GF256multiply(coeff, gf256_x_raised_to_x_exp);
                sum ^= gf256num;
            }
        }
        return sum;
    }//polycalc
    
    
    public void poly_long_division(int[] top_poly, int[] bottom_poly)
    {
        //top_poly[0] == coeff for x^0
        //bottom_poly[0] == coeff for x^0
        
        //top_poly[1] == coeff for x^1
        //bottom_poly[1] == coeff for x^1
        
        //top_poly[top_poly.length - 1] == coeff for x^(top_poly.length - 1)
        //bottom_poly[bottom_poly.length - 1] == coeff for x^(bottom_poly.length - 1)
        
        Tools tool = new Tools();
        
        ArrayList<Integer> numerator = new ArrayList<Integer>();
        for (int ii=top_poly.length-1; ii >= 0; ii--)
            numerator.add(top_poly[ii]);
        
        ArrayList<Integer> denominator = new ArrayList<Integer>();
        for (int ii=bottom_poly.length-1; ii >= 0; ii--)
            denominator.add(bottom_poly[ii]);
        
        ArrayList<Integer> long_division_answer = new ArrayList<Integer>();
        while (numerator.size() >= denominator.size())
        {
            int mult = 0;
            int mult_alphaexp = 0;
            for (int ii=0; ii < denominator.size(); ii++)
            {
                if (ii==0)
                {
                    mult = numerator.get(ii);
                    long_division_answer.add(mult);
                    numerator.set(ii, 0);
                }
                else
                {
                    if (mult != 0)
                    {
                        mult_alphaexp = tool.Table_Integer_To_Exponent_Of_Alpha()[mult];
                        int temp_top = numerator.get(ii);                    
                        int temp_bottom = denominator.get(ii);
                        if (temp_bottom != 0)
                        {
                            int temp_bottom_alphaexp = tool.Table_Integer_To_Exponent_Of_Alpha()[temp_bottom];
                            int total_alphaexp = temp_bottom_alphaexp + mult_alphaexp;
                            if (total_alphaexp >= 256)
                                total_alphaexp %= 255;

                            int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_alphaexp];
                            int xornum  = temp_top ^ temp;
                            numerator.set(ii, xornum);
                        }
                    }
                    //else numerator[ii] remains same
                }       
            }//for
            
            numerator.remove(0);
        }//while
        
        poly_long_division_answer = new int[long_division_answer.size()];
        //System.out.println();
        //System.out.print("long_division answer = ");
        for (int ii=0; ii < long_division_answer.size(); ii++)
        {
            //System.out.print(long_division_answer.get(ii) + "x^" + (long_division_answer.size() - ii - 1));
            //if (ii != (long_division_answer.size() - 1))
            //     System.out.print(" + ");
            
            poly_long_division_answer[long_division_answer.size() - ii - 1] = long_division_answer.get(ii); //because coefficients are reversed between the arrays
        }
        //System.out.println();
        
        poly_long_division_remainder = new int[numerator.size()];
        //System.out.print("long_division remainder = ");
        for (int ii=0; ii < numerator.size(); ii++)
        {
            //System.out.print(numerator.get(ii) + "x^" + (numerator.size() - ii - 1));
            //if (ii != (numerator.size() - 1))
            //    System.out.print(" + ");
            
            poly_long_division_remainder[numerator.size() - ii - 1] = numerator.get(ii); //because coefficients are reversed between the arrays 
        }
    }//poly_long_division
    
    
    public void test_POLYROOTS()
    {
        //Compute the coefficients of
        //p(x) = PRODUCT( αi ∈ ~α (x - αi))
        //page121
        //Algorithm A.3.4 Constructing a polynomial from its roots
        
        int[] a_vec = { 0, 1, 2, 3, 4, 5};
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
        
        //GF(256) (x^5 + 4x^4 + 7x^3 + 26x^2 + 24x)(x + 5) == x^6 + (5 ^ 4)x^5 + (4*5 + 7)x^4      + (7*5 + 26)x^3        + (26*5 + 24)x^2       + (24*5)x
        //                                                 == x^6 +  1x^5      + (2^(2+50) ^ 7)x^4 + (2^(198+50) ^ 26)x^3 + (2^(105+50) ^ 24)x^2 + (2^(28+50)x
        //                                                 == x^6 +  1x^5      + (20 ^ 7)x^4       + (2^248  ^  26)x^3    + (2^155  ^  24)x^2    + (2^78)x
        //                                                 == x^6 +  1x^5      +  19x^4            + (27     ^  26)x^3    + (114    ^  24)x^2    +  120x
        //                                                 == x^6 +  1x^5      +  19x^4            +    1x^3              +  106x^2              +  120x
        int[] pdeg5 = POLYROOTS(a_vec);
        Debug_Print(pdeg5); 
    }
}//class


