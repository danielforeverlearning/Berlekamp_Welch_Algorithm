
package Berlekamp_Welch_Algorithm.test;

import java.util.ArrayList;

public class PolynomialFunctions {
    
    public PolynomialFunctions()
    {
        
    }//constructor
    
    
    
    
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
    
    
}//class
