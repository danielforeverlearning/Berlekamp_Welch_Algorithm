import java.util.ArrayList;
import java.util.Scanner;


/********************************************
The Berlekamp–Welch algorithm, also known as the Welch–Berlekamp algorithm, 
is named for Elwyn R. Berlekamp and Lloyd R. Welch. 
This is a decoder algorithm that efficiently corrects errors in Reed–Solomon codes for an RS(n, k), 
code based on the Reed Solomon original view where a message 
m1, ..., mk is used as coefficients of a polynomial F(ai)

The goal of the decoder is to recover the original encoding polynomial 

 ********************************************/



public class TrueGF256_Berlekamp_Welch_algo_additiveinverseisitself {
    
    //n == 7 == message data bytes + E.C. bytes
    //k == 3 == message data bytes
    //so floor((n-k)/2)== maximum number of errors that can be corrected is 2
    
    //b is the received bytes, do not use b[0] so we will not be confused, b will have length==8
    //recv_max_index = b.length - 1;
    
    //e_count is the number of errors that can be corrected and detected
    //so if e_count==2 then we will have e0 and e1 and constrain e2==1 
    //so we will do
    //int[] e = new int[2]; // size is 2 for e[0] and e[1] and e[2]==1 constraint
    
    //a is the inputs, ai = i-1 : {0,1,2,3,4,5,6 ..... i-1}, do not use a[0] so we will not be confused
    //a will be same size of b
    //so
    //int[] a = new int[b.length];
    //it also uses recv_max_index
    
    //answer_matrix is rows==b.length x columns=1 represented by an array where answer_matrix[0] is not used so we will not be confused
    //int[] answer_matrix = new int[b.length];
    //it holds values -b[1]a[1]^2 .....-b[recv_max_index]a[recv_max_index]^2
    
    //int q_max_index = recv_max-index - e_count - 1
    //int[] q = new int[recv_max-index - e_count] we will use q[0]
    
    //so if e_count is 2
    //we calculate rows=recv_max_index x columns=1  matrix
    // e[0]
    // e[1]
    // q[0]
    // q[1]
    // q[2]
    // q[3]
    // q[4]
    //i guess we can call this matrix question_matrix
    
    //int[][] matrix is the matrix we try to get to the identity-matrix so that
    //
    //matrix x question_matrix == answer_matrix
    //
    //
    //matrix = b[1]  b[1]a[1]  -1  -a[1]  -a[1]^2  -a[1]^3  -a[1]^4
    //         b[2]  b[2]a[2]  -1  -a[2]  -a[2]^2  -a[2]^3  -a[2]^4
    //         b[3]  b[3]a[3]  -1  -a[3]  -a[3]^2  -a[3]^3  -a[3]^4
    //         b[4]  b[4]a[4]  -1  -a[4]  -a[4]^2  -a[4]^3  -a[4]^4
    //         b[5]  b[5]a[5]  -1  -a[5]  -a[5]^2  -a[5]^3  -a[5]^4
    //         b[6]  b[6]a[6]  -1  -a[6]  -a[6]^2  -a[6]^3  -a[6]^4
    //         b[7]  b[7]a[7]  -1  -a[7]  -a[7]^2  -a[7]^3  -a[7]^4
    //
    //int[][] matrix = new int[b.length][b.length]
    //and when we access matrix
    //we use
    //for (int row=1; row <= recv_max_index; row++)
    //for (int col=1; col <= recv_max_index; col++)
    //so do not use matrix[0][anything] or matrix[anything][0]
    
    int GF;

    int[] b = null;          //DO NOT USE b[0]
    int recv_max_index = 0;
    
    int e_count = 0;
    int[] e = null;
    
    int[] a = null;          //DO NOT USE a[0]
    
    int[] answer_matrix = null;
    
    int[][] matrix = null;
    
    int q_max_index = 0;
    int[] q = null;
    
    //F is calculated after identity-matrix solved if solvable
    ArrayList<Integer> E = null;
    ArrayList<Integer> F = null;
    ArrayList<Integer> F_Remainder = null;
    
    boolean Switch_aii;
    
    
    public TrueGF256_Berlekamp_Welch_algo_additiveinverseisitself(int gf_val, int[] temp_b, int temp_e_count, boolean switch_aii) throws Exception
    {
        Switch_aii = switch_aii;
        GF = gf_val;
        b = temp_b; //do not use b[0]
        recv_max_index = b.length - 1;
        
        e_count = temp_e_count;
        e = new int[e_count];
        
        q_max_index = recv_max_index - e_count - 1;
        q = new int[recv_max_index - e_count];// we will use q[0]
        
        Tools tool = new Tools();
        a = new int[b.length]; //DO NOT USE a[0]
        a[0] = Integer.MAX_VALUE; //DO NOT USE a[0] just try to mark it as not used
        for (int ii=1; ii <= recv_max_index; ii++)
        {
            if (Switch_aii)
            {
                int num = tool.Table_Exponent_Of_Alpha_To_Integer()[ii-1];
                a[ii] = num;
            }
            else
            {
                a[ii] = ii-1; //0,1,2,3,4,5,6,7,8 ..... ii-1
            }
        }
        
        //answer_matrix is rows==b.length x columns=1 represented by an array 
        //where answer_matrix[0] is not used so we will not be confused
        answer_matrix = new int[b.length];
        //it holds values -b[1]a[1]^2 .....-b[recv_max_index]a[recv_max_index]^2
        for (int ii=1; ii <= recv_max_index; ii++)
        {
            //GF(256)
            //no such thing as 2^x==0
            if (a[ii]==0 || b[ii]==0)
                answer_matrix[ii]=0;
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[ii]];
                //alpha^(alpha_exp) * alpha^(alpha_exp) ==
                //alpha^(alpha_exp + alpha_exp)
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[ii]];
                int total_exp = b_as_alpha_to_something_exp + (alpha_exp * e_count);
                if (total_exp >= 256)
                    total_exp %= 255;

                int num = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                answer_matrix[ii] = num;
            }
        }
        
        
        matrix = new int[b.length][b.length];
        //so do not use matrix[0][anything] or matrix[anything][0] just try to mark it as not used
        for (int ii=0; ii <= recv_max_index; ii++)
        {
            matrix[0][ii] = Integer.MAX_VALUE;
            matrix[ii][0] = Integer.MAX_VALUE;
        }
        
        
        //matrix = b[1]  b[1]a[1]  -1  -a[1]  -a[1]^2  -a[1]^3  -a[1]^4
        //         b[2]  b[2]a[2]  -1  -a[2]  -a[2]^2  -a[2]^3  -a[2]^4
        //         b[3]  b[3]a[3]  -1  -a[3]  -a[3]^2  -a[3]^3  -a[3]^4
        //         b[4]  b[4]a[4]  -1  -a[4]  -a[4]^2  -a[4]^3  -a[4]^4
        //         b[5]  b[5]a[5]  -1  -a[5]  -a[5]^2  -a[5]^3  -a[5]^4
        //         b[6]  b[6]a[6]  -1  -a[6]  -a[6]^2  -a[6]^3  -a[6]^4
        //         b[7]  b[7]a[7]  -1  -a[7]  -a[7]^2  -a[7]^3  -a[7]^4
        
        switch (e_count)
        {
            case 14: Fill_Matrix_14();
                     break;
            case 13: Fill_Matrix_13();
                     break;
            case 12: Fill_Matrix_12();
                     break;
            case 11: Fill_Matrix_11();
                     break;
            case 10: Fill_Matrix_10();
                     break;
            case 9: Fill_Matrix_9();
                    break;
            case 8: Fill_Matrix_8();
                    break;
            case 7: Fill_Matrix_7();
                    break;
            case 6: Fill_Matrix_6();
                    break;
            case 5: Fill_Matrix_5();
                    break;
            case 4: Fill_Matrix_4();
                    break;
            case 3: Fill_Matrix_3();
                    break;
            case 2: Fill_Matrix_2();
                    break;
            case 1: Fill_Matrix_1();
                    break;
            default: throw new Exception("constructor: unsupported e_count");
        }
    }//constructor
    
    
    private void Fill_Matrix_14()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            //e3
            //e4
            //e5
            //e6
            //e7
            //e8
            //e9
            //e10
            //e11
            //e12
            //e13
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
                matrix[row][4] = 0; //e3
                matrix[row][5] = 0; //e4
                matrix[row][6] = 0; //e5
                matrix[row][7] = 0; //e6
                matrix[row][8] = 0; //e7
                matrix[row][9] = 0; //e8
                matrix[row][10] = 0; //e9
                matrix[row][11] = 0; //e10
                matrix[row][12] = 0; //e11
                matrix[row][13] = 0; //e12
                matrix[row][14] = 0; //e13
            }   
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 2);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
                
                //e3
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 3);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][4] = gf256val;
                
                //e4
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 4);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][5] = gf256val;
                
                //e5
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 5);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][6] = gf256val;
                
                //e6
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 6);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][7] = gf256val;
                
                //e7
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 7);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][8] = gf256val;
                
                //e8
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 8);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][9] = gf256val;
                
                //e9
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 9);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][10] = gf256val;
                
                //e10
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 10);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][11] = gf256val;
                
                //e11
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 11);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][12] = gf256val;
                
                //e12
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 12);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][13] = gf256val;
                
                //e13
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 13);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][14] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][15] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][16] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][17] = 0;
                matrix[row][18] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][17] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][18] = temp;
            }


            for (int col=19; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 15);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_14
    
    
    private void Fill_Matrix_13()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            //e3
            //e4
            //e5
            //e6
            //e7
            //e8
            //e9
            //e10
            //e11
            //e12
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
                matrix[row][4] = 0; //e3
                matrix[row][5] = 0; //e4
                matrix[row][6] = 0; //e5
                matrix[row][7] = 0; //e6
                matrix[row][8] = 0; //e7
                matrix[row][9] = 0; //e8
                matrix[row][10] = 0; //e9
                matrix[row][11] = 0; //e10
                matrix[row][12] = 0; //e11
                matrix[row][13] = 0; //e12
            }   
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 2);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
                
                //e3
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 3);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][4] = gf256val;
                
                //e4
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 4);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][5] = gf256val;
                
                //e5
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 5);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][6] = gf256val;
                
                //e6
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 6);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][7] = gf256val;
                
                //e7
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 7);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][8] = gf256val;
                
                //e8
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 8);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][9] = gf256val;
                
                //e9
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 9);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][10] = gf256val;
                
                //e10
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 10);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][11] = gf256val;
                
                //e11
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 11);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][12] = gf256val;
                
                //e12
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 12);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][13] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][14] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][15] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][16] = 0;
                matrix[row][17] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][16] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][17] = temp;
            }


            for (int col=18; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 14);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_13
    
    
    private void Fill_Matrix_12()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            //e3
            //e4
            //e5
            //e6
            //e7
            //e8
            //e9
            //e10
            //e11
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
                matrix[row][4] = 0; //e3
                matrix[row][5] = 0; //e4
                matrix[row][6] = 0; //e5
                matrix[row][7] = 0; //e6
                matrix[row][8] = 0; //e7
                matrix[row][9] = 0; //e8
                matrix[row][10] = 0; //e9
                matrix[row][11] = 0; //e10
                matrix[row][12] = 0; //e11
            }   
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 2);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
                
                //e3
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 3);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][4] = gf256val;
                
                //e4
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 4);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][5] = gf256val;
                
                //e5
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 5);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][6] = gf256val;
                
                //e6
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 6);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][7] = gf256val;
                
                //e7
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 7);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][8] = gf256val;
                
                //e8
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 8);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][9] = gf256val;
                
                //e9
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 9);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][10] = gf256val;
                
                //e10
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 10);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][11] = gf256val;
                
                //e11
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 11);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][12] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][13] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][14] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][15] = 0;
                matrix[row][16] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][15] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][16] = temp;
            }


            for (int col=17; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 13);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_12
    
    
    private void Fill_Matrix_11()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            //e3
            //e4
            //e5
            //e6
            //e7
            //e8
            //e9
            //e10
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
                matrix[row][4] = 0; //e3
                matrix[row][5] = 0; //e4
                matrix[row][6] = 0; //e5
                matrix[row][7] = 0; //e6
                matrix[row][8] = 0; //e7
                matrix[row][9] = 0; //e8
                matrix[row][10] = 0; //e9
                matrix[row][11] = 0; //e10
            }   
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 2);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
                
                //e3
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 3);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][4] = gf256val;
                
                //e4
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 4);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][5] = gf256val;
                
                //e5
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 5);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][6] = gf256val;
                
                //e6
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 6);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][7] = gf256val;
                
                //e7
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 7);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][8] = gf256val;
                
                //e8
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 8);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][9] = gf256val;
                
                //e9
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 9);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][10] = gf256val;
                
                //e10
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 10);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][11] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][12] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][13] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][14] = 0;
                matrix[row][15] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][14] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][15] = temp;
            }


            for (int col=16; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 12);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_11
    
    
    private void Fill_Matrix_10()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            //e3
            //e4
            //e5
            //e6
            //e7
            //e8
            //e9
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
                matrix[row][4] = 0; //e3
                matrix[row][5] = 0; //e4
                matrix[row][6] = 0; //e5
                matrix[row][7] = 0; //e6
                matrix[row][8] = 0; //e7
                matrix[row][9] = 0; //e8
                matrix[row][10] = 0; //e9
            }   
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 2);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
                
                //e3
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 3);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][4] = gf256val;
                
                //e4
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 4);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][5] = gf256val;
                
                //e5
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 5);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][6] = gf256val;
                
                //e6
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 6);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][7] = gf256val;
                
                //e7
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 7);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][8] = gf256val;
                
                //e8
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 8);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][9] = gf256val;
                
                //e9
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 9);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][10] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][11] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][12] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][13] = 0;
                matrix[row][14] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][13] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][14] = temp;
            }


            for (int col=15; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 11);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_10
    
    
    private void Fill_Matrix_9()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            //e3
            //e4
            //e5
            //e6
            //e7
            //e8
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
                matrix[row][4] = 0; //e3
                matrix[row][5] = 0; //e4
                matrix[row][6] = 0; //e5
                matrix[row][7] = 0; //e6
                matrix[row][8] = 0; //e7
                matrix[row][9] = 0; //e8
            }   
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 2);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
                
                //e3
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 3);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][4] = gf256val;
                
                //e4
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 4);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][5] = gf256val;
                
                //e5
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 5);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][6] = gf256val;
                
                //e6
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 6);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][7] = gf256val;
                
                //e7
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 7);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][8] = gf256val;
                
                //e8
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 8);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][9] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][10] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][11] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][12] = 0;
                matrix[row][13] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][12] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][13] = temp;
            }


            for (int col=14; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 10);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_9
    
    
    private void Fill_Matrix_8()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            //e3
            //e4
            //e5
            //e6
            //e7
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
                matrix[row][4] = 0; //e3
                matrix[row][5] = 0; //e4
                matrix[row][6] = 0; //e5
                matrix[row][7] = 0; //e6
                matrix[row][8] = 0; //e7
            }   
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 2);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
                
                //e3
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 3);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][4] = gf256val;
                
                //e4
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 4);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][5] = gf256val;
                
                //e5
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 5);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][6] = gf256val;
                
                //e6
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 6);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][7] = gf256val;
                
                //e7
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 7);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][8] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][9] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][10] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][11] = 0;
                matrix[row][12] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][11] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][12] = temp;
            }


            for (int col=13; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 9);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_8
    
    
    private void Fill_Matrix_7()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            //e3
            //e4
            //e5
            //e6
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
                matrix[row][4] = 0; //e3
                matrix[row][5] = 0; //e4
                matrix[row][6] = 0; //e5
                matrix[row][7] = 0; //e6
            }   
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 2);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
                
                //e3
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 3);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][4] = gf256val;
                
                //e4
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 4);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][5] = gf256val;
                
                //e5
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 5);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][6] = gf256val;
                
                //e6
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 6);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][7] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][8] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][9] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][10] = 0;
                matrix[row][11] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][10] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][11] = temp;
            }


            for (int col=12; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 8);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_7
    
    
    private void Fill_Matrix_6()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            //e3
            //e4
            //e5
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
                matrix[row][4] = 0; //e3
                matrix[row][5] = 0; //e4
                matrix[row][6] = 0; //e5
            }   
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + alpha_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
                
                //e3
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 3);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][4] = gf256val;
                
                //e4
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 4);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][5] = gf256val;
                
                //e5
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 5);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][6] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][7] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][8] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][9] = 0;
                matrix[row][10] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][9] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][10] = temp;
            }


            for (int col=11; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 7);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_6
    
    
    private void Fill_Matrix_5()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            //e3
            //e4
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
                matrix[row][4] = 0; //e3
                matrix[row][5] = 0; //e4
            }
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + alpha_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
                
                //e3
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 3);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][4] = gf256val;
                
                //e4
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 4);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][5] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][6] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][7] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][8] = 0;
                matrix[row][9] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][8] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][9] = temp;
            }


            for (int col=10; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 6);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_5
    
    
    private void Fill_Matrix_4()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            //e3
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
                matrix[row][4] = 0; //e3
            }
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + alpha_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
                
                //e3
                total_exp = b_as_alpha_to_something_exp + (alpha_exp * 3);
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][4] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][5] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][6] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][7] = 0;
                matrix[row][8] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][7] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][8] = temp;
            }


            for (int col=9; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 5);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_4
    
    
    private void Fill_Matrix_3()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            
            //e0
            matrix[row][1] = b[row];

            //e1
            //e2
            if (b[row]==0 || a[row]==0)
            {
                matrix[row][2] = 0; //e1
                matrix[row][3] = 0; //e2
            }
            else
            {
                //e1
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = b_as_alpha_to_something_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
                
                //e2
                total_exp = b_as_alpha_to_something_exp + alpha_exp + alpha_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][3] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q0
            matrix[row][4] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            //q1
            matrix[row][5] = a[row];

            //q2 and q3
            if (a[row]==0)
            {
                matrix[row][6] = 0;
                matrix[row][7] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][6] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][7] = temp;
            }


            for (int col=8; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 4);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_3
    
    
    private void Fill_Matrix_2()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            //e0
            matrix[row][1] = b[row];

            //e1
            if (b[row]==0 || a[row]==0)
                matrix[row][2] = 0;
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int b_as_alpha_to_something_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[b[row]];
                int total_exp = alpha_exp + b_as_alpha_to_something_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256val = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                matrix[row][2] = gf256val;
            }

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            matrix[row][3] = 1;

            //in true GF(256) additive inverse is itself which is why XOR for addition works
            matrix[row][4] = a[row];

            if (a[row]==0)
            {
                matrix[row][5] = 0;
                matrix[row][6] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][5] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][6] = temp;
            }


            for (int col=7; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 3);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_2
    
    
    private void Fill_Matrix_1()
    {
        Tools tool = new Tools();
        for (int row=1; row <= recv_max_index; row++)
        {
            //GF(256)
            //e0
            matrix[row][1] = b[row];

            //q0
            //in true GF(256) additive inverse is itself which is why XOR for addition works
            matrix[row][2] = 1;

            //q1
            //in true GF(256) additive inverse is itself which is why XOR for addition works
            matrix[row][3] = a[row];

            //q4
            //q5
            if (a[row]==0)
            {
                matrix[row][4] = 0;
                matrix[row][5] = 0;
            }
            else
            {
                int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                int total_exp = alpha_exp * 2;
                if (total_exp >= 256)
                    total_exp %= 255;
                int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][4] = temp;


                total_exp = alpha_exp * 3;
                if (total_exp >= 256)
                    total_exp %= 255;
                temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                matrix[row][5] = temp;
            }


            for (int col=6; col <= recv_max_index; col++)
            {
                if (a[row] == 0)
                {
                    matrix[row][col] = 0;
                }
                else
                {
                    int alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[a[row]];
                    int total_exp = alpha_exp * (col - 2);
                    if (total_exp >= 256)
                         total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row][col] = temp;
                }
            }
        }//for row
    }//Fill_Matrix_1
    
    
    public void Debug_Print()
    {
        System.out.println();
        for (int row=1; row <= recv_max_index; row++)
        {
            System.out.print(" |");
            for (int col=1; col <= recv_max_index; col++)
            {
                String tempstr = Integer.toString(matrix[row][col]);
                if (tempstr.length() == 3)
                    tempstr = "  " + tempstr;
                else if (tempstr.length() == 2)
                    tempstr = "   " + tempstr;
                else if (tempstr.length() == 1)
                    tempstr = "    " + tempstr;
                System.out.print(tempstr);
            }
            
            System.out.printf("|   |%3d|", answer_matrix[row]);
            System.out.println();
        }
    }//Debug_Print
    
    
    public void GF256_MultiplyAndAdd(int row_mult, int mult, int row_add) throws Exception
    {
        Tools tool = new Tools();
        if (GF == 256)
        {
            for (int col=1; col <= recv_max_index; col++)
            {
                int temp = matrix[row_mult][col];
                if (temp != 0)
                {
                    int temp_exp  = tool.Table_Integer_To_Exponent_Of_Alpha()[temp];
                    int mult_exp  = tool.Table_Integer_To_Exponent_Of_Alpha()[mult];
                    int total_exp = temp_exp + mult_exp;
                    if (total_exp >= 256)
                        total_exp %= 255;
                    int gf256num = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //in true GF(256) additive inverse is itself which is why XOR for addition works
                    matrix[row_add][col] ^= gf256num;
                }
                //else matrix[row_add][col] remains the same
            }

            int temp = answer_matrix[row_mult];
            if (temp != 0)
            {
                int temp_exp  = tool.Table_Integer_To_Exponent_Of_Alpha()[temp];
                int mult_exp  = tool.Table_Integer_To_Exponent_Of_Alpha()[mult];
                int total_exp = temp_exp + mult_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256num = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                //in true GF(256) additive inverse is itself which is why XOR for addition works
                answer_matrix[row_add] ^= gf256num;
            }
            //else answer_matrix[row_add] remains the same
        }
        else
            throw new Exception("GF256_MultiplyAndAdd only for GF(256)!!!!!");
    }//GF256_MultiplyAndAdd
    
    
    public void GF256_Multiply_Row(int row, int mult) throws Exception
    {
        Tools tool = new Tools();
        if (GF == 256)
        {
            for (int col=1; col <= recv_max_index; col++)
            {
                int temp = matrix[row][col];
                if (temp != 0)
                {
                    int temp_exp  = tool.Table_Integer_To_Exponent_Of_Alpha()[temp];
                    int mult_exp  = tool.Table_Integer_To_Exponent_Of_Alpha()[mult];
                    int total_exp = temp_exp + mult_exp;
                    if (total_exp >= 256)
                        total_exp %= 255;
                    int gf256num = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    matrix[row][col] = gf256num;
                }
                //else matrix[row][col] remains 0
            }

            int temp = answer_matrix[row];
            if (temp != 0)
            {
                int temp_exp  = tool.Table_Integer_To_Exponent_Of_Alpha()[temp];
                int mult_exp  = tool.Table_Integer_To_Exponent_Of_Alpha()[mult];
                int total_exp = temp_exp + mult_exp;
                if (total_exp >= 256)
                    total_exp %= 255;
                int gf256num = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                answer_matrix[row] = gf256num;
            }
            //else answer_matrix[row] remains 0
        }
        else //GF(256)
            throw new Exception("GF256_Multiply_Row only used for GF(256)!!!!!");
    }//GF256_Multiply_Row
    
    
    public void GF256_Add_RowA_By_RowB(int RowA, int RowB) throws Exception
    {
        if (GF == 256)
        {
            for (int col=1; col <= recv_max_index; col++)
            {
                int Aval = matrix[RowA][col];
                int Bval = matrix[RowB][col];
                int temp = Aval ^ Bval;
                matrix[RowA][col] = temp;
            }

            int Aval = answer_matrix[RowA];
            int Bval = answer_matrix[RowB];
            int temp = Aval ^ Bval;
            answer_matrix[RowA] = temp;
        }
        else
            throw new Exception("GF256_Add_RowA_By_RowB only for GF(256)!!!!!");
    }//GF256_Add_RowA_By_RowB
    
    
    public void Manual_Console_Solve()
    {   
        Scanner scanner = new Scanner(System.in);
           
        String str = "";
        while (str.equals("x") == false)
        {
            //Debug_Print();
            
            System.out.println();
            System.out.println("Enter m,x,y: Multiply row-x by y (GF7)");
            System.out.println("Enter a,x,y: Add row-x by row-y (GF7)");
            System.out.println("Enter ma,x,y,z: Multiply row-x by y then add it to row-z(GF7)");
            System.out.println("Enter x: Exit");
            System.out.println("Enter:");
            
            str = scanner.nextLine();
            
            try
            {
                String[] temp = str.split(",");
                if (temp[0].equals("m"))
                    GF256_Multiply_Row(Integer.parseInt(temp[1]), Integer.parseInt(temp[2]));
                else if (temp[0].equals("a"))
                    GF256_Add_RowA_By_RowB(Integer.parseInt(temp[1]), Integer.parseInt(temp[2]));
                else if (temp[0].equals("ma"))
                    GF256_MultiplyAndAdd(Integer.parseInt(temp[1]), Integer.parseInt(temp[2]), Integer.parseInt(temp[3]));
                else if (temp[0].equals("x") == false)
                {
                    System.out.println("BAD STRING ENTERED!!!!!");
                    System.out.println();
                }
            }
            catch (Exception ex)
            {
                ex.printStackTrace();
                return;
            }
        }//while
    }//Manual_Console_Solve
    
    private boolean Check_Identity_Matrix()
    {
        for (int row=1; row <= recv_max_index; row++)
        {
            for (int col=1; col <= recv_max_index; col++)
            {
                if (matrix[row][col] != 1 && row==col)
                    return false;
                else if (matrix[row][col] != 0 && row!=col)
                    return false;
            }
        }
        return true;
    }//Check_Identity_Matrix    
    
    public boolean Robot_Solve()
    {   
        int ii=1;
        //Debug_Print();
        while (true)
        {   
            boolean good = false;
            if (GF==256)
            {
                try
                {
                    good = GF256_Make_Column_Good_1_And_0s(ii);
                    //Debug_Print();
                }
                catch (Exception ex)
                {
                    ex.printStackTrace();
                    return false;
                }
                
            }
            else
            {
                System.out.println("Robot_Solve: only for GF(256)!!!!!");
                return false;
            }
            if (!good)
                return false;
            
            ii++;
            if (ii==(recv_max_index + 1))
                break;
        }
        
        boolean tempbool = Check_Identity_Matrix();
        return tempbool;
    }//Robot_Solve
    
    
    private boolean GF256_Make_Column_Good_1_And_0s(int ii)
    {
        Tools tool = new Tools();
        boolean good = false;
        for (int row=1; row <= recv_max_index; row++)
        {
            if (row==ii && matrix[row][ii] == 1)
                good = true;
            else if (row!=ii && matrix[row][ii] == 0)
                good = true;
            else
            {
                good = false;
                break;
            }
        }
        
        if (good) //this column is good already
            return true;
        
        //if column all 0s, can not make it good
        int zero_count=0;
        for (int row=1; row <= recv_max_index; row++)
        {
            if (matrix[row][ii] == 0)
                zero_count++;
        }
        if (zero_count == recv_max_index)
            return false; //can not make it good
        
        
        //make matrix[ii][ii]==1
        if (matrix[ii][ii] != 1)
        {
            if (matrix[ii][ii] != 0)
            {
                int temp_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[matrix[ii][ii]];
                int mult_exp = 255 - temp_exp;
                int mult = tool.Table_Exponent_Of_Alpha_To_Integer()[mult_exp];
                try
                {
                    this.GF256_Multiply_Row(ii, mult);
                }
                catch (Exception ex)
                {
                    ex.printStackTrace();
                    return false;
                }
            }
            else //matrix[ii][ii] == 0
            {
                //find another row after this row==ii
                //whose matrix[row][ii] != 0 and row <= recv_max_index
                int row = ii;
                while (matrix[row][ii] == 0 && row <= recv_max_index)
                {
                    row++;
                    if (row > recv_max_index)
                        return false;
                }
                
                try
                {
                    if (matrix[row][ii] == 1)
                        this.GF256_Add_RowA_By_RowB(ii, row);
                    else //non-0 and non-1
                    {
                        int temp_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[matrix[row][ii]];
                        int mult_exp = 255 - temp_exp;
                        int mult     = tool.Table_Exponent_Of_Alpha_To_Integer()[mult_exp];
                        this.GF256_MultiplyAndAdd(row, mult, ii);
                    }//non-0 and non-1
                }
                catch (Exception ex)
                {
                    ex.printStackTrace();
                    return false;
                }
            }//matrix[ii][ii] == 0
        }//make matrix[ii][ii]==1
        
        //ok matrix[ii][ii]==1
        //make cell above or below diag to 0
        for (int row=1; row <= recv_max_index; row++)
        {
            if (row != ii)
            {
                if (matrix[row][ii] != 0)
                {
                    //make it 0
                    //mult diag-row and then XOR because XOR is additive-inverse
                    int temp = matrix[row][ii];
                    try
                    {
                        this.GF256_MultiplyAndAdd(ii, temp, row);
                    }
                    catch (Exception ex)
                    {
                        ex.printStackTrace();
                        return false;
                    }
                }
                //else already 0
            }
        }
        
        return true;
    }//GF256_Make_Column_Good_1_And_0s
    
    
    public void Fill_QuestionMatrix_With_AnswerMatrix()
    {
        int ii=1; //answer_matrix[0] is not used
        
        for (int err=0; err < e_count; err++)
        {
            e[err] = answer_matrix[ii];
            ii++;
        }
        for (int qii=0; qii < q.length; qii++)
        {
            q[qii] = answer_matrix[ii];
            ii++;
        }
    }//Fill_QuestionMatrix_With_AnswerMatrix
    
    public void Debug_Print_Q_And_E_Functions()
    {
        
        System.out.print("Q(ai) = ");
        for (int qii=q_max_index; qii >= 0; qii--) //q uses q[0], q[0] is coefficient for x^0
        {
            System.out.print(q[qii] + "x^" + qii);
            if (qii != 0)
                System.out.print(" + ");
        }
        
        System.out.println();
        System.out.print("E(ai) = ");
        for (int eii=e_count; eii >= 0; eii--) //e uses e[0], e[0] is coefficient for x^0, e[e_count] is constrained to 1 and is the coefficient for x^e_count
        {
            if (eii == e_count)
                System.out.print("1x^" + eii + " + ");
            else
            {
                System.out.print(e[eii] + "x^" + eii);
                if (eii != 0)
                    System.out.print(" + ");
            }
        }
        System.out.println();
    }//Debug_Print_Q_And_E_Functions
    
    
    public void GF_Polynomial_Long_Division_To_Find_F_Function()
    {
        //q uses q[0], q[0] is coefficient for x^0
        //e uses e[0], e[0] is coefficient for x^0, e[e_count] is constrained to 1 and is the coefficient for x^e_count
        //F(ai) = Q(ai) / E(ai)
        //
        //example:
        //Q(ai) = 3x^4 + 1x^3 + 3x^2 + 3x + 4
        //E(ai) = 1x^2 + 2x + 4
        //
        //
        //F(ai) = 3x^2  2x    1
        //        --------------------------------
        //        3x^4  1x^3  3x^2  3x  4
        //        3x^4  6x^3  12x^2
        //        --------------------------------
        //              2x^3  5x^2  3x  4
        //              2x^3  4x^2  8x
        //              --------------------------
        //                    1x^2  2x  4
        //                    1x^2  2x  4
        //                    --------------------
        //                          
        E           = new ArrayList<Integer>();
        F           = new ArrayList<Integer>();
        F_Remainder = new ArrayList<Integer>();
        
        Tools tool = new Tools();
        
        ArrayList<Integer> dividend = new ArrayList<Integer>();
        
        for (int qii=q_max_index; qii >= 0; qii--) //q uses q[0], q[0] is coefficient for x^0
            dividend.add(q[qii]);
        
        for (int eii=e_count; eii >= 0; eii--) //e uses e[0], e[0] is coefficient for x^0, e[e_count] is constrained to 1 and is the coefficient for x^e_count
        {
            if (eii == e_count)
                E.add(1);
            else
                E.add(e[eii]);
        }
        
        int E_term_length = e_count + 1;
        
        while (dividend.size() >= E_term_length)
        {
            int mult = 0;
            int mult_alphaexp = 0;
            for (int ii=0; ii < E.size(); ii++)
            {
                if (ii==0)
                {
                    mult = dividend.get(ii);
                    F.add(mult);
                    dividend.set(ii, 0);
                }
                else
                {
                    if (mult != 0)
                    {
                        mult_alphaexp = tool.Table_Integer_To_Exponent_Of_Alpha()[mult];
                        int temp_dividend = dividend.get(ii);                    
                        int temp_E = E.get(ii);
                        if (temp_E != 0)
                        {
                            int temp_E_alphaexp = tool.Table_Integer_To_Exponent_Of_Alpha()[temp_E];
                            int total_alphaexp = temp_E_alphaexp + mult_alphaexp;
                            if (total_alphaexp >= 256)
                                total_alphaexp %= 255;

                            int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_alphaexp];
                            int xornum  = temp_dividend ^ temp;
                            dividend.set(ii, xornum);
                        }
                    }
                    //else dividend[ii] remains same
                }       
            }//for E
            
            dividend.remove(0);
        }//while
        
        System.out.println();
        System.out.print("F(ai) = ");
        for (int ii=0; ii < F.size(); ii++)
        {
            System.out.print(F.get(ii) + "x^" + (F.size() - ii - 1));
            if (ii != (F.size() - 1))
                 System.out.print(" + ");
        }
        System.out.println();
        System.out.print("Remainder = ");
        for (int ii=0; ii < dividend.size(); ii++)
        {
            F_Remainder.add(dividend.get(ii));
            System.out.print(dividend.get(ii) + "x^" + (dividend.size() - ii - 1));
            if (ii != (dividend.size() - 1))
                System.out.print(" + ");
        }
    }//GF_Polynomial_Long_Division_To_Find_F_Function
    
    public boolean Check_Remainder_All_0()
    {
        System.out.println();
        System.out.print("a[ii] = ");
        for (int ii=1; ii < this.a.length; ii++) //DO NOT USE a[0]
            System.out.print(a[ii] + " "); //DO NOT USE a[0]
        System.out.println();
        
        for (int ii=0; ii < F_Remainder.size(); ii++)
        {
            if (F_Remainder.get(ii) != 0)
                return false;
        }
        return true;
    }//Check_Remainder_All_0
    
    
    public ArrayList<Integer> GF256_Find_Error_Locations() throws Exception
    {
        if (GF != 256)
            throw new Exception("GF256_Find_Error_Locations: Only for GF(256)");
        
        Tools tool = new Tools();
        
        //error location found if E(ai)==0 and i==location of error in received-bytes
        //and correction is F(ai)
        ArrayList<Integer> loc = new ArrayList<Integer>();
        
        for (int ii=1; ii <= recv_max_index; ii++)
        {
            //calculate E(ai)
            int sum=0;
            for (int eterm=0; eterm < E.size(); eterm++)
            {
                int coeff = E.get(eterm);
                
                int x_exp = E.size() - eterm - 1;
                int x = a[ii];
                
                //a(i) = 2^0, 2^1, 2^2, 2^3 ..... 2^(i-1)
                //so you can use tables
                //that means x is a power-of-2 in GF(256)
                if (x!=0 && coeff!=0)
                {
                    int x_as_alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[x];
                    int coeff_as_alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[coeff];
                    int total_exp = (x_as_alpha_exp * x_exp) + coeff_as_alpha_exp;
                    if (total_exp >= 256)
                        total_exp %= 255;
                    int temp = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    //because additive-inverse in GF(256) is itself so XOR is addition/subtraction in GF(256)
                    sum ^= temp;
                }
            }
            
            if (sum==0) //found error location == ii
                loc.add(ii);
        }//ai loop
        
        return loc;
    }//GF256_Find_Error_Locations
    
    
    /*****
    private void later()
    {
                //find correction
                sum = 0;
                for (int f_term=0; f_term < F.size(); f_term++)
                {
                    int coeff = F.get(f_term);
                    int x_exp = F.size() - f_term - 1;
                    int x = a[ii];
                    if (GF == 7)
                    {
                          int temp = coeff * (int)Math.pow(x, x_exp);
                          sum += temp;
                    }
                    else if (GF == 256)
                    {
                          //make sure a(i) = 2^0, 2^1, 2^2, 2^3 ..... 2^(i-1)
                          //so you can use tables
                          //that means x==2
                          throw new Exception("Find_Error_Locations_And_Corrections: codearea2: did not finish code for GF(256)");
                    }
                    else
                    {
                          throw new Exception("Find_Error_Locations_And_Corrections: codearea2: did not finish code for unknown GF because Math.pow may get too big");
                    }
                }
                int correction = GF_value(sum);
                locandcorrect.add(ii + "," + correction);
    }
    *****/
    
    
}//class