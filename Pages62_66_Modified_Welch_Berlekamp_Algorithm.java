


/***********************************************
Modified-Welch-Berlekamp
see arizona-paper

l == 8 , 8 bits per byte
systematic encoding == received-bytes is message-bytes followed by E.C-bytes
n == total number of bytes     (ex. 28 QR-code version 6H)
k == number of message bytes   (ex. 15 QR-code version 6H)

~r == received-bytes vector
~c == correct-bytes vector

f = number of erased-bytes among the E.C-byte-indices
F = indices that were erased as a set

mj = µ(αj) for 0 ≤ j < k 
 
u(x) = s(x) / q(x)

(so if you do not erase some E.C-bytes then tˆ==6 for QRcode version 6H
 so if you erase 1 E.C-bytes then tˆ==6 for QRcode version 6H
 so if you erase 2 E.C-bytes then tˆ==5 for QRcode version 6H
 so if you erase 3 E.C-bytes then tˆ==5 for QRcode version 6H
 so if you erase 4 E.C-bytes then tˆ==4 for QRcode version 6H
 so if you erase 5 E.C-bytes then tˆ==4 for QRcode version 6H
 so if you erase 6 E.C-bytes then tˆ==3 for QRcode version 6H
 so if you erase 7 E.C-bytes then tˆ==3 for QRcode version 6H
 so if you erase 8 E.C-bytes then tˆ==2 for QRcode version 6H
 so if you erase 9 E.C-bytes then tˆ==2 for QRcode version 6H
 so if you erase 10 E.C-bytes then tˆ==1 for QRcode version 6H
 so if you erase 11 E.C-bytes then tˆ==1 for QRcode version 6H
 for QRcode version 6H you can not erase 12 or more E.C.bytes to use this algorithm)
 *****************************************************************************************/

/****************************************************************************************
4.2 Error Detection
(see arizona-paper)

Example 4.2.1. 
Let l = 8, k = 5, n = 8, ~α = (0, 1, 2, 3, 4, 5, 6, 7)T
and
~m = (233, 211, 0, 7, 18)T

Encoding systematically yields
~c = (233, 211, 0, 7, 18, 166, 14, 135)T

~c is then transmitted and an error occurs resulting in a received word
~r = (233, 117, 0, 7, 18, 166, 14, 135)T

Since this encoding is systematic, we can recover a message from ~r by simply selecting the
first k characters, so we have
~m' = (233, 117, 0, 7, 18)T

Encoding ~m' yeilds
~c' = (233, 117, 0, 7, 18, 243, 87, 45)T

Clearly ~r != ~c
so we correctly conclude that ~r != ~c and an error has occurred.

Example 4.2.2. This example illustrates how error detection can fail. 
Let l = 8, k = 5, n = 6, ~α = (0, 1, 2, 3, 4, 5)T
~m = (233, 211, 0, 7, 18)T

Encoding systematically yields
~c = (233, 211, 0, 7, 18, 166)T

~c is then transmitted and two errors occur resulting in a received word
~r = (233, 117, 0, 7, 18, 243)T

Since this encoding is systematic, we can recover a message from ~r by simply selecting the
first k characters, so we have
~m' = (233, 117, 0, 7, 18)T

Encoding ~m' yeilds
~c' = (233, 117, 0, 7, 18, 243)T

~r == ~c'
so we incorrectly conclude that ~r == ~c and there have been no errors.

4.3.1 Distance
The notion of distance between words is important in coding theory, and will provide an
appropriate tool to determine when errors are detectable and correctable. In this section
we define the Hamming distance, and prove an important result about the distance between
two codewords in a Reed-Solomon code. For generality, we assume that all arithmetic is
performed over some GF, although in applications in this paper we will always have
GF(2^l) (where l==8)

Definition 4.3.1. Let ~a and ~b be vectors in GF(256)
Then the Hamming distance or just distance between ~a and ~b is
d(~a,~b)
(see arizona-paper for exact series-equation)
d(~a,~b) is the count of the indices where ~a and ~b differ.

Example 4.3.3. The Hamming distance between the vectors 
~a = (3, 7, 2, 4, 5) and ~b = (3, 1, 2, 8, 5) is 2 
since they differ in 2 locations.

Remark 4.3.4. Let ~c be a codeword is a Reed-Solomon code, and ~r be a received word.
Then d(~c, ~r) gives the number of errors in the received word.

Definition 4.3.5. The distance d of a code is the minimum distance between any 2
codewords in the code C:
d = min( d(~c1, ~c2) ) where  ~c1,~c2 ∈ C

Theorem 4.3.6. Let ~c1 and ~c2 be two codewords in a Reed-Solomon code. 
Then d(~c1, ~c2) ≥ n − k + 1 and so Reed-Solomon codes have distance d = n − k + 1.
(see arizona-paper for proof)


Remark 4.3.8. n − k + 1 errors must occur to transform one codeword into another.
Therefore, we are guaranteed to be able to detect n − k or fewer errors in a received word.
If n−k+1 or more errors occurred, they might go undetected. We make a base assumption
that the probability of an error affecting any particular character in a codeword is smaller
than one half, and therefore it is more probable that a small number of errors accrue in a
received word than that a large number of errors have accrued. By making n sufficiently
larger than k, we can design the code so that an undetected error is a sufficiently low
probability event. Of course, this requires greater redundancy in the code which makes the
code less efficient.

Definition 4.3.9. A code of length n (number of characters in a codeword), 
dimension k (number of characters in the message), and distance d is called an (n, k, d) code.

Corollary 4.3.10. Reed-Solomon codes are (n, k, n − k + 1) codes.

Example 4.3.18. The Hamming distance between the vectors ~a = (3, 7, 2, 4, 5) and
~b = (3, 1, 2, 8, 5) restricted by I = {0, 1, 2} is dI(~a,~b) = 1 
since they differ in 1 of the first 3 locations.

4.3.2 (i think typo should be 4.3.19) Assumptions and Conditions for Error Correction
Given a received word ~r in which errors have been detected and no outside information
about the nature of the errors accumulated, it is impossible to recover the intended message
~m with absolute certainty. Two different codewords could be transformed into the same
received word through the accumulation of various errors, and by examining only the
received word, it is impossible to tell which of the codewords we started with. Error
correction is necessarily an uncertain business and the best we can do is to seek to recover ~m
with high probability. To do so, we always assume that the codeword closest in Hamming
distance to the received word is the correct codeword. This is equivalent to choosing
the codeword corresponding the fewest possible accumulated errors in the received word.
This assumption is called maximum likelihood decoding. The problem of decoding (error
correction) will be to identify the closest codeword to a given received word as efficiently as
possible. Maximum likelihood decoding will occasionally result in an incorrect decoding. In
this section we derive conditions under which maximum likelihood decoding is guaranteed
to succeed and see that success can be achieved with high probability.

Theorem 4.3.19. (i think typo should be 4.3.20)
(Error bound) In a Reed-Solomon code, if t errors have accumulated in
the received word ~r, and t < d/2,
then the closest codeword to ~r is the correct codeword.
Proof. Let ~c be the original codeword and ~r be the received word in which t errors have
accumulated. Then d(~c, ~r) = t. In search of contradiction, assume there exits some other
codeword ~c'
that is as close as or closer to ~r than ~c. 
Then 
d(~r, ~c') ≤ t and 
d(~c, ~c') ≤ d(~c, ~r) + d(~r, ~c') ≤ 2t < d. 
This violates the distance bound for Reed-Solomon codes, thus no such ~c' may exist.
If t ≥ d/2
errors are present, then some other codeword may be closer to or as close to
the received word as the original codeword. If so, the decoding may be incorrect.
It is also sometimes possible to recover a message from a received word in the presence
of both erasures (in known locations) and errors (in unknown locations). In this case we
seek the closest codeword to the received word not counting errors in the known erasure
locations.

Corollary 4.3.20. (i think typo should be 4.3.21)
(Error bound with erasures) In a Reed-Solomon, if f erasures in the
known locations F ⊂ {0, 1, . . . , n−1} and t errors (in unknown locations) have accumulated
in the received word ~r, and 2t+f < d, then the closest codeword to ~r in Hamming distance
restricted by Fcomplement is the correct codeword. 
Here Fcomplement is the complement of F: 
Fcomplement = {0, 1, . . . n − 1} \ F.


Chapter 5  Error Correction
The most challenging aspect to error correcting codes is efficient error correction when the
errors are in unknown locations. In the remainder of this dissertation, we explore a few
algorithms for error correction for Reed-Solomon codes. These error correction algorithms
all follow a common three-phase structure:
• In phase 1, a “key equation” is derived based on properties of the code. A particular
instance of this key equation depends on the received word, and the various encoding
parameters. The solution of the key equation will give information about the most
probable codeword corresponding to the given received word, the locations of the
errors, and/or their magnitudes.
• In phase 2, a solution to this key equation is computed.
• In phase 3, the solution of the key equation is used to recover the original message
or codeword.
In the following chapters, I will present multiple approaches to each of these phases, and
compare the efficiency of each.

5.1 Inputs, Assumptions, and Notation for ECC Algorithms
Each of these algorithms will accept the following inputs:
• the parameters, k, n, and ~α,
• a received word ~r in which errors have been detected (Section 4.2),
• the number of erasures f in ~r along with the set of indices F at which those erasures
occurred,
• the type of encoding (non-systematic, or systematic),
• (optional) the generating matrix for the code.
l, the size of a character (l==8==8bitsperbyte), 
will be treated as a system wide constant, and it will often be
assumed that the encoding was systematic. 
These inputs are subject to the constraints
k < n ≤ 2^l, f < n − k, |F| = f, |~α| = n,

( f = number of erased-bytes among the E.C-byte-indices
  F = indices that were erased as a set )

and all αi unique. Given these inputs, define the following intermediary items:

• r = n − k, the number of redundant characters in a codeword,
(not to be confused with ~r==received-bytes as vector, i do not think i should have pasted this but i did,
it is not necessary)
• d = r + 1 = n − k + 1, the distance of the code, 
(if n==28 and k==15 then d==14)
• W = Fcomplement = {0, 1, . . . , n − 1} \ F, the set of indices of characters in ~r that were not
erased,
(W is the message-byte-indices, you only try to erase in the E.C-byte-indices)
• tˆ= floor( (n−k−f)/2 ), the maximum number of errors within the error bound
(so if you do not erase some E.C-bytes then tˆ==6 for QRcode version 6H
 so if you erase 1 E.C-bytes then tˆ==6 for QRcode version 6H
 so if you erase 2 E.C-bytes then tˆ==5 for QRcode version 6H
 so if you erase 3 E.C-bytes then tˆ==5 for QRcode version 6H
 so if you erase 4 E.C-bytes then tˆ==4 for QRcode version 6H
 so if you erase 5 E.C-bytes then tˆ==4 for QRcode version 6H
 so if you erase 6 E.C-bytes then tˆ==3 for QRcode version 6H
 so if you erase 7 E.C-bytes then tˆ==3 for QRcode version 6H
 so if you erase 8 E.C-bytes then tˆ==2 for QRcode version 6H
 so if you erase 9 E.C-bytes then tˆ==2 for QRcode version 6H
 so if you erase 10 E.C-bytes then tˆ==1 for QRcode version 6H
 so if you erase 11 E.C-bytes then tˆ==1 for QRcode version 6H
 for QRcode version 6H you can not erase 12 or more E.C.bytes to use this algorithm)

In addition, in deriving the algorithms it will be useful to have notation for the following
unknown items:
• ~c, the unique codeword such that dW(~r, ~c) ≤ tˆ,
(the correct-bytes , dW means distance between ~r and ~c looking at only the message-byte-indices)
• µ(x), the message polynomial under systematic encoding, µ(~α) = ~c, deg(µ) = k − 1,
• ~m, the message corresponding to the codeword ~c under the given encoding,
• t, the number of errors present in ~r, t ≤ tˆ,
• δ = tˆ− t
• ~e, the error vector, ~r = ~c + ~e, (in GF(256) world)
• E, the set of indices of errors, E = {i : ei != 0}, |E| = t,
• g(x), the error locator polynomial, (see arizona-paper for equation)

We assume that t ≤ tˆ and therefore the number of errors is within the error bounds
of Section 4.3.2. Under this assumption, maximum likelihood decoding is guaranteed to
succeed. We may think of the problem of decoding equivalently as identifying ~c, ~m, m(x),
~e, E, or g(x) as given any one of these items we can easily calculate all the others.

5.2 The Modified Welch-Berlekamp Algorithm
We start our discussion of error correcting codes by considering a conceptually simple
decoding algorithm with suboptimal computational efficiency. This algorithm is based on
polynomial interpolation, and the presentation follows the work of Gemmel and Sudan in
[9]. I will to refer to this algorithm as the modified Welch-Berlekamp algorithm.
5.2.1 Algorithm Outline
The problem of decoding is to find the unique degree k − 1 or less polynomial m(x), such
that
m(αi) = ri

for all but t values of i ∈ W, where t < (d−f)/2
There exists a non-zero, degree t polynomial q such that
q(αi)m(αi) = q(αi)ri  ,  for all i ∈ W.
In particular, if E is the set of the (unknown) indices of the errors in ~r, and q(x) is the
monic error locator polynomial g(x) = PRODUCTSERIES( i∈E (x − αi) )     ,
then the above equation will hold for
q(x) = cg(x) for any constant c. 
Let s(x) = q(x)m(x). deg(s) ≤ k + t − 1, and
s(αi) = q(αi)ri,  for all i ∈ W.

Let 
s(x) = SERIES(from i=0 to k+t-1) (si * x^i)    (in GF(256) world)
and
q(x) = SERIES(from i=0 to t) (qi * x^i)        (in GF(256) world)
Then there are k + 2t + 1 unknowns that
define s(x) and q(x), and (5.1) yields n − f independent equations. Since 2t < d − f =
n − k − f + 1, the number of unknowns is at most k + n − k − f + 1 = n − f + 1 and there
is at most one more unknown than there are equations. By requiring that q is non-zero,
we are assured a solution that is unique up to a constant factor.

5.2.2 Computing m(x) and Recovering ~m
Once a solution s(x), q(x) is found, we can recover m(x) by m(x) = s(x)/q(x). 
If the encoding was non-systematic, then the coefficients of m(x) form the message vector ~m. 
If the encoding was systematic, then we must compute
~m = m( ~α:0 to k−1 ).
Notice that in previous sections
we have referred to the polynomial m(x) found here as µ(x) in the systematic encoding
case; in the algorithm outlined in Section 5.2.6 we will use this notation to emphasize the
assumption of systematic encoding.

5.2.3 t is Unknown
The technique just outlined would lead immediately to an algorithm for finding m(x)
if t were known a priori. However, we generally do not know how many errors are present
in ~r  before decoding, just that there are some errors. 
We will assume that t < (d−f)/2
and
therefore that maximum likelihood decoding will succeed. 
Let tˆ be the maximum allowable value of t, 
tˆ = floor( (d−f−1)/2 ) == floor( (n−k−f)/2 ).
We can rewrite the above technique in terms of the
assumption that t ≤ tˆ without assuming that we know t:
We seek the unique degree k-1 or less polynomial m(x), such that
m(αi) = ri
for all but at most tˆ values of i ∈ W. 
There exists a non-zero polynomial q such that
deg(q) ≤ tˆ, and
q(αi)m(αi) = q(αi)ri  ,
for all i ∈ W.
Let s(x) = q(x)m(x). Then

∀i ∈ W, s(αi) = q(αi)ri with deg(s) ≤ k + tˆ− 1 and deg q ≤ tˆ . (5.2)

Let s(x) = SERIES(from i=0 to (k + tˆ − 1)) (si * x^i)
and
q(x) = SERIES(from i=0 to tˆ) (qi * x^i)   .
Then there are k + 2tˆ + 1 = n−f+1 unknowns
that define s(x) and q(x), and (5.2) yields n − f equations which are not necessarily
independent. This system is underdetermined and thus many solutions exist. Nonetheless,
we shall see that any solution to (5.2) yields the same value for m(x), and that we may
use standard linear algebra techniques to find such a solution.

Theorem 5.2.1. Let s, q be a solution to (5.2). Then s(x)/q(x) = m(x) where m(x) is the
unique polynomial such that deg m ≤ k − 1 and m(αi) = ri for all but at most tˆ values of
i ∈ W.


*************************************************/


package Berlekamp_Welch_Algorithm.test;

import java.util.ArrayList;


public class Pages62_66_Modified_Welch_Berlekamp_Algorithm {
    
    private int k=0; //number of message-bytes
    private int[] r_vector = null;//received-bytes
    private int n = 0; //number of received-bytes systematic-encoding, same length as r_vector.length
    //inputs to x in polynomials, same length as r_vector.length
    private int[] a_vector = null;
    
    //F = set of byte-indices in received-bytes which you "erase" to do algorithm, can only erase in the E.C-bytes because it is systemic encoding
    private ArrayList<Integer> F = null;
    //fˆ == f == length of F
    private ArrayList<Integer> r_vector_ECbytes_erased = null;
    
    private int t_hat = 0; //maximum-errors that can be detected and corrected
    
    private int[][] tempA=null;
    private int[][] A=null;
    private int[] answer_matrix = null; //DO NOT USE answer_matrix[0]
    private int q_count = 0;
    private int[] q = null;
    private int s_count = 0;
    private int[] s = null;
    
    //Calculated after identity-matrix solved if solvable
    ArrayList<Integer> u = null;
    ArrayList<Integer> remainder = null;
    ArrayList<Integer> m = null;
    
    public Pages62_66_Modified_Welch_Berlekamp_Algorithm(Integer[] param_r_vector, int param_k)
    {
        k = param_k;
        
        r_vector = new int[param_r_vector.length];
        for (int ii=0; ii < param_r_vector.length; ii++)
            r_vector[ii] = param_r_vector[ii];
        
        n = r_vector.length; //number of received-bytes systematic-encoding, same length as r_vector.length
        
        a_vector = new int[n];
        for (int ii=0; ii <  n; ii++)
            a_vector[ii] = ii; //0, 1, 2, 3, 4, 5 .....
        
        F = new ArrayList<Integer>(); //F = set of byte-indices in received-bytes which you "erase" to do algorithm, can only erase in the E.C-bytes because it is systemic encoding
        
        r_vector_ECbytes_erased = new ArrayList<Integer>();
        for (int ii=0; ii < r_vector.length; ii++)
            r_vector_ECbytes_erased.add(r_vector[ii]);
        
    }//constructor
    
    
    public Pages62_66_Modified_Welch_Berlekamp_Algorithm(int[] param_r_vector, int param_k)
    {
        k = param_k;
        
        r_vector = param_r_vector;
        
        n = r_vector.length; //number of received-bytes systematic-encoding, same length as r_vector.length
        
        a_vector = new int[n];
        for (int ii=0; ii <  n; ii++)
            a_vector[ii] = ii; //0, 1, 2, 3, 4, 5 .....
        
        F = new ArrayList<Integer>(); //F = set of byte-indices in received-bytes which you "erase" to do algorithm, can only erase in the E.C-bytes because it is systemic encoding
        
        r_vector_ECbytes_erased = new ArrayList<Integer>();
        for (int ii=0; ii < r_vector.length; ii++)
            r_vector_ECbytes_erased.add(r_vector[ii]);
        
    }//constructor
    
    /**********************************************
    1. Lines: 4 - 10
    (F = set of byte-indices in received-bytes which you "erase" to do algorithm, in the E.C-bytes because it is systemic encoding)
    (f = length of F)
    (Fcomplement = set of byte-indices in received-bytes which you "do NOT erase" to do algo)
    Purpose: Let fˆ = f.                    
    If n−k−f is odd, add the last non-erased index in ~r to F
    and increment fˆ.
    fˆ is the number “psychological” number of erasures. 
    Let tˆ = (n−k−fˆ)/2
    Let nˆ = n − ˆf
    Let W = Fcomplement = set of byte-indices in received-bytes which are NOT "erased" to do algorithm
    Computational complexity: No Galois field computation required. 
      
    **********************************************/
    
    private void Lines_4_to_10()
    {
        /*************************************************
        4: if n − k − f odd then
        5: Add the last non-erased index in ~r to F
        6: ˆf ← ˆf + 1
        7: end if
        8: tˆ← (n−k−fˆ)/2
        9: nˆ ← n − ˆf
        10: W ← Fcomplement
        *************************************************/
        int temp = n - k - F.size();
        if ((temp % 2) == 1) //odd
        {
            int EC_byte_index_erased = this.r_vector_ECbytes_erased.size() - 1;
            F.add( EC_byte_index_erased );
            this.r_vector_ECbytes_erased.remove(EC_byte_index_erased);
        }
        
        temp = n - k - F.size();
        this.t_hat = temp / 2; //maximum-errors that can be detected and corrected
        
        answer_matrix = new int[this.r_vector_ECbytes_erased.size() + 1]; //DO NOT USE answer_matrix[0]
        q_count = this.t_hat + 1;
        q = new int[q_count];
        s_count = this.k + this.t_hat;
        s = new int[s_count];
        s[s_count - 1] = 1;
    }//Lines_4_to_10
    
    
    /****************************************************
    2. Lines: 12 - 14
    Purpose: Construct the matrix A from (5.3)
    
    A = matrix( diag(~rW) * Vandermonde(n−f,tˆ+1 (~αW)) − Vandermode(n−f,k+tˆ (~αW)) )   (5.3)
    
    Let 
    ~s = (s[0], . . . , s[k + tˆ − 1])T 
    and 
    ~q = (q[0], . . . , q[tˆ])T
    
    then when
    A
    can be trans formed to identity-matrix or identity-matrix with some all-0-rows at the bottom
    then you know
    [ q[0] ... q[t^] s[0] ... s[k + tˆ − 1]
    then you know
    polynomials q(x) and s(x)
    and then if
    u(x) = s(x) / q(x) and remainder is all-0 then m(x) == u(x) using inputs 0,1,2,3,4,5, ..... ii-1  remember in GF(256) world

    *****************************************************/
    public void Lines_12_to_14() throws Exception
    {
        //12: B ← Vand(nˆ,tˆ+ 1, ~αW)
        //13: C ← Vand(nˆ, k + tˆ, ~αW)
        //14: A ← [diag(~rW)B − C]
        
        Matrix mm = new Matrix();
        int[][] B = mm.Vandermonde(this.r_vector_ECbytes_erased.size(), this.t_hat + 1, this.a_vector);
        int[][] C = mm.Vandermonde(this.r_vector_ECbytes_erased.size(), this.k + this.t_hat, this.a_vector);
        //mm.Debug_Print(C);
        int[][] diag_r_vec_erased = mm.Diag(this.r_vector_ECbytes_erased);
        //mm.Debug_Print(diag_r_vec_erased);
        //mm.Debug_Print(B);
        int[][] matrix=null;
        
        try
        {
            matrix = mm.Matrix_Multiply(diag_r_vec_erased, B);
            //mm.Debug_Print(matrix);
            //mm.Debug_Print(C);
            tempA = mm.Matrix_Concatenate(matrix, C);
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
        
    }//Lines_12_to_14
    
    
    
    
    public boolean go()
    {
        try
        {
            Lines_4_to_10();
            Lines_12_to_14();
            Matrix mm = new Matrix();
            this.A = mm.MoveLastColumn(this.tempA, this.answer_matrix);
            boolean identitymatrix = mm.Robot_Solve(this.A, this.answer_matrix);
            System.out.println("identitymatrix = " + identitymatrix);
            if (identitymatrix)
            {
                this.Debug_Print();
                this.Calculate_u();
                boolean remainder_is_0_vector = this.Remainder_Is_0_Vector();
                if (remainder_is_0_vector)
                {
                    this.Calculate_m();
                    return true;
                }
            }
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
        
        return false;
    }//go
    
    
    public void Calculate_u() throws Exception
    {
        int aa=0; 
        for (int ii=0; ii < this.q_count; ii++)
        {
            q[ii] = answer_matrix[aa];
            aa++;
        }
        for (int ii=0; ii < this.s_count; ii++) //s[this.s_max_index] == 1
        {
            if (ii == (this.s_count - 1))
                s[ii] = 1;
            else
            {
                s[ii] = answer_matrix[aa];
                aa++;
            }
        }
        
        ArrayList<Integer> divisor = new ArrayList<Integer>();
        //s has a 1 at the end but we added it already
        //i thought s(x)/q(x) to get u(x)
        //would be simple polynomial-long-division
        //but i may be wrong
        for (int ii=0; ii < s_count; ii++)
            divisor.add(s[ii]);
        
        u = new ArrayList<Integer>();
        Tools tool = new Tools();
        
        while (divisor.size() >= this.q_count)
        {
            int mult=0;
            int mult_exp=0;
            int dd=0;
            for (int qq=0; qq < this.q_count; qq++)
            {
                if (qq==0)
                {
                    int num_to_kill = divisor.get(dd);
                    dd++;
                    int num_to_divide_by = q[qq];
                    int num_to_kill_as_alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[num_to_kill];
                    int num_to_divide_by_as_alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[num_to_divide_by];
                    mult_exp = num_to_kill_as_alpha_exp + 255 - num_to_divide_by_as_alpha_exp;
                    if (mult_exp >= 256)
                        mult_exp %= 255;
                    int check_exp = mult_exp + num_to_divide_by_as_alpha_exp;
                    if (check_exp >= 256)
                        check_exp %= 255;
                    if ((check_exp != num_to_kill_as_alpha_exp) && (check_exp!=255 && num_to_kill_as_alpha_exp!=0))
                        throw new Exception("Calculate_u: check_exp != num_to_kill_as_alpha_exp");
                    mult = tool.Table_Exponent_Of_Alpha_To_Integer()[mult_exp];
                    u.add(mult);
                }
                else
                {
                    int top = divisor.get(dd);
                    int bottom = q[qq];
                    int bottom_as_alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[bottom];
                    int total_exp = bottom_as_alpha_exp + mult_exp;
                    if (total_exp >= 256)
                        total_exp %= 255;
                    int num = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    int gfsubtractnum = top ^ num;
                    divisor.set(dd, gfsubtractnum);
                    dd++;
                }
            }//for
            
            divisor.remove(0);
        }//while
        
        
        remainder = new ArrayList<Integer>();
        for (int dd=0; dd < divisor.size(); dd++)
            remainder.add(divisor.get(dd));
        
        System.out.println();
        System.out.print("remainder = ");
        for (int rr=0; rr < remainder.size(); rr++)
        {
            System.out.print(remainder.get(rr) + " ");
        }
        
        System.out.println();
        System.out.print("u = ");
        for (int uu=0; uu < u.size(); uu++)
        {
            System.out.print(u.get(uu) + "x^" + uu);
            if (uu != (u.size()-1))
                System.out.print(" + ");
        }
        System.out.println();
    }//Calculate_u
    
    
    public boolean Remainder_Is_0_Vector()
    {
        for (int rr=0; rr < remainder.size(); rr++)
        {
            if (remainder.get(rr) != 0)
                return false;
        }
        return true;
    }//Remainder_Is_0_Vector
    
    public void Calculate_m()
    {
        Tools tool = new Tools();
        m = new ArrayList<Integer>();
        
        for (int ii=0; ii < this.a_vector.length; ii++)
        {
            //x == 0
            if (this.a_vector[ii] == 0)
                m.add(u.get(0)); //apparently u[0] is coefficient*x^0
            //x == 1 == a[2]
            else if (this.a_vector[ii] == 1)
            {
                int sum=0;
                for (int uu=0; uu < u.size(); uu++)
                    sum ^= u.get(uu);
                m.add(sum);
            }
            //x >= 2
            //x == a[ii]
            //only use some a[ii] because some vectors thrown away
            else if (ii < u.size())
            {
                int sum=0;
                int x = this.a_vector[ii];
                int x_alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[x];
                for (int u_exp=0; u_exp < u.size(); u_exp++)
                {
                    int coeff = u.get(u_exp);
                    int coeff_alpha_exp = tool.Table_Integer_To_Exponent_Of_Alpha()[coeff];
                    int total_exp = coeff_alpha_exp + (x_alpha_exp * u_exp);
                    if (total_exp >= 256)
                        total_exp %= 255;
                    int term = tool.Table_Exponent_Of_Alpha_To_Integer()[total_exp];
                    sum ^= term;
                }
                m.add(sum);
            }
            else
                break;
        }
        
        System.out.println();
        System.out.print("correct message = ");
        for (int mm=0; mm < m.size(); mm++)
            System.out.print(m.get(mm) + " ");
        System.out.println();
    }//Calculate_m
    
    public void Debug_Print()
    {
        System.out.println();
        for (int row=0; row < A.length; row++)
        {
            System.out.print(" |");
            for (int col=0; col < A[0].length; col++)
            {
                String tempstr = Integer.toString(A[row][col]);
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
    
}//class
