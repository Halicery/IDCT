# Scaled LLM[#LLM 1989 (C. Loeffler, A. Ligtenberg, and G. S. Moschytz, Practical fast 1-D DCT algorithms with 11 multiplications] 1-D IDCT with 6 multiplications 

&copy; A. Tarpai   

Source code for this scaled integer IDCT can be found on [BitBucket"https://bitbucket.org/Halicery/" target="_blank"], as part of an ITU T.81 JPEG decoder I have written. 

## Results

Using any of these scaling vectors the number of multiplications can be reduced from 11 to 6 in a modified 1-D transform based on the LLM method:

**v** = (1, &eta;, &beta;, &gamma;&eta;, 1, &gamma;&eta;, &alpha;, &eta;)

**v** = (1, &delta;, &beta;, &gamma;&delta;, 1, &gamma;&delta;, &alpha;, &delta;) 

**v** = (1, &theta;, &beta;, &gamma;&theta;, 1, &gamma;&theta;, &alpha;, &theta;) 

**v** = (1, &epsilon;, &beta;, &gamma;&epsilon;, 1, &gamma;&epsilon;, &alpha;, &epsilon;) 

The motivation for these algorithms was to implement an acceptable and fast 8x8 integer IDCT for my JPEG and video decoders, without borrowing code. The focus is therefore on the inverse transform. 


## Principles of scaled IDCT

The inverse DCT is the very last stage of the image decoding process, which reconstructs the original image sample values. After coefficient decoding each value of the block is multiplied by the quantization table values (*de-quantization*). The transform itself is quite computation intensive, the 2-D formula is an iteration of 8x8x8x8 = 4096 multiplications with irrational numbers, based on the equation 

["formula_2d.png" width="300em"]: 

	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .     x Q-block     .  .  .  .  .  .  .  .       IDCT        .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .    --------->     .  .  .  .  .  .  .  .    --------->     .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	    coeff block                            de-quantized block                          image samples

By clever *pre-multiplication* (scaling) of the 8x8 input block before the IDCT, some calculations can be eliminated, others gets simpler by using a modified IDCT equation. Scaling is only efficient when it's performed together with the de-quantization step, that is the Q-block is already multiplied by the scaling block before the transform, the scaling matrix *SQ = S x Q*: 

	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .      simple       .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .    x SQ-block     .  .  .  .  .  .  .  .       IDCT        .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .    --------->     .  .  .  .  .  .  .  .    --------->     .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	    coeff block                            de-quantized                                image samples
	                                              and scaled block

## 1-D transform

The 2-D IDCT above is tremendously slow to implement. Most of the fast DCT algorithms are based on the separability of the 2-D transform into successive 1-D row- and column transforms (or column- and row transform). The 8-point 1-D transform computes 8 outputs from 8 inputs based on ["formula_1d.png" width="200em"] and can be implemented as matrix multiplication. An iteration of 8x8 = 64 multiplications with irrational numbers. The transform is linear (this will be important later): 

		                                   1-D transform
		X¨0¨  X¨1¨  X¨2¨  X¨3¨  X¨4¨  X¨5¨  X¨6¨  X¨7¨     --------------->    x¨0¨  x¨1¨  x¨2¨  x¨3¨  x¨4¨  x¨5¨  x¨6¨  x¨7¨

This transform function is applied 16 times for the 8x8 block as 8 row- and 8 column-transforms (16 x 64 = 1024 multiplication). Usually implementations use only one 1-D routine working on input and output arrays. In order to use the same function, the result is written transposed through a temp-block:


	0  1  2  3  4  5  6  7                   0  .  .  .  .  .  .  .                   0  1  2  3  4  5  6  7 
	.  .  .  .  .  .  .  .                   1  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .     1D IDCT       2  .  .  .  .  .  .  .     1D IDCT       .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .    --------->     3  .  .  .  .  .  .  .    --------->     .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .       row         4  .  .  .  .  .  .  .     column        .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   5  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   6  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .
	.  .  .  .  .  .  .  .                   7  .  .  .  .  .  .  .                   .  .  .  .  .  .  .  .


## Fast 1-D DCT: the LLM method

Among the endless number of fast algorithms I've been fascinated by Löffler's paper from 1989, its clean and simple solution. It gives a method the compute the 8-point 1-D transform with 11 multiplications:

["LLMIDCT.GIF" width="50%"] 
§Image source: Reznik, Hindsy, Zhangz, Yuz, and Ni, Efficient Fixed-Point Approximations of the 8x8 Inverse Discrete Cosine Transform§

As you can see there are no multiplications needed for X¨0¨ and X¨4¨. This is because LLM *uniformly scales* the orthonormal DCT matrix by &radic;8, i.e. the DCT matrix's vectors are of lenght &radic;8 (instead of 1 as for an orthonormal matrix). This trick makes v¨0¨ and v¨4¨ all having coordinates of 1. Applying this DCT matrix twice (row/column method) gives results multiplied by 8 = (&radic;8)^2^ - which is easy to implement as right shift to get the correct values after the transform. 

LLM uses these 7 irrational constants:

&gamma; = &radic;2 &ap; 1.414213562

&alpha; = &radic;2cos(3&pi;/8) &ap; 0.5411961
&beta; = &radic;2sin(3&pi;/8) &ap; 1.306562965

&eta; = cos(3&pi;/16) &ap; 0.831469612
&theta; = sin(3&pi;/16) &ap; 0.555570233

&delta; = cos(&pi;/16) &ap; 0.98078528
&epsilon; = sin(&pi;/16) &ap; 0.195090322


After rearranging the order of inputs 2 basic structures reveal: 

- an adder (blue)
- and a butterfly-multiplier (red):  


["LLMIDCT_TA.gif" width="50%"]

The adder takes 4 inputs (a, b, c and d) and computes 4 new outputs according to: 

	p  = a + b
	n  = a - b
	a' = p + d
	b' = n + c
	c' = n - c
	d' = p - d

The multiplier is the butterfly or rotation. It takes 2 inputs and 2 irrational constants. The constants are cos/sin pairs of &alpha;/&beta;, &eta;/&theta; or &delta;/&epsilon;; denoted here as **K/S**. Two outputs are computed according to the equation: 

	a' = Ka + Sb
	b' = Kb - Sa

The equation needs 4 multiplications to evaluate. As in the original paper, it can be reduced to 3 using intermediates. There are 6 possibilities:

	c  = K(a+b)                 c  = S(a+b)       
	a' = c - b(K-S)             a' = a(K-S) + c 
	b' = c - a(K+S)             b' = b(K+S) - c 
	                                            
	c  = K(a-b)                 c  = S(a-b)       
	a' = b(K+S) + c             a' = a(K+S) - c 
	b' = a(K-S) - c             b' = b(K-S) - c 
	                                            
	c  = K(b-a)                 c  = S(b-a)       
	a' = b(K+S) - c             a' = a(K+S) + c 
	b' = a(K-S) + c             b' = b(K-S) + c 

These are mathematically equivalent but chosing one of them have an impact when designing the scaling matrix and the integer IDCT implementation. The first column is *K-type*, where the 3 new constants are K, K+S and K-S, while the second column the *S-type*, there the 3 new constants are S, K+S and K-S. The adder and the multiplier can be implemented as macros (`XADD`, `KROT/SROT`). 


## How the scaling matrix works

First lets see how to eliminate the &radic;2 multiplication of F¨3¨ and F¨5¨ inside a simplified LLM 1D-transform. 


		    +-------------+                                     +----------+         
		F¨0¨ -|---->       -|-->  f¨0¨                       F¨0¨ ----|->       -|-->  f¨0¨    
		F¨1¨ -|---->       -|-->  f¨1¨                       F¨1¨ ----|->       -|-->  f¨1¨    
		F¨2¨ -|---->       -|-->  f¨2¨                       F¨2¨ ----|->       -|-->  f¨2¨    
		F¨3¨ -|-&radic;2->       -|-->  f¨3¨   --->              &radic;2F¨3¨ ----|->       -|-->  f¨3¨
		F¨4¨ -|---->       -|-->  f¨4¨                       F¨4¨ ----|->       -|-->  f¨4¨    
		F¨5¨ -|-&radic;2->       -|-->  f¨5¨                     &radic;2F¨5¨ ----|->       -|-->  f¨5¨
		F¨6¨ -|---->       -|-->  f¨6¨                       F¨6¨ ----|->       -|-->  f¨6¨    
		F¨7¨ -|---->       -|-->  f¨7¨                       F¨7¨ ----|->       -|-->  f¨7¨    
		    +-------------+                                     +----------+         
		        original                                          modified


The solution uses scaled inputs: that is F¨3¨ and F¨5¨ entering the 1-D transform is already multiplied by &radic;2. Because if the 2-D transform is successive row- the column transforms with transposed output (using the same function for both), understanding the scaling matrix is a little complicated. So lets go backwards:


1. The modified 1-D transform of the second block works correctly when F¨3¨ and F¨5¨ of each row is already &radic;2-multiplied: 

	                                                                 0  1  2  3  4  5  6  7
	                                                                 
	                                                                 .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                                                                 .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                                                                 .  .  .  +  .  +  .  .   8 x 1D  .  .  .  .  .  .  .  .
	                                                         ----->  .  .  .  +  .  +  .  .   ----->  .  .  .  .  .  .  .  .
	                                                                 .  .  .  +  .  +  .  .    tr.p   .  .  .  .  .  .  .  .
	                                                                 .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                                                                 .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                                                                 .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                                                                                                           output

2. These two columns are the result of the modified 1-D transform of row 3 and 5 of the previous block (transposed output). Because the transform is linear, obtaining all outputs scaled by &radic;2 can be achieved by scaling all inputs by &radic;2 before the transform, using &radic;2f(**x**) = f(&radic;2**x**), where **x** is the (row) vector:  

	                                                                 0  1  2  3  4  5  6  7
	                                                                   
	                           0    .  .  .  .  .  .  .  .           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                           1    .  .  .  .  .  .  .  .           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                           2    .  .  .  .  .  .  .  .   8 x 1D  .  .  .  +  .  +  .  .   8 x 1D  .  .  .  .  .  .  .  .
	                           3    +  +  +  +  +  +  +  +   ----->  .  .  .  +  .  +  .  .   ----->  .  .  .  .  .  .  .  .
	                           4    .  .  .  .  .  .  .  .           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                           5    +  +  +  +  +  +  +  +           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                           6    .  .  .  .  .  .  .  .           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                           7    .  .  .  .  .  .  .  .           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                                                                                                           output

3. But the first block transform using the modified 1-D transform also expects F¨3¨ and F¨5¨ of each row multiplied by &radic;2. In row 3 and 5, where all inputs are already scaled by &radic;2, F¨3¨ and F¨5¨ will be *double scaled*, ie. (&radic;2)&sup2; = 2, the intersection points marked by **#**:  

	                                                                 0  1  2  3  4  5  6  7
	                                                                   
	                           0    .  .  .  +  .  +  .  .           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                           1    .  .  .  +  .  +  .  .           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                           2    .  .  .  +  .  +  .  .   8 x 1D  .  .  .  +  .  +  .  .   8 x 1D  .  .  .  .  .  .  .  .
	                           3    +  +  +  #  +  #  +  +   ----->  .  .  .  +  .  +  .  .   ----->  .  .  .  .  .  .  .  .
	                           4    .  .  .  +  .  +  .  .           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                           5    +  +  +  #  +  #  +  +           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                           6    .  .  .  +  .  +  .  .           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                           7    .  .  .  +  .  +  .  .           .  .  .  +  .  +  .  .           .  .  .  .  .  .  .  .
	                                                                                                           output

This gives our scaling matrix, **S**, to multiply values of the 8x8 coefficient block before applying the modified IDCT twice, first on rows then on colums, thus eliminating the &radic;2 multiplication of F¨3¨ and F¨5¨ inside the modified LLM 1-D transform: 


	          1   1   1  &radic;2   1  &radic;2   1   1
	          1   1   1  &radic;2   1  &radic;2   1   1
	          1   1   1  &radic;2   1  &radic;2   1   1
	   S =   &radic;2  &radic;2  &radic;2   2  &radic;2   2  &radic;2  &radic;2
	          1   1   1  &radic;2   1  &radic;2   1   1
	         &radic;2  &radic;2  &radic;2   2  &radic;2   2  &radic;2  &radic;2
	          1   1   1  &radic;2   1  &radic;2   1   1
	          1   1   1  &radic;2   1  &radic;2   1   1

**S** is also the result of a matrix multiplication of the column vector **v** and its transposed vector **v**^T^: 

**S** = **v v**^T^, 

where **v** = (1, 1, 1, &radic;2, 1, &radic;2, 1, 1): 

		      |    1   1   1  &radic;2   1  &radic;2   1   1     (**v**^T^)
		   ---|---------------------------------------------          
		      |                                          
		   1  |    1   1   1  &radic;2   1  &radic;2   1   1                 
		   1  |    1   1   1  &radic;2   1  &radic;2   1   1                 
		   1  |    1   1   1  &radic;2   1  &radic;2   1   1               
		  &radic;2  |   &radic;2  &radic;2  &radic;2   2  &radic;2   2  &radic;2  &radic;2       =  S
		   1  |    1   1   1  &radic;2   1  &radic;2   1   1                 
		  &radic;2  |   &radic;2  &radic;2  &radic;2   2  &radic;2   2  &radic;2  &radic;2                 
		   1  |    1   1   1  &radic;2   1  &radic;2   1   1                 
		   1  |    1   1   1  &radic;2   1  &radic;2   1   1                 
		  (**v**) |                                          

This essentially means that we have to focus only on finding a *good* scaling vector for a modified LLM transform, then we can simply obtain the final scaling matrix by computing **v v^T^** for the 2-D row/column method. Integrating this scaling matrix (*S*) into the quantization matrix (*Q*) by multiplying each value pairs yields the (real) *SQ* matrix - and a faster 1-D transform. 


## Scaling the even part butterfly

["even.gif" width="400em"]

For X¨2¨ and X¨6¨ the even part starts with a butterfly multiplication. The 2 irrational constants (K and S) are:

	K = &radic;2cos(3&pi;/8)
	S = &radic;2sin(3&pi;/8)

The outputs are computed according to the equation: 

		x¨6¨ = K X¨6¨ + S X¨2¨ 
		x¨2¨ = K X¨2¨ - S X¨6¨ 

We can obtain x¨6¨ without multiplication inside a modified LLM transform when the inputs are scaled: let W¨6¨ = K X¨6¨ and W¨2¨ = S X¨2¨. Then x¨6¨ is simply 

		x¨6¨ = W¨6¨ + W¨2¨

while for x¨2¨ we get:

		x¨2¨ = K/S W¨2¨ - S/K W¨6¨

The constants of K/S and S/K are trigonometric identities: 

	K/S = &radic;2cos(3&pi;/8) / &radic;2sin(3&pi;/8) = cot(3&pi;/8) = &radic;2 - 1
	S/K = &radic;2sin(3&pi;/8) / &radic;2cos(3&pi;/8) = tan(3&pi;/8) = &radic;2 + 1

then

		x¨2¨ = (&radic;2 - 1) W¨2¨  - (&radic;2 + 1) W¨6¨  = &radic;2 W¨2¨ - W¨2¨ - &radic;2 W¨6¨ - W¨6¨  = &radic;2 ( W¨2¨ - W¨6¨ ) - W¨2¨ - W¨6¨

**That is, the scaled even part needs only one multiplication of &radic;2.** 

§I discovered this accidently by looking at values in the Excel sheet I use to calculate factors. Indeed sin(3&pi;/8) = &half;&radic;(2+&radic;2) and cos(3&pi;/8) = &half;&radic;(2-&radic;2). Then K/S = &radic;( (2-&radic;2)/(2+&radic;2) ) = &radic;( (2-&radic;2)&sup2;/(4-2) ) = &radic; ( (4-4&radic;2+2) / 2 ) = &radic; ( 2-2&radic;2+1 ) = &radic; ( (&radic;2-1)&sup2; ) = &radic;2-1. S/K is similar.§ 

Now we can include &alpha; and &beta; in the scaling vector **v**, and use simplified equations for computing the rotation for the even part:

**v** = (1, 1, &beta;, &gamma;, 1, &gamma;, &alpha;, 1) = ( 1, 1, &radic;2sin(3&pi;/8), &radic;2, 1, &radic;2, &radic;2cos(3&pi;/8), 1 )



## Scaling the odd part butterfly

The odd-part starts with an adder: that is all 4 inputs must be equally scaled. That's a problem for simplifying the two butterfly multiplications for &eta;/&theta; and &delta;/&epsilon; that follows. The task is to find a common factor, *r*, which may simplify these two butterfly multiplications. When all 4 inputs of X¨1¨, X¨7¨, X¨5¨ and X¨3¨ are pre-scaled by *r*, then the equations become: 

x¨1¨ = &eta;/r X¨1¨ + &theta;/r X¨7¨
x¨7¨ = &eta;/r X¨7¨ - &theta;/r X¨1¨

x¨5¨ = &delta;/r X¨5¨ + &epsilon;/r X¨3¨
x¨3¨ = &delta;/r X¨3¨ - &epsilon;/r X¨5¨

["odd.gif" width="400em"]

The only thing I could do is to chose one of the constants for pre-scaling - thus reducing the number of multiplications in one of the butterfly to 2. F. ex. when r = &eta; we get: 

x¨1¨ = X¨1¨ + &theta;/&eta; X¨7¨
x¨7¨ = X¨7¨ - &theta;/&eta; X¨1¨

x¨5¨ = &delta;/&eta; X¨5¨ + &epsilon;/&eta; X¨3¨
x¨3¨ = &delta;/&eta; X¨3¨ - &epsilon;/&eta; X¨5¨

Furthermore, using one of the 3-mul equations for x¨5¨/x¨3¨, the total number of multiplications will be 5. F. ex. using the S(b-a) intermediate:

c = &epsilon;/&eta; (x¨3¨ - x¨5¨) 
x¨5¨ = c + ( &delta;/&eta; + &epsilon;/&eta; ) X¨5¨
x¨3¨ = c + ( &delta;/&eta; - &epsilon;/&eta; ) X¨3¨

Including r = &eta; in the scaling vector **v** for the odd part leads one of the the final solutions: 

**v** = (1, &eta;, &beta;, &gamma;&eta;, 1, &gamma;&eta;, &alpha;, &eta;)

and a modified 1-D LLM transform with 6 multiplications. 


### Detour

The last equations reveal something promising: &delta; + &epsilon; and &delta; - &epsilon; are the sum and difference of cos(x) and sin(x) of the same angle and are trigonometric identities! I was hoping these equations might help the IDCT calculation, somehow, especially in the scaled version. Not really. The complexity remains the same, but it is an interesting investigation:  

["pi16.png" width="400em"]

The odd part *runs* on the &pi;/16 line and there seem to be almost no relationship between sin/cos(&pi;/16) and sin/cos(3&pi;/16). Almost. 

I read this somewhere: 

sin(x) + cos(x) = &radic;2 sin( x + &pi;/4) 

I always wonder about trigonometric identities, how come is this true? It's pretty hard to see this with geometry, drawing triangles, so lets plot these 3 functions using [Fooplot"http://fooplot.com/" target="blank"] for example: 

**Plotting sin(x) + cos(x):**

["sinx+cosx.png"]

The fig above shows sin(x), cos(x) and sin(x)+cos(x). Ok, periodic, same as sin/cos on 2&pi;. Shifted on the x axis by +/- &pi;/4 compared to sin and cos. The maximum is at &pi;/4, where both sin(&pi;/4) = cos(&pi;/4) = &radic;2/2, ergo max = &radic;2. 

So sin(x) + cos(x) is something like:

sin(x) + cos(x) = &radic;2 sin(x + &pi;/4) or

sin(x) + cos(x) = &radic;2 cos(x - &pi;/4) or

sin(x) + cos(x) = &radic;2 sin(x - 7&pi;/4) or

sin(x) + cos(x) = &radic;2 cos(x + 7&pi;/4) or

sin(x) + cos(x) = &radic;2 cos(-x + &pi;/4).. and so on. 



Why is this interesting? Because in the LLM transform both &delta;/&epsilon; (&pi;/16) and &eta;/&theta; (3&pi;/16) are multiplies of &pi;/16 angles. Adding &pi;/4 to &pi;/16 is 5&pi;/16, which is *symmetric* to 3&pi;/16, sin(3&pi;/16)=cos(5&pi;/16) and cos(3&pi;/16)=sin(5&pi;/16). Now using the above equations we can state the following identity for the sum of &delta; + &epsilon;: 

**&delta; + &epsilon;** = cos(&pi;/16) + sin(&pi;/16) = &radic;2 sin(&pi;/16 + &pi;/4) = &radic;2 sin(5&pi;/16) = &radic;2 cos(3&pi;/16) = **&radic;2 &eta;**

and similarly for the sum of &eta; + &theta;:

**&eta; + &theta;** = cos(3&pi;/16) + sin(3&pi;/16) = &radic;2 sin(3&pi;/16 + &pi;/4) = &radic;2 sin(7&pi;/16) = &radic;2 cos(&pi;/16) = **&radic;2 &delta;**


**Plotting cos(x) - sin(x):**

For the other multiplier in the 3-mul butterfly version K-S appears, lets make identities for these two in regard of the constants of the LLM transform. 

["cosx-sinx.png"]

cos(x) - sin(x) = &radic;2 sin(x + 3&pi;/4) 

cos(x) - sin(x) = &radic;2 cos(x + &pi;/4) 

cos(x) - sin(x) = &radic;2 sin(x - 5&pi;/4) 

cos(x) - sin(x) = &radic;2 cos(x - 7&pi;/4) 


Now using the above we can state the following: 

**&delta; - &epsilon;** = cos(&pi;/16) - sin(&pi;/16) = &radic;2 cos(&pi;/16 + &pi;/4) = &radic;2 cos(5&pi;/16) = &radic;2 sin(3&pi;/16) = **&radic;2 &theta;**

and 

**&eta; - &theta;** = cos(3&pi;/16) - sin(3&pi;/16) = &radic;2 cos(3&pi;/16 + &pi;/4) = &radic;2 cos(7&pi;/16) = &radic;2 sin(&pi;/16) = **&radic;2 &epsilon;**

Giving som useful identities for the LLM transform: 

** &delta; + &epsilon; = &radic;2 &eta; **
** &delta; - &epsilon; = &radic;2 &theta; **
** &epsilon; - &delta; = -&radic;2 &theta;  **

and 

** &eta; + &theta; = &radic;2 &delta; **
** &eta; - &theta; = &radic;2 &epsilon; **
** &theta; - &eta; = -&radic;2 &epsilon; **

Again, I just could't find a better solution even with these identities for the scaled LLM transform. 


## Summary

Using scaled coefficients before the IDCT can be a significant improvement in speed. 

The idea is that scaling happens once for the whole 8x8 block together with dequantization (both are multiplications), simplifying computations for the *modified* 1D transform function. This is called a scaled-quantization table. 

For the integer version, only an extra post-scale stage is required in the form of right shift: 


	Scaled dequantization --> modified IDCT --> right shift


#### Using integer CPU arithmetics

When using integer arithmetics, a good approximation of **v** is used, where **V** = INT (2^n^ **v**). After the second 1D-transform, all values are shifted right by *n*: 

- pre-scale (multiplication)
- row-transform
- column-transform
- post-scale (right shift)

And an example integer matrix for scaled x3 and x5 by &radic;2. Here n=7 and &radic;2 &times; 128 = 181,019336 seems like a good integer approximation:

	128	128	128	181	128	181	128	128
	128	128	128	181	128	181	128	128
	128	128	128	181	128	181	128	128
	181	181	181	256	181	256	181	181
	128	128	128	181	128	181	128	128
	181	181	181	256	181	256	181	181
	128	128	128	181	128	181	128	128
	128	128	128	181	128	181	128	128

As for the modified LLM transform the final post-scale will be 7 + 3 = 10. 

