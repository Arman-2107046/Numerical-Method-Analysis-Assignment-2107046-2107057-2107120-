# Numerical-Method-Analysis-Assignment-2107046-2107057-2107120-

A program to solve equations using numerical methods

<h3>Author-Arman Rahman Rafi.</h3>
<h3>Roll-2107046</h3>
<h5>Methods of solving system of linear equations</h5>
<h4>Included methods--</h4>

1.Jacobi Iterative Method<br>

<p>The code applies the Jacobi Iterative Method to approximate solutions for a system of linear equations. First, the user inputs the size n and values for an augmented matrix.
Each row represents an equation, with the last number being the constant. The code begins with zero as an initial guess for each variable and then iteratively refines these guesses. Two settings control the process: a tolerance of 1×10^−6, defining how close the solution needs to be, and a limit of 1000 iterations.
In each iteration, a new value for each variable is calculated based on the latest guesses. If any change exceeds the tolerance, another iteration is needed. The process stops once all changes are within tolerance, and the code outputs the approximated values as the solution.</p>

2.Gauss Seidel Iterative Method<br>

<p>
 Gauss-Seidel Iterative Method tries to solve a system of linear equations by refining guesses for each variable step-by-step. After entering the size and values of an augmented matrix, including both the coefficients of each variable and the constants, the code starts with all variables set to zero. It then updates each variable by plugging in the most recent values from the same iteration, making each guess more accurate as it progresses. The process stops once changes between guesses are 1×10^−6 or after 1000 tries if it hasn’t converged. Finally, it outputs the closest possible values for each variable, giving an approximate solution to the system.</p>
 
3.Gauss Elimination Method<br>
<p>
The GaussElimination method helps solve a system of linear equations by transforming the equations into an upper triangular matrix and then solving for each variable. First, it takes the system size and the augmented matrix entries from the user. The function then applies a process called forward elimination, where each row is adjusted using pivot elements to zero out the values below the diagonal, creating a "triangular" shape in the matrix. To avoid errors, it checks if any pivot (diagonal element) is zero before proceeding, stopping if a division by zero would occur. Once the matrix is in upper triangular form, the function uses back substitution to determine each variable's value, starting from the last row and moving upward. Finally, it displays the solution. This function could be further improved by showing the matrix after each step or formatting the output to be clearer, especially for larger systems.</p>

4.Gauss Jordan Elimination method<br>

<p>
The GaussJordanElimination method solves a system of linear equations by transforming the augmented matrix into a reduced row echelon form, where each variable has a coefficient of 1, and all other entries in its column are zero. After taking the system size and matrix entries as input, the function begins the elimination phase, where it uses each row as a pivot to zero out all other entries in that column, both above and below the pivot. This process continues until the matrix has only 1s on the diagonal and 0s elsewhere in each column, making it easy to read off the solution directly from the last column. The function carefully avoids division by zero by checking each pivot element before dividing, reporting an error if a zero pivot is encountered. Once the matrix is fully simplified, the solution is displayed. This approach provides a clear, direct solution and could be further enhanced by displaying each matrix step for easier understanding, especially in larger systems.</p>

5.LU Factorization Method<br>

<p>
The LUFactorization method solves a system of linear equations by decomposing the matrix into a lower triangular matrix (L) and an upper triangular matrix (U) such that A=LU. This process starts by taking the matrix as input, then successively eliminating elements below the pivots to form U while storing the multipliers in L. This factorization allows solving the system in two steps: first solving Ly=b through forward substitution, then solving ,Ux=y through back substitution to find the solution vector x. This approach is efficient and avoids modifying the original matrix, making it useful for large systems or multiple right-hand sides.</p>
<br>

# Numerical-Method-Analysis-Assignment-2107046-2107057-2107120-

A program to solve equations using numerical methods

<h3>Author-Megha Tania</h3>
<h3>Roll-2107057</h3>
<h5>Methods of solving polynomial equations</h5>
<h4>Included methods--</h4>
1.Bisection Method<br>
<p> 
Step 1: Firstly, we calculated an interval of the equation f(x), in which, the roots will be present.
Step 2: Then we ran a loop to find two points a and b, for which f(a)*f(b)<0. If we find a point x in the interval for which f(x)=0 then x is a root.
Step 3: we calculated midpoint c=(a+b)/2. . If f(c)=0 or less than tolerance (which is close to 0, or 0.00001)then c is a root.
Step 4: if f(a)*f(c)<0 then set b=c, otherwise set a=c. and repeat step 3 and 4.
</p>

2.False Position Method<br>
<p> 
Step 1: Firstly, we calculated an interval of the equation f(x), in which, the roots will be present.
Step 2: Then we ran a loop to find two points a and b, for which f(a)*f(b)<0. If we find a point x in the interval for which f(x)=0 then x is a root.
Step 3: we calculated c = (a * f(b) - b * f(a)) / (f(b) – f(a)).
 If f(c) or less than tolerance (which is close to 0, or 0.00001)then c is a root.
If prev_c-c is closer to 0 then c is a root.
Step 4: if f(a)*f(c)<0 then set b=c, else set a=c.
Set prev_c=c; repeat step 3 and 4.
</p>

3.Newton Raphson Method<br>
<p>
Step 1: initially guess a number x. Then calculate nextx=x-f(x)/f'(x).
Step 2: if (nextx-x) is equal to 0, or 0.0001, then nextx is the root.
Step 3: set x=nextx.
Repeat step 2 and 3.
Here, to get all the roots of the equation, we used "deflate" function to store all the roots, and then used it to print all the roots. 
</p>
4.Secant method<br>
<p>
step 1:Start with two initial guesses, x1 and x2,close to the expected root.
step 2: Calculate the next approximation using,
x3 = (x1 * fx(x2, arr) - x2 * fx(x1, arr)) / (fx(x2, arr) - fx(x1, arr));

step 3: If (abs(x1 - x2) < tolerance), then x3 is the root.
step 4:  Set x1=x2 and x2=x3 and repeat step 2 and step 3.

​
  using: </p>

# Numerical-Method-Analysis-Assignment-2107046-2107057-2107120-

A program to solve equations using numerical methods

<h3>Author-Siyam Khan</h3>
<h3>Roll-2107120</h3>
<h5>Runge Kutta method and Newton forward-backward Interpolation method</h5>
<h4>Included methods--</h4>
1.Runge Kutta Method<br>
<p>The Runge-Kutta method is a straightforward and effective way to solve ordinary differential equations when you can’t find exact solutions. It works by taking several measurements of how the function is changing at different points during each step, then combining those measurements to get a better estimate of the function’s next value. The fourth-order version, often called RK4, is especially popular because it strikes a good balance between accuracy and ease of use. It’s commonly used in various fields like physics, engineering, and finance to model complex systems.</p>

2.Newton Forward Backward Interpolation Method<br>

<p>
Newton's Forward and Backward Interpolation Methods are handy tools for estimating the values of a function when you only have a few known data points. Forward Interpolation is great for when you want to find a value within the range of your data, using a polynomial based on the starting points. On the other hand, Backward Interpolation works well when you're looking for a value near the end of your data set, relying on the last points to create its polynomial. Both methods use divided differences to build these polynomials, making them effective for filling in gaps in your data and getting a better sense of the function's behavior.</p>
