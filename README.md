Shamir Secret Reconstruction

Reconstruct the secret from shared points using Shamir-style polynomial interpolation. Given n shares (x, y) with mixed bases for y, and a threshold k = m + 1 (m is degree), the task is:
- Parse JSON shares.
- Convert base-encoded y values to integers.
- Reconstruct the polynomial’s constant term P(0) using Lagrange interpolation.
- Detect any wrong (tampered) shares by inlier checking.

## Problem summary
- Each share is a point (x, y) where x is the JSON object key (e.g., "1", "2", …) and y is the value parsed from its given base.  
- With any k correct shares, the unique degree-(k−1) polynomial P(x) can be reconstructed.  
- The output secret is the constant term: P(0). Wrong shares are those that don’t lie on the interpolated polynomial.  

## Approach
1. Parse JSON input: read keys.n, keys.k, then each share’s base and value.  
2. Base conversion: convert each value string from its base (2–16+) to a big integer without using built-in interpolation helpers.  
3. Exact interpolation: compute P(0) using Lagrange interpolation with exact rational arithmetic:  
   - P(0) = Σ y_i · Π_{j≠i} [(-x_j)/(x_i − x_j)]  
   - Also support evaluating P(x) at arbitrary x via Lagrange form for verification.  
4. Inlier detection: iterate over all k-sized subsets, interpolate, then evaluate at all n x’s; choose the subset with the most exact matches. Shares that don’t match are flagged wrong.  

This guarantees an exact integer secret and robustly identifies any manipulated shares.

## Build and run
- Language: C++17 (no Python as per constraints)  
- Dependencies: uses standard library and big-integer rational arithmetic (document what your code uses, e.g., Boost.Multiprecision for big ints and a simple Fraction class).

Build:
- g++ -std=gnu++17 -O2 -pipe -static -s main.cpp -o solve

Run:
- ./solve < input.json

Input format:
- Provide the JSON exactly as in the assignment prompt.

Output:
- secret=<integer_constant_term>  
- wrong_share_indices=[comma-separated x indices]

## Correctness notes
- Uses exact rationals (numerator/denominator reduced by gcd) to avoid floating-point error.  
- Lagrange basis is computed directly at 0 for P(0) and at each x for verification.  
- For n up to 10 and k up to 7, checking all combinations is fast.

## Time complexity
- Combinations: C(n, k).  
- Each interpolation/evaluation: O(k^2).  
- Feasible for the provided test sizes.

## Sample results
- Output for TestCase-1: secret = 3  
- Wrong Data Set Points for TestCase-1: []  
- Output for TestCase-2: secret = 79836264049851000  
- Wrong Data Set Points for TestCase-2:[2][8]

## Repository structure
- main.cpp — complete solution (parsing, base conversion, interpolation, inlier detection).  
- README.md — this file.  
- testcases/ — optional folder with the two JSON files from the prompt.

## How to verify
- Cross-check by selecting any k inlier shares, reconstruct P(0), then evaluate P(x) for all x; correct shares match exactly and wrong ones don’t.  
- For TestCase-2, indices 2 and 8 fail the verification against the best-fitting polynomial.

## Notes
- If a prime modulus were specified, this would be over a finite field; here it’s plain integers/rationals with exact arithmetic.  
- Ensure base parsing supports digits a–f/A–F for bases ≥ 11.

## Author
Naveenkumar J — CSE student.

