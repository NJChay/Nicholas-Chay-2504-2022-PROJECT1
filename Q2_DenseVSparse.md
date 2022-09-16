Dense

Dense polynomial representations store term values for every degree up to the degree of the leading term, with 0 specified as the value for unspecified terms.
This clearly mean that dense polynomials take up more space in memory as a polynomial such as x^500 +1 would have 501 stored terms as opposed to 2 
This allows addition operations between dense polynomials to be done in linear time as the coefficient of each degree is in the same place for each degree in every dense polynomial.
For multiplication however, dense representation is less effective as the size of the output grows in a squared relation with the inputs. Thus as dense polynomials store significantly more inputs the output size grows much faster.
Additionally, functions such as modulo that require on iterating through each stored term are slower, again for the example above requiring an additional 500 iterations 
Therefore, more complex operations that rely on multiplication and modulo like GCD and factorization will become much slower for large polynomials with dense representation compared to sparse

Sparse

Sparse polynomials only store the terms that are specified when constructing the polynomial. 
Therefore, storage of sparse polynomials takes up less space in memory as only the relevant terms are stored.
However, operations like addition are less efficient as each term will have a unique location for each polynomial unlike dense where degree 0 is stored in the first index and degree one in the second and so forth. Therefore requiring additional lookup to find the coefficients of matching degrees when performing these operations.
As mentioned above, modulo and multiplication issues will be significantly faster as no lookup is required during polynomial multiplication and sparse polynomials simply have less terms to iterate through for modulo
