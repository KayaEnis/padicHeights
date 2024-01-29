This is SageMath (and some magma) code for the examples in "Algorithms for p-adic Heights on Hyperelliptic Curves of Arbitrary Reduction".

**Dependencies**

* Jennifer Balakrishnan's code for computing Coleman integrals of differentials of the third kind is required. Get it from https://github.com/jbalakrishnan/AWS and follow the instructions there.
* The third-named author's and Stevan Gajović's code for computing Coleman integrals of differentials of the third kind is required. Get it from https://github.com/StevanGajovic/heights_above_p and follow the instructions there.
* The first-named author's code for computing p-adic Néron functions is required. Get it from https://github.com/bianchifrancesca/padic_heights_g2 and move the folder p-adic_heights_g2 to the present folder. Also change the loading paths as follows:
  - In p-adic_heights_g2/heights_and_abelian_integrals.sage change
    
      ```load("./formal_group_expansions.sage")``` 

to

  ```load("../padic_heights_g2/formal_group_expansions.sage")```
  
and

  ```load("./division_polynomials.sage")```
  
to

  ```load("../padic_heights_g2/division_polynomials.sage")```
  
  - In p-adic_heights_g2/division_polynomials.sage change 

```load("./phi345.sage")```

to  

```load("../padic_heights_g2/phi345.sage")```

**Authors**
Francesca Bianchi, Enis Kaya and J. Steffen Müller

**Files**

* Ex5.1.py, Ex5.3.py and Ex5.4.py contain the Coleman-Gross height computations for, respectively, Examples 5.1, 5.3 and 5.4 from the paper. Ex5.4(fast).py is the same as Ex5.4.py, except that we use the code of the third-named author and Gajović instead of Balakrishnan's code for computing Coleman integrals of differentials of the third kind; it's faster than Ex5.4.py.
* Blakestad's canonical space can be computed using canonical_subspace.py.
* The quadratic Chabauty Example 5.5 can be computed by running qchabauty_example.sage from the subfolder QCExample in SageMath. The final part of the computation requires Magma.
* The local heights away from p in Examples 5.3 and 5.4 can be computed using local Neron functions by running example_5_3.sage and example_5_4.sage in the folder NeronExamples.

