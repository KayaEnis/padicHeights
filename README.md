This SageMath code for the examples in XXXXXXXX.

**Dependencies**

* Jennifer Balakrishnan's code for computing Coleman integrals of differentials of the third kind is required. Get it from https://github.com/jbalakrishnan/AWS and follow the instructions there.
* The third-named author's and Stevan Gajović's code for computing Coleman integrals of differentials of the third kind is required. Get it from https://github.com/StevanGajovic/heights_above_p and follow the instructions there.
* The first-named author's code for XXXXXXXXXXXXXXXXXXXXXX is required. Get it from https://github.com/bianchifrancesca/padic_heights_g2 and move the folder p-adic_heights_g2 to the present folder. Also change the loading paths as follows:
  - In p-adic_heights_g2/heights_and_abelian_integrals.sage change 
  load("./formal_group_expansions.sage")
  load("./division_polynomials.sage")
  to 
  load("../padic_heights_g2/formal_group_expansions.sage")
  load("../padic_heights_g2/division_polynomials.sage")
  - In p-adic_heights_g2/division_polynomials.sage change 
load("./phi345.sage") 
to  
load("../padic_heights_g2/phi345.sage")



**Authors**
Francesca Bianchi, Enis Kaya and J. Steffen Müller


**Files**

* The quadratic Chabauty example 5.5 can be computed by running qchabauty_example.sage from the subfolder QCExample in Sage. The final part of the computation requires magma.

