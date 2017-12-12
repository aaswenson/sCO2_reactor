MCNP6 base input file for leakage study. 
c Cell Card
1 1 -6.25 -1 -2 3 imp:n=1  $ cylindrical reactor
2 0 2:-3:(1 -2 3) imp:n=0  $ outside world
 
c Surface Card
1 CZ 50
2 PZ 25
3 PZ -25

c Data Card
m1
     1001 -7.4958e-03
     1002 -1.7229e-06
     6012 -1.6247e-05
     6013 -1.9042e-07
     8016 -1.6486e-01
     8017 -6.6741e-05
     14028 -1.8876e-04
     14029 -9.9318e-06
     14030 -6.7804e-06
     15031 -9.4518e-06
     16032 -5.8384e-06
     16033 -4.7539e-08
     16034 -2.7753e-07
     16036 -6.9144e-10
     24050 -3.2588e-04
     24052 -6.5352e-03
     24053 -7.5531e-04
     24054 -1.9156e-04
     25055 -4.1095e-04
     26054 -1.6280e-03
     26056 -2.6502e-02
     26057 -6.2299e-04
     26058 -8.4362e-05
     28058 -2.5544e-03
     28060 -1.0178e-03
     28061 -4.4982e-05
     28062 -1.4578e-04
     28064 -3.8318e-05
     92235 -4.7189e-02
     92238 -7.3929e-01
kcode 5000 1 40 25
ksrc 0 0 0
     1 1 1
     -1 -1 -1 
     1 -1 1
     -1 1 1
     1 1 -1
mode n
print
