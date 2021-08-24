#!/bin/bash
cd test/cv_file/umb_2.2/
  ../../../bin/Probability_analysis.x ,\
 -T0 300                 ,\
 -T 1000                 ,\
 -bias_fact 1500         ,\
 -tmin 5000              ,\
 -ncv 5                  ,\
 -UCV 1                  ,\
 -MTD n                  ,\
 -MCV 0                  ,\
 -tool probT		 ,\
 -Prob_nD 1              ,\
 -CV_num 1 2               ,\
 -grid 1.0 5.0 0.02 1.0 10.0 0.02 1.0 9.0 0.02 3.0 5.0 0.02 1.0 6.0 0.05 ,\
 -pfrqMD 1 ,\
 -dtMTD 200
