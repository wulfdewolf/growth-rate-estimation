# Growth Rate Estimation
The matlab file in this repository contains some functions to estimate growth-rate from the `.txt` file outputted by the OGI-BIO Bioreactor.

# 1. Usage
The main function `compute_growth_rate` takes as input:
* the path to the `.txt` file
* the column name of the OD you want to compute growth rate for:
  - tends to be `OD_A`, `OD_B`, `OD_C`, etc.
* the length of the smoothing window:
  - an integer, e.g. 5 represents 5 values in your time resolution
* the smoothing method:
  - "gaussian": gaussian smoothing
  - "avg": rolling average smoothing
  - "optimise": try out both and see which yields the highest R2
* the maximum number of change points:
  - an integer representing the maximum number of slope change points the function will look for
  - "optimise": try out a bunch of values and use the on that yields the highest R2 (recommended)

Calling the main function will return estimates of:
* doubling rate
* doubling time
* growth rate  
  
A figure representing the estimation will also appear.

# 2. Algorithm
The algorithm uses the matlab [`ischange`](https://uk.mathworks.com/help/matlab/ref/ischange.html) function to look for changes in linear slope in the smoothed log(OD).  
It fits a linear line to the points between the first two change points and uses the slope of that line to estimate doubling rate, doubling time, and growth rate.
OD values smaller than 0.05 are dropped, since this represents the detection threshold of the Bioreactor.


# 3. Questions
Feel free to email me if you have questions: [wolf.de.wulf@ed.ac.uk](mailto:wolf.de.wulf@ed.ac.uk)