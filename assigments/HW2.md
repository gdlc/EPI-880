**Homework 2**

In Example 1.2 we evaluate the effects of the number of markers used on prediction accuracy.  In this example, models are build by 
first ranking markers based on GWAS and then estimating effects using OLS for the top 1, 2,….,300 markers. 
Each of these models are then evaluated based on their prediction accuracy in the testing set. 

We observed that prediction accuracy grows up to a point and then declines, probably due to increase variance of the OLS estimates.
In our second module we discussed Ridge Regression as a way to reduce the variance of estimates. In this HW we looked at combining both variable selection
and shrinkage estimation. 

Assignment:

- Replace the OLS estimation with RR estimates (for that you can use the code below)
- Produce a plot with R-squared versus number of markers with results obtained using both OLS and Ridge Regression. 
- Submit by Wednesday July 20th, 10:00am  a ppt or a pdf with your final plot and your code.

Note: I’ll randomly draw a name and this person will need to explain what he/she did and discuss the results.
