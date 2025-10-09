# sskmeans
This is my Masters Thesis project on _Accommodating Positive & Negative Constraints in the K-Means Algorithm_. 

ABSTRACT:  K-Means is one of the most widely used clustering algorithms, but its standard form is an unsupervised learning method that assumes no prior information about data labels. In practice, partial knowledge often exists in the form of constraints: some points must be grouped together (positive constraints), while others must be kept apart (negative constraints). Positive constraints are relatively easy to incorporate due to their transitive nature, but negative constraints introduce complex dependencies that are much more difficult to accommodate. This project explores modifications to the K-Means algorithm that integrate both types of constraints, with the goal of improving clustering accuracy in semi-supervised settings.

The result of my efforts will be an `R` package with a main function of performing this modified calculation of the K-Means algortihm, called `sskmeans`. 

So far, data must be entered into the function with the following format:

 - data must be numeric, with the exception of columns `pos.eq` (positive equivalence constraints) and `neg.eq` (negative equivalence constraints)
 - `pos.eq` and `neg.eq` must be specified as names of existing columns of the data, otherwise defaults to _NULL_
 - `pos.eq` is an integer column of block labels, indicating membership of INDIVIDUAL OBSERVATIONS to the same cluster
 - `neg.eq` must be a list column indicating INDIVIDUAL OBSERVATIONS that should be in a separate cluster than the row indicated

I have included a toy example dataframe to illustrate this set-up at the bottom of the function file.
