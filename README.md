# sskmeans
This is my Masters Thesis project on _Accommodating Positive & Negative Constraints in the K-Means Algorithm_. 

ABSTRACT:  K-Means is one of the most widely used clustering algorithms, but its standard form is an unsupervised learning method that assumes no prior information about data labels. In practice, partial knowledge often exists in the form of constraints: some points must be grouped together (positive constraints), while others must be kept apart (negative constraints). Positive constraints are relatively easy to incorporate due to their transitive nature, but negative constraints introduce complex dependencies that are much more difficult to accommodate. This project explores modifications to the K-Means algorithm that integrate both types of constraints, with the goal of improving clustering accuracy in semi-supervised settings.

The result of my efforts will be an `R` package with a main function of performing this modified calculation of the K-Means algortihm, called `sskmeans(data, k, pos.eq = NULL, neg.eq = NULL, ...)`. 

So far, data must be entered into the function with the following format:

 - data entered must be numeric
 - `pos.eq` (positive equivalence constraints) is a BLOCK-LEVEL list indicating membership of INDIVIDUAL OBSERVATIONS (row numbers) to the same cluster
 - `neg.eq` (negative equivalence constraints) is an adjacency matrix at the BLOCK-LEVEL indicating `1` for pairwise negative constraints; `0` otherwise (i.e. no negative constraint)
 - `pos.eq` and `neg.eq` both default to _NULL_; function will run with either one, both, or neither

I have included a toy example dataframe to illustrate this set-up at the bottom of the function file.
