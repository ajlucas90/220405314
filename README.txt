README for MATLAB codes for 2204.05314.

======
DISCLAIMER: This code was written for personal use.  I would arguably say that the codes needed to simulate this model are simple enough it would be beneficial for a user to write them for themselves. In case you want my codes (e.g. to inspect if I made a mistake when coding this model) I am providing them here.  To avoid proliferating comments through a bunch of different files, many of which (out of laziness) define and call the same or similar functions, I will just remark on important variable names and summarize function purposes here.

======
V9 VS. V10:

v10 codes are updated to handle the modeling of production costs (SI Appendix 6) and buyer/seller noise (SI Appendix 7).  If you look at the v9 codes, you will see that some of these ideas are implemented already in them...but in a lousy way (which I ended up not using).   It is recommended that you use the v10 codes.  

I include some v9 codes which have slightly modified scripts or calls which I have not yet included in v10 codes because I never used them again.  I have explicitly checked (that to the best of my knowledge) all the scripts and methods included here work as intended - the models have been run many times and give sensible results.  While it should be straightforward to cut/paste from v9 into v10, in case that process introduced any unintended bugs I haven't done this.

======
COMMON VARIABLE NAMES (not in all scripts):

N = number of buyers

M = number of sellers (be careful setting M = 1 or 2, certain normalizations in some scripted files can break!)

Nt = number of time steps the simulation will run

mu = parameter from main text of the same letter

sigma = parameter from main text of the same letter

J = the parameter \bar J (up to an overall factor of sigma) from the main text in the manyruns scripts (also called Jrel there); in onerun script, the unnormalized parameter J. 

kappa = \bar\rho from the main text.

waitprob = the probability that a buyer does NOT act in a given turn. in onerun scripts this is used to define kappa directly (as should be obvious from the way it is coded).

a = parameter in SI Appendix 6 of the same letter, in manyruns scripts called arel.  In onerun script, a has been multiplied by M (but it should be obvious from context). (v10 only: set to 0 if using v9)

b = parameter in SI Appendix 6 of the same letter. (v10 only: ignore in v9)

beta = parameter of the same letter in SI Appendix 7. (v10 only)

sellernoise = parameter called eta in SI Appendix 7. (v10 only: set to 0 if using v9)

sellerprob = parameter called f in SI Appendix 7. (v10 only)

numruns = in manyruns files, the number of independent simulations to run. 

======
FUNCTIONS:

manyruns_fast_v10:  Bread and butter runs with all parameters except beta.  Returns 8 arrays with data from each run on aggregate market properties.  "meanQarray" is the average over t of the parameter M/Q defined in the main text.  "varQarray" is the standard deviation over t of M/Q.  "varmeanqarray" is the standard deviation over seller i of the mean in t of q_i(t).  "meanvarqarray" is the mean over seller i of the standard deviation over t in q_i(t). "avgprofitarray" reports the i/t averaged pi_i(t).  "avgparray" reports the i/t averaged p_i(t).  "avgqarray" reports the i/t averaged q_i(t) (i.e. 1-q_0).  "fliprate" reports the parameter gamma defined in the main text.  

manyruns_fast_v10beta: Same as v10 except also includes beta.

manyruns_fast_v10shocks: Used for the exogenous model in SI Appendix 8.  "shocktime" is the number of time steps between each shock.  "tt0" is the t parameter defined in SI Appendix 8 for correlation functions.  "meanq0"/"stdq0" are mean/std of C(t). "rmeanq0"/"rstdq0" are mean/std of R(t).

manyruns_fast_v9_distq: Used to generate Fig. 4 and related in the main text.  [qi,Pqi] is log-log-scale histogram for the probability distribution function P(q_i), normalized as a (discretized) continuum PDF. In this script, data is only collected at times when the seller with highest market share changes.

manyruns_fast_v9_distqtyp: Same as above, except data is collected at all times.

manyruns_fast_v9_SBtimescale: Used to generate Figure S5 in the SI.  Calculate the time for a monopolist who overprices good to lose market share (assuming they can never re-choose price).  This script is used to demonstrate that the instability of SB to NE is associated with early time dynamics, not late time dynamics.

onerun calls do what they say -- just run one instance of the simulation!  The "selleroptimizationplot" is used to visualize the seller's optimization problem as in Figure S7 in the SI.  These are useful for getting quick results, and the terminology follows the scripts to run many times.  You can often see a few different ways of plotting things or collecting extra data commented out here.