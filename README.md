# C57104wMCM
## Title: Transformer DGA Screening - IEEE C57.104-2019 with Propagated Uncertainty via a Monte Carlo Method
### Overview:
Please refer to the [provided images](Plots_Example_output.pdf) to look at the validation plots for a better understanding of what is happening.
I will link to my Thesis draft once it gets submitted for added details - not sure I'm allowed to do so now.

The general flow internally is: C57104wMCM -> preprocess -> sample -> screening -> diagnosis.
I tried to provide documentation even for the internal functions, but they are not really intended for direct use.

The function will take in all given gas samples, generate random values based on the statistical description of the uncertainty provided, and then apply a simplified implementation of the IEEE C57.104-2019 Methodology. It is assumed only appropriate samples meant for comparison against the relevant IEEE Tables are provided. Similarly, it is assumed that the appropriate limits for said Tables are provided. This function ought to be "wrapped" around with additional logic to ensure these assumptions are met, else it is baseless. 

Using 'validate = TRUE' will output various plots that can help with validation. These are done using the packages: [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html), [GGally](https://cran.r-project.org/web/packages/GGally/index.html), and [ternary](https://cran.r-project.org/web/packages/Ternary/index.html). Please note, plotting can be quite slow.

### Sampling:
The function 'sample()' is responsible for the generation of random samples. This is primarily done via R's "rnorm" and "mvrnorm" from the package, [MASS](https://cran.r-project.org/web/packages/MASS/index.html). There is an allowance for "rho" or a correlation matrix to be inputted to link distributions across gases for the samples. There is no such allowance for distributions to be linked across samples for the gases. This is quite contentious (for anyone that might care) as there is a strong case for its inclusion. However, I felt underskilled to handle to statistical implications as there are many problematic edge-cases to consider. If a correlation matrix is provided, the "flex" argument allows for it to be adjusted to ensure it is a positive definite, otherwise, the function cannot run if given a problematic correlation matrix. This is done via the "make.positive.definite" function from the package, [corpcor](https://cran.r-project.org/web/packages/corpcor/index.html). The supported distribution shapes are normal, triangular, and uniform. I use the "dtriangle" function from the package, [triangle](https://cran.r-project.org/web/packages/triangle/index.html) to help with some plots. To implement correlations for the latter two, I used the idea sourced from [Glen_b](stats.stackexchange.com/q/143280). A multivariate normal distribution with given correlations is first generated, and then converted to either a triangular or uniform distribution. This is an approximation. I am unconvinced that correlation across gases should be added given the uncertainty surrounding the appropriate values. However, I would suggest investigating "rho" = 0 and "rho" = 1 as a form of sensitivity analysis of the outputs.

Uncertainty is represented primarily via two arguments: "accuracy" and "min_accuracy". The first is multiplied by the given gas sample to provide a +/- estimate. If the chosen distribution shape is normal, it is interpreted as a 95% probability level. Otherwise, it is interpreted as a 100% probability level. This driven by my interpretation of the Normative Standards feeding into IEEE C57.104-2019, although it does make things less intuitive, unfortunately. The "min_accuracy" is treated as +/- PPM directly, and will override the calculated value for "accuracy" if greater. Either of these can be omitted to force the use of the other. If unsure, perhaps use 0.15 for "accuracy", 1 for "min_accuracy", and a normal distribution. However, know that these are crude estimates. A better implementation would have allowed for scaling uncertainty (the "IEC Specification" would suggest 0.15 -> 0.30 as the gas values approach their analytical limit). Similarly, it would have been better to allow for "min_accuracy" to scale up to a given value to avoid any negative samples being generated.

If using 'validate = TRUE', then the plots associated with the 'sample()' function include:
* 1 x Plot of gas values (PPM) versus dates since the oldest sample was taken. Plots are seperated into H2+CnHn gases and COn gases.
* 1 x Plot of gas values (PPM) versus dates since the oldest sample was taken, with error bars indicated Accuracy. Plots are seperated into individual gases. If a given error bar is in red, then it indicates the 'min_accuracy' was superceded the 'accuracy'.
* d x Correlation Plots of generated gas values (PPM), where there is one plot for each day, d. These will display a histogram of the generated samples with a superimposed ideal distribution shape to ensure marginal distributions are expected. It also includes the correlation values to ensure the correlations are as expected.

### Screening:
The function 'screening()' is responsible for applying the IEEE C57.104-2019 Methodology to each of the generated random samples. The frequency of given occurances, such as the number of times a gas fails its Table 1 relative to the total number of attempts, is used as a proxy for the probability of given events. If only one day's samples are available, then only IEEE Tables 1 and 2 are applicable. If two days' samples, then Tables 1, 2, and 3 are applicable. If more, then all are used for all Tables. The methodology stipulates only 3-6 samples be used. This is ignored! If you wish this enforced, simply input only the samples that are deemed appropriate for the Screening. The estimated likelihood for each the DGA Status for each gas, and combined, are outputted from this function. This can be used as an indicator for confidence in the Screening output. 

If using 'validate = TRUE', then the plots associated with the 'screening()' function include:
* g x Correlation Plots of IEEE metrics derived from generated gas values, where there is one plot for each gas, g. These will display a histogram of the generated samples colour coded based on the given gas's DGA Status. It also includes the correlation values between the metrics, and boxplots of the metric divided by the given gas's DGA Status. Lastly, it includes a swarmplot of the values divided by the given gas's DGA Status. In the future, the limits should perhaps be drawn explictly here for added clarity. Nevertheless, this visual should provide good insight as to what metric / Table is driving the given gas's DGA Status.
* 1 x 'Correlation' Plot of Per-Gas DGA Statuses and Combined DGA Status. This shows the proportion of DGA Statuses and their co-occurances. Also highlights one of the issues I have with the Methodology: using worst-case Per-Gas DGA Status as the Combined DGA Status is very onerous. 

### Diagnosis:
The function 'diagnosis()' is responsible for applying some basic diagnostics to each of the generated random samples. From the IEEE C57.104-2019, Duval Triangles 1, 4, and 5 are included, as well as Rogers Ratio Method. For additional reference, IEC 60599:2022 Tables 1 and 2 are also included. The 'diagnosis' argument can be used to select which of this (if any) to run. If choosing to apply diagnostics, H2 + all CnHn gases must be present. A simplified interpretation of the IEEE C57.104-2019 Methodology is that Diagnostics should be applied only if the outputted DGA Status level is greater than 1 (i.e. 2 or 3). The current implementation provides several variations of the outputs:
* Results using all samples.
* Results using only samples with DGA Status 1.
* Results using only samples with DGA Status not 1.
* (The source code has commented out the required code to also process DGA Status 2, and DGA Status 3 seperately).

If using 'validate = TRUE', then there will be plots associated with the 'diagnosis()' function include. Plots are only included if the relevant method was selected via the 'diagnosis' argument. Similarly, the variant of a given plot of will be omitted if no samples are applicable. The potentially-included plots are as follows:
* Duval Triangles (1|4|5) Ternary Plots. The plot contains black dots representing each of the generated gas values, and a white asterix representing the measured gas sample. These can unfortunately be quite hard to see and should be developed further for clearer visuals. A contour plot is also included to indicate the distribution densities.
* Duval Triangles (1|4|5) Summary Plots. The plot shows the relative frequency of each given diagnosis. Please note: the data here is different to that of the above! Only when applicable are the outputs from Duval Triangles 4 and 5 included here. This discrepency is a potential cause for confusion, sorry!
* Rogers Ratio Method Projected Plot. This is a composite plot where the 3D space is projected onto 2x2D space. The plot(s) contain black dots representing each of the generated gas values, and a white asterix representing the measured gas sample. Logarthmic scales are used, and the space is truncated for clarity despite really being unbound.  
* Rogers Ratio Method Summary Plot. The plot shows the relative frequency of each given diagnosis.
* IEC 60599:2022 Tables (1|2) Projected Plot. This is a composite plot where the 3D space is projected onto 2x2D space. The plot(s) contain black dots representing each of the generated gas values, and a white asterix representing the measured gas sample. Logarthmic scales are used, and the space is truncated for clarity despite really being unbound.
* IEC 60599:2022 Tables (1|2) Summary Plots. The plot shows the relative frequency of each given diagnosis.
