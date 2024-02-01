This is an R implementation of the stable Sr reduction derivation outlined in Krabenhoft et al. (2009) for measurements made on a TIMS with an 87Sr/84Sr double spike. 

reduce_functions.R is a script containing two functions: one which performs the actual calculations described in Krabbenhoft et al. (2009), and a second which applies the first function to a batch of measurements for fast data processing.

Batch_Reduction_Template.Rmd is an RMarkdown document demonstrating how to apply the two functions described above on sample data.

Two notes of caution:
- the tidyverse package must be installed for these functions to work
- the efficacy of these functions is very dependent on proper organization of the input data

For a full explanation of the derviation please refer to:
Krabbenhoft, A., Fietzke, J., Eisenhauer, A., Liebetrau, V., Bohm, F., Vollstaedt, H. (2009) Determination of radiogenic and stable strontium isotope ratios (87Sr/86Sr; δ88/86Sr) by thermal ionization mass spectrometry applying an 87Sr/84Sr double spike. Journal of Analytical Atomic Spetrometry 24(9):1267-1271. doi: 10.1039/B906292K.











