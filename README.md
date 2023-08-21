# Unfolding

For properly filling of response matrix, gen-level and reco-level hist
Check Anal_Leptop_PROOF.C code

The sequence in which you can try is

Loop over all gen-level objects
  - Fill the gen-level hist for all gen-level object
 - Fill the response matrix first for the cases for which you find gen level obejct but not reco level object (Take gen object quantity as -1 or -100) (Check Line 1896)
   
Loop over all reco-level objects
 - Fill the reco  for all reco-level objects
- Fill the response matrix first for the cases for which you find reco obejct and gen object (Check Line 2620)
 -Fill the response matrix first for the cases for which you find reco obejct but not gen object - (Take gen object quantity as -1 or -100) (Check Line 2623)


You can also check for the simple case with no fake rate and efficiency correction required first. In this case, just fill all 3 hists for the cases you have both reco and gen object.


For efficiency and fake rate correction

Check make_input_histforunfolding function in unfolding_function.C (https://github.com/S-ritik/Unfolding/blob/main/unfolding_function.C#L1543-L1613)

