# KnownDatedAll.csv

This dataframe is an output where each row is known recombinant
block assesed and dated by HybRIDS.

This dataframe contains the following important columns:

**Triplet:**
Names of the three sequences contained in the triplet the
block was found in.

**Sequence_Pair:**
The two sequences in the triplet between which the
recombination block was found.

**Sequence_Similarity_Threshold:**
The seuqence similarity threshold the block beat to be detected as
recombinant by HybRIDS.

**First_BP_Position:**
Approximate base position where the block starts.

**Last_BP_Position:**
Approximate base position where the block ends.

**Approximate_Length_BP:**
Approximate length of the recombinant block in bp.

**Number_of_SNPs:**
Number of mutations in the block between the two sequences.

**Corrected_Number_of_SNPs:**
Same as above but corrected, for example by a JC model.

**ninetyFiveAge:**
The 95% age estimate for the given block.

**fiftyAge:**
The 50% age estimate for the given block.

**fiveAge:**
The 5% age estimate for the given block.

**P_Value:**
The p-value of the given block.

**P_Thresh:**
The p-value threshold the block beat to be confirmed as
significant and reported by HybRIDS.

**ActualAge:**
The actual amount of time (in mu * generations) since the recombination
event that generated the block occured.

**ActualAgeBefore:**
The actual amount of time (in mu * generations) that the parental sequences
in the triplet diverged before coming together to form the recombinant sequence.