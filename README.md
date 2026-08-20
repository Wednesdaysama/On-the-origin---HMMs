# On the origin - HMMs
We provide the hidden Markov models (HMMs) corresponding to the different functions we defined based on the two phylogenetic trees.

### Workflow
Eight protein sequence sets were collected. These included three **inferred** putative sodium-transporting P-type ATPases, MrpA, and MrpD (the first ring of each phylogenetic tree). We also included five alkaline-enriched subfamilies we **selected**: putative sodium-transporting P-type ATPases, MrpA, MrpA*, MrpD, and MrpD* (the third ring of each phylogenetic tree).

For each protein set, we clustered these sequences using MMSeq2 with (-c 0.8 --cov-mode 0 --min-seq-id 0.5) to reduce redundancy. Representative sequences were kept and aligned using Clustal Omega. Then, *hmmbuild* was used to construct HMMs. Finally, we checked these models and report the results.

### HMMs
[inferred_MrpA.hmm](): HMM for the entire inferred MrpA section of the tree.



