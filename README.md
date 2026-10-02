# On the origin - HMMs
We provide the hidden Markov models (HMMs) corresponding to the different functions we defined based on the two phylogenetic trees.

### Workflow
15 protein sequence sets (in FASTA format) were collected. These included 10 **inferred** functions: putative sodium-transporting P-type ATPases, Potassium-transporting_ATPase_ATP-binding_subunit, Magnesium-transporting_ATPase, Zinc_cadmium_lead_cobalt-transporting_P-type_ATPase, Copper-exporting_P-type_ATPase, MrpA, and MrpD (the first ring of each phylogenetic tree). We also **selected** five alkaline-enriched subfamilies, namely putative sodium-transporting P-type ATPases, MrpA, MrpA, MrpD, and MrpD (from the third ring of each phylogenetic tree).


#### 1. Remove incomplete sequences from the raw protein sets.

    mkdir -p complete
    for fasta in ./*.fasta
    do
        name=$(basename "$fasta")
        awk '
        /^>/ {
            keep = (tolower($0) !~ /partial/)
        }
        keep {
            print
        }
        ' "$fasta" > "./complete/$name"
        echo "Finished: $name"
    done
    


#### 2. For each .fasta in *complete*, align them with Clustal Omega, then build HMMs with hmmbuild (took a while). 

    for fasta in ./*.fasta
    do
        name=$(basename "$fasta" .fasta)
        echo "===== processing $name ====="
        clustalo -i "$fasta" -o "./${name}.aligned.fasta" --force -v
        hmmbuild "./hmm/${name}.hmm" "./${name}.aligned.fasta"
        echo "===== finished $name ====="
    done

#### 3. Run HMMs against all sequences
    mkdir -p ./hmm_results
    for fasta in ./*.fasta
    do
        [[ "$fasta" == *aligned* ]] && continue
        test_name=$(basename "$fasta" .fasta)
        for hmm in ./hmm/*.hmm
        do
            hmm_name=$(basename "$hmm" .hmm)
            echo "===== $test_name vs $hmm_name ====="
            hmmsearch \
                --tblout "./hmm_results/${test_name}__vs__${hmm_name}.tbl" \
                "$hmm" "$fasta" > "./hmm_results/${test_name}__vs__${hmm_name}.out"
        done
    done

#### 4. Results summary and cutoff distribution

    python ./hmm_results.py

Based on the true-positive and false-positive scores, select a cutoff for each HMM. Make sure the selected cutoff has the lowest false positive rate (FPR) plus the false negative rate (FNR) , and that the cutoff > 0, except for the selected_MrpA_ast HMM. Please refer to the [HMM report](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/Results/hmm_cutoff_analysis.xlsx) for the complete cutoff distribution.

| model | score_type | N_positive | N_negative | cutoff | TP | FP | TN | FN | TPR | FPR | FNR | FPR+FNR |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| selected_MrpA | relative | 363 | 18735 | 1.6 | 362 | 48 | 18687 | 1 | 99.72% | 0.26% | 0.28% | 0.53% |
| selected_MrpA_ast | relative | 430 | 18668 | 0.4 | 413 | 208 | 18460 | 17 | 96.0% | 1.1% | 4.0% | 5.07% |
| selected_MrpD | relative | 357 | 18741 | 85.1 | 348 | 2 | 18739 | 9 | 97.5% | 0.0% | 2.5% | 2.53% |
| selected_MrpD_ast | relative | 309 | 18789 | 61.3 | 308 | 17 | 18772 | 1 | 99.7% | 0.1% | 0.3% | 0.41% |
| selected_Putative_sodium-transporting_P-type_ATPase | relative | 455 | 18643 | 28.6 | 455 | 103 | 18540 | 0 | 100.0% | 0.6% | 0.0% | 0.55% |


For each sequence, the relative score was calculated as the difference between the bit score of the target HMM and the highest bit score among all competing HMMs. A sequence was classified as the corresponding selected subfamily when its relative score was greater than or equal to the model-specific relative score cutoff.


