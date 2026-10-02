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

Based on the true-positive and false-positive scores, select a cutoff for each HMM. Make sure the true positive rate (TPR) >= 80%, the false negative rate (FNR) <= 20%, and the precision rate >= 90%, except for the selected_MrpA_ast HMM. Please refer to the [HMM report](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/Results/hmm_cutoff_analysis.xlsx) for the complete cutoff distribution.

| model | score_type | N_positive | N_negative | selected_cutoff | TP | FP | TN | FN | TPR | FPR | FNR | precision |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| selected_MrpA | relative | 363 | 18735 | 167.2 | 314 | 0 | 18735 | 49 | 86.5% | 0.0% | 13.5% | 100.0% |
| selected_MrpA_ast | relative | 430 | 18668 | 76.8 | 341 | 38 | 18630 | 89 | 79.3% | 0.2% | 20.7% | 90.0% |
| selected_MrpD | relative | 357 | 18741 | 162 | 310 | 1 | 18740 | 47 | 86.8% | 0.0% | 13.2% | 99.7% |
| selected_MrpD_ast | relative | 309 | 18789 | 279.2 | 275 | 6 | 18783 | 34 | 89.0% | 0.0% | 11.0% | 97.9% |
| selected_Putative_sodium-transporting_P-type_ATPase | relative | 455 | 18643 | 268.9 | 383 | 0 | 18643 | 72 | 84.2% | 0.0% | 15.8% | 100.0% |


For each sequence, the relative score was calculated as the difference between the bit score of the target HMM and the highest bit score among all competing HMMs. A sequence was classified as the corresponding selected subfamily when its relative score was greater than or equal to the model-specific relative score cutoff.


