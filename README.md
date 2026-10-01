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

#### 4. According to the true positive and false positive sequences, calculate the cutoff for each HMM


#### 5. Data summary

    python ../hmm_results_new.py


[Results](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/Results/hmm_validation_scores_0.8.xlsx)

| model | TP | FP | TN | FN | precision | recall | false_positive_rate | false_negative_rate | lowest_true_score | highest_false_score | lowest_true_relative_score | highest_false_relative_score | relative_score_separation | relative_score_cutoff |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| selected_MrpA | 37 | 4 | 1876 | 0 | 90.2% | 1 | 0.2% | 0 | 581.9 | 1011.7 | 182.9 | 205 | -22.1 | |
| selected_MrpA_ast | 43 | 13 | 1861 | 0 | 76.8% | 1 | 0.7% | 0 | 452.4 | 711.2 | 45.4 | 55.4 | -10 | |
| selected_MrpD | 36 | 12 | 1869 | 0 | 75.0% | 1 | 0.6% | 0 | 655.6 | 659.2 | 152.4 | 70.3 | 82.1 | 111.35 |
| selected_MrpD_ast | 31 | 1 | 1885 | 0 | 96.9% | 1 | 0.1% | 0 | 434.8 | 481.8 | 204.1 | 48.7 | 155.4 | 126.4 |
| selected_Putative_sodium-transporting_P-type_ATPase | 46 | 17 | 1854 | 0 | 73.0% | 1 | 0.9% | 0 | 1470.6 | 1522.5 | 189.4 | 119.9 | 69.5 | 154.65 |

For each sequence, the relative score was calculated as the difference between the bit score of the target HMM and the highest bit score among all competing HMMs. A sequence was classified as the corresponding selected subfamily when its relative score was greater than or equal to the model-specific relative score cutoff.

There is no clear cutoff for distinguishing selected_MrpA and selected_MrpA_ast from inferred_MrpA. 
Notably, all false positives for selected_MrpA and selected_MrpA_ast belonged to inferred_MrpA. 
Therefore, we performed a sliding-window optimization to identify a region that can discriminate selected_MrpA (and selected_MrpA_ast) from inferred_MrpA.

The optimization was performed separately for selected_MrpA versus inferred_MrpA and selected_MrpA_ast versus inferred_MrpA. For each comparison, sequences in the HMM-building sets of the target and inferred_MrpA classes were combined and aligned using Clustal Omega. Sliding windows of 100, 150, 200, 250, and 300 alignment columns were evaluated at 25-column intervals. Windows with excessive gaps or insufficient sequence coverage were excluded before HMM construction. For each retained window, separate profile HMMs were constructed from the target and inferred_MrpA training sequences using HMMER and evaluated against the corresponding independent test sequences.

[Results of selected_MrpA_ast vs inferred_MrpA](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/Results/selected_MrpA_ast_vs_MrpA_optimization.xlsx)

| rank | region | window_size | cutoff | TP | FP | TN | FN | precision | recall | FPR | F1 | lowest_true_relative_score | highest_false_relative_score | relative_score_separation |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | 1251_1550 | 300 | 14.05 | 43 | 0 | 205 | 0 | 1 | 1 | 0 | 1 | 21.2 | 6.9 | 14.3 |
| 2 | 1326_1575 | 250 | 15.6 | 43 | 0 | 205 | 0 | 1 | 1 | 0 | 1 | 22.7 | 8.5 | 14.2 |
| 3 | 1301_1550 | 250 | 14.2 | 43 | 0 | 205 | 0 | 1 | 1 | 0 | 1 | 21.2 | 7.2 | 14 |
| 4 | 1301_1600 | 300 | 14.4 | 43 | 0 | 205 | 0 | 1 | 1 | 0 | 1 | 19 | 9.8 | 9.2 |
| 5 | 1276_1525 | 250 | 16.45 | 43 | 0 | 205 | 0 | 1 | 1 | 0 | 1 | 20.3 | 12.6 | 7.7 |

[Results of selected_MrpA vs inferred_MrpA](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/Results/selected_MrpA_vs_MrpA_optimization.xlsx)

| rank | region | window_size | cutoff | TP | FP | TN | FN | precision | recall | FPR | F1 | lowest_true_relative_score | highest_false_relative_score | relative_score_separation |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | 1526_1675 | 150 | 21.95 | 37 | 0 | 205 | 0 | 1 | 1 | 0 | 1 | 22.1 | 21.8 | 0.3 |
| 2 | 1526_1825 | 300 | 25.15 | 37 | 1 | 204 | 0 | 0.973684211 | 1 | 0.004878049 | 0.986666667 | 25.4 | 26.7 | -1.3 |
| 3 | 1551_1650 | 100 | 11.35 | 37 | 1 | 204 | 0 | 0.973684211 | 1 | 0.004878049 | 0.986666667 | 12.3 | 15 | -2.7 |
| 4 | 1676_1925 | 250 | 17.95 | 37 | 1 | 204 | 0 | 0.973684211 | 1 | 0.004878049 | 0.986666667 | 21.8 | 24.8 | -3 |
| 5 | 1526_1625 | 100 | 11.25 | 37 | 1 | 204 | 0 | 0.973684211 | 1 | 0.004878049 | 0.986666667 | 12.2 | 15.3 | -3.1 |


For selected_MrpA_ast, the optimal window corresponded to alignment columns 1251–1550. The lowest relative score among selected_MrpA_ast sequences was 21.2, whereas the highest relative score among inferred_MrpA sequences was 6.9, giving a relative-score cutoff of **14.05**. For selected_MrpA, the optimal window corresponded to alignment columns 1526–1675, with a lowest true relative score of 22.1 and a highest false relative score of 21.8, giving a cutoff of **21.95**. Sequences with relative scores greater than or equal to the corresponding cutoff were classified as selected_MrpA_ast or selected_MrpA, respectively; sequences below the cutoff were classified as inferred_MrpA.
