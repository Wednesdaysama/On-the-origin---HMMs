# On the origin - HMMs
We provide the hidden Markov models (HMMs) corresponding to the different functions we defined based on the two phylogenetic trees.

### Workflow
15 protein sequence sets (in FASTA format) were collected. These included 10 **inferred** functions: putative sodium-transporting P-type ATPases, Potassium-transporting_ATPase_ATP-binding_subunit, Magnesium-transporting_ATPase, Zinc_cadmium_lead_cobalt-transporting_P-type_ATPase, Copper-exporting_P-type_ATPase, MrpA, and MrpD (the first ring of each phylogenetic tree). We also **selected** five alkaline-enriched subfamilies, namely putative sodium-transporting P-type ATPases, MrpA, MrpA, MrpD, and MrpD (from the third ring of each phylogenetic tree).

For each protein set, we reduced redundancy using MMseqs2 with parameters *-c 0.8 --cov-mode 0 --min-seq-id 0.5*. 

  mkdir mmseq2
  cd mmseq2

Loop:

  for fasta in ../*.fasta
  do
      name=$(basename "$fasta" .fasta)
      echo "===== Processing $name ====="
      mmseqs createdb "$fasta" "./${name}.mmseqdb"
      mmseqs cluster \
          -c 0.8 \
          --cov-mode 0 \
          --min-seq-id 0.5 \
          "./${name}.mmseqdb" \
          "./${name}.clustering" \
          "./tmp_${name}"
      mmseqs createtsv \
          "./${name}.mmseqdb" \
          "./${name}.mmseqdb" \
          "./${name}.clustering" \
          "./${name}.clustering.tsv"
      mv "./${name}.clustering.tsv" ../
      echo "===== Finished $name ====="
  done



We retained representative sequences and aligned them with Clustal Omega, then built HMMs with hmmbuild. We then evaluated the resulting models and report the summary statistics.

### HMMs
[inferred_sodium_ATPase.hmm](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/inferred_sodium_ATPase.hmm): HMM for the entire inferred putative sodium-transporting P-type ATPases branch of the tree (Figure 3, first ring, orange).

[inferred_MrpA.hmm](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/inferred_MrpA.hmm): HMM for the entire inferred MrpA branch of the tree (Figure 4, first ring, orange).

[inferred_MrpD.hmm](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/inferred_MrpD.hmm): HMM for the entire inferred MrpD branch of the tree (Figure 4, first ring, red).

[selected_sodium_ATPase.hmm](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/selected_sodium_ATPase.hmm): HMM for the entire selected putative sodium-transporting P-type ATPases branch of the tree (Figure 3, third ring, blue).

[selected_MrpA.hmm](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/selected_MrpA.hmm): HMM for the alkaline-enriched MrpA subfamilies branch of the tree (Figure 4, third ring, light blue).

[selected_MrpA_asterisk.hmm](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/selected_MrpA_asterisk.hmm): HMM for the alkaline-enriched MrpA* subfamilies branch of the tree (Figure 4, third ring, dark blue).

[selected_MrpD.hmm](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/selected_MrpD.hmm): HMM for the alkaline-enriched MrpD subfamilies branch of the tree (Figure 4, third ring, brown).

[selected_MrpD_asterisk.hmm](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/selected_MrpD_asterisk.hmm): HMM for the alkaline-enriched MrpD* subfamilies branch of the tree (Figure 4, third ring, purple).

### HMM report

| Model | seq | Eff_nseq | M | relent | info | p |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| inferred_MrpA.hmm | 477 | 49.38 | 1508 | 0.59 | 0.63 | 0.37 |
| inferred_MrpD.hmm | 453 | 18.54 | 1225 | 0.59 | 0.63 | 0.44 |
| inferred_sodium_ATPase.hmm | 275 | 15.30 | 1712 | 0.59 | 0.63 | 0.50 |
| selected_MrpA.hmm | 61 | 3.01 | 777 | 0.59 | 0.61 | 0.52 |
| selected_MrpA_asterisk.hmm | 22 | 5.34 | 1141 | 0.59 | 0.61 | 0.51 |
| selected_MrpD.hmm | 29 | 2.71 | 498 | 0.59 | 0.62 | 0.54 |
| selected_MrpD_asterisk.hmm | 106 | 6.14 | 1058 | 0.59 | 0.63 | 0.49 |
| selected_sodium_ATPase.hmm | 25 | 3.60 | 1912 | 0.59 | 0.62 | 0.52 |

### Test
Demo data can be found [here](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/sample.fasta). 

Results can be found [here](https://github.com/Wednesdaysama/On-the-origin---HMMs/blob/main/results.txt). 

