import os
import glob
import numpy as np
import pandas as pd

input_dir = "./test/hmm_results"
output_file = "./test/hmm_validation_scores.xlsx"

records = []


for tbl_file in glob.glob(os.path.join(input_dir, "*.tbl")):
    basename = os.path.basename(tbl_file)
    name = basename.removesuffix(".tbl")
    if "__vs__" not in name:
        print(f"Skipping unexpected filename: {basename}")
        continue
    true_class, hmm_model = name.split("__vs__", 1)
    with open(tbl_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            sequence = parts[0]
            evalue = float(parts[4])
            score = float(parts[5])
            records.append({
                "sequence": sequence,
                "true_class": true_class,
                "hmm_model": hmm_model,
                "score": score,
                "evalue": evalue
            })

df = pd.DataFrame(records)
print(f"Total HMM hits read: {len(df)}")
if df.empty:
    raise ValueError(
        f"No HMM hits were found in {input_dir}"
    )
score_matrix = df.pivot_table(
    index=["sequence", "true_class"],
    columns="hmm_model",
    values="score",
    aggfunc="max"
)

comparison_matrix = score_matrix.fillna(
    float("-inf")
)
hmm_models = list(score_matrix.columns)
results = []
for (sequence, true_class), row in comparison_matrix.iterrows():
    sorted_scores = row.sort_values(ascending=False)
    best_model = sorted_scores.index[0]
    best_score = sorted_scores.iloc[0]
    if len(sorted_scores) > 1:
        second_model = sorted_scores.index[1]
        second_score = sorted_scores.iloc[1]
        delta_score = (best_score - second_score)
    else:
        second_model = None
        second_score = np.nan
        delta_score = np.nan
    if np.isneginf(best_score):    # If everything is -inf, no HMM actually produced a hit
        best_model = None
        best_score = np.nan
        second_model = None
        second_score = np.nan
        delta_score = np.nan
        correct = False
    else:
        if np.isneginf(second_score):
            second_score = np.nan
            delta_score = np.nan
        correct = (best_model == true_class)
    results.append({
        "sequence":sequence,
        "true_class":true_class,
        "best_model":best_model,
        "best_score":best_score,
        "second_model":second_model,
        "second_score":second_score,
        "delta_score":delta_score,
        "correct_classification":correct
    })
classification_df = pd.DataFrame(results)
score_output = score_matrix.reset_index()
final_df = classification_df.merge(score_output, on=["sequence", "true_class"], how="left")

for model in hmm_models:
    relative_col = (
        f"{model}_relative_score"
    )
    model_score = comparison_matrix[
        model
    ]
    other_models = [
        m for m in hmm_models
        if m != model
    ]
    if len(other_models) == 0:
        relative_score = pd.Series(
            np.nan,
            index=comparison_matrix.index
        )
    else:
        best_other_score = (       # Highest score from any competing model
            comparison_matrix[
                other_models
            ].max(axis=1)
        )
        relative_score = (
            model_score
            - best_other_score
        )

        relative_score = (
            relative_score.replace(
                [np.inf, -np.inf],
                np.nan
            )
        )

    relative_df = (
        relative_score
        .rename(relative_col)
        .reset_index()
    )

    final_df = final_df.merge(relative_df, on=["sequence", "true_class"], how="left")

accuracy = (classification_df["correct_classification"].mean())

print()
print("========================================")
print("HMM VALIDATION")
print("========================================")

print(
    f"Total test sequences: "
    f"{len(classification_df)}"
)

print(
    f"Correctly classified: "
    f"{classification_df['correct_classification'].sum()}"
)

print(
    f"Accuracy: "
    f"{accuracy * 100:.2f}%"
)

print()


def calculate_relative_cutoff(
    dataframe,
    model
):

    relative_col = (
        f"{model}_relative_score"
    )


    if relative_col not in dataframe.columns:

        return {

            "lowest_true_relative_score":
                np.nan,

            "highest_false_relative_score":
                np.nan,

            "relative_score_separation":
                np.nan,

            "relative_score_cutoff":
                np.nan
        }

    true_relative_scores = dataframe.loc[
        dataframe["true_class"] == model,
        relative_col
    ].dropna()

    false_relative_scores = dataframe.loc[
        dataframe["true_class"] != model,
        relative_col
    ].dropna()


    if len(true_relative_scores) > 0:
        lowest_true_relative_score = (
            true_relative_scores.min()
        )
    else:
        lowest_true_relative_score = np.nan
    if len(false_relative_scores) > 0:
        highest_false_relative_score = (
            false_relative_scores.max()
        )
    else:
        highest_false_relative_score = np.nan

    if (
        pd.notna(lowest_true_relative_score)
        and
        pd.notna(highest_false_relative_score)
    ):
        relative_score_separation = (
            lowest_true_relative_score
            -
            highest_false_relative_score
        )
    else:
        relative_score_separation = np.nan

    if (
        pd.notna(relative_score_separation)
        and
        relative_score_separation > 0
    ):
        relative_score_cutoff = (
            lowest_true_relative_score
            +
            highest_false_relative_score
        ) / 2
    else:
        relative_score_cutoff = np.nan
    return {
        "lowest_true_relative_score":lowest_true_relative_score,
        "highest_false_relative_score":highest_false_relative_score,
        "relative_score_separation":relative_score_separation,
        "relative_score_cutoff":relative_score_cutoff
    }

models = sorted(
    set(
        classification_df[
            "true_class"
        ].dropna()
    )
    |
    set(
        classification_df[
            "best_model"
        ].dropna()
    )
)
stats = []
for model in models:
    actual_positive = (classification_df["true_class"] == model)
    predicted_positive = (classification_df["best_model"] == model)
    TP = (actual_positive&predicted_positive).sum()
    FP = (~actual_positive&predicted_positive).sum()
    FN = (actual_positive&~predicted_positive).sum()
    TN = (~actual_positive&~predicted_positive).sum()
    precision = (TP / (TP + FP)
        if (TP + FP) > 0
        else np.nan)
    recall = (TP / (TP + FN)
        if (TP + FN) > 0
        else np.nan)
    false_positive_rate = (FP / (FP + TN)
        if (FP + TN) > 0
        else np.nan)
    false_negative_rate = (FN / (FN + TP)
        if (FN + TP) > 0
        else np.nan)

    if model in final_df.columns:
        true_scores = final_df.loc[
            final_df["true_class"] == model,
            model
        ].dropna()
        false_scores = final_df.loc[
            final_df["true_class"] != model,
            model
        ].dropna()
        if len(true_scores) > 0:
            lowest_true_score = (
                true_scores.min()
            )
        else:
            lowest_true_score = np.nan
        if len(false_scores) > 0:
            highest_false_score = (
                false_scores.max()
            )
        else:
            highest_false_score = np.nan
    else:
        lowest_true_score = np.nan
        highest_false_score = np.nan

    relative_info = (calculate_relative_cutoff(final_df, model))

    stats.append({
        "model":model,
        "TP":TP,
        "FP":FP,
        "TN":TN,
        "FN":FN,
        "precision":precision,
        "recall":recall,
        "false_positive_rate":false_positive_rate,
        "false_negative_rate":false_negative_rate,
        "lowest_true_score":lowest_true_score,
        "highest_false_score":highest_false_score,
        "lowest_true_relative_score":relative_info["lowest_true_relative_score"],
        "highest_false_relative_score":relative_info["highest_false_relative_score"],
        "relative_score_separation":relative_info["relative_score_separation"],
        "relative_score_cutoff":relative_info["relative_score_cutoff"]
    })

stats_df = pd.DataFrame(stats)
basic_columns = [
    "sequence",
    "true_class",
    "best_model",
    "best_score",
    "second_model",
    "second_score",
    "delta_score",
    "correct_classification"
]

model_columns = []
for model in hmm_models:
    if model in final_df.columns:
        model_columns.append(
            model
        )
    relative_col = (
        f"{model}_relative_score"
    )
    if relative_col in final_df.columns:
        model_columns.append(
            relative_col
        )
final_df = final_df[basic_columns+model_columns]

print()
print("========================================")
print("MODEL STATISTICS")
print("========================================")

display_columns = [
    "model",
    "TP",
    "FP",
    "precision",
    "recall",
    "lowest_true_score",
    "highest_false_score",
    "lowest_true_relative_score",
    "highest_false_relative_score",
    "relative_score_separation",
    "relative_score_cutoff"
]

print(stats_df[display_columns].to_string(index=False))

with pd.ExcelWriter(
    output_file
) as writer:
    final_df.to_excel(
        writer,
        sheet_name="score_matrix",
        index=False
    )
    stats_df.to_excel(
        writer,
        sheet_name="model_statistics",
        index=False
    )
    df.to_excel(
        writer,
        sheet_name="raw_scores",
        index=False
    )
print()
print(
    f"Results saved to: "
    f"{output_file}"
)
