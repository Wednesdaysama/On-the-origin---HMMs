import os
import glob
import numpy as np
import pandas as pd

# true_class__vs__hmm_model.tbl
input_dir = "./hmm_results"
output_file = "./hmm_cutoff_analysis.xlsx"
records = []
tbl_files = glob.glob(
    os.path.join(input_dir, "*.tbl")
)
print(f"Found {len(tbl_files)} tbl files.")
for tbl_file in tbl_files:
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
            if len(parts) < 6:
                continue
            sequence = parts[0]
            evalue = float(parts[4])
            score = float(parts[5])
            records.append(
                {
                    "sequence": sequence,
                    "true_class": true_class,
                    "hmm_model": hmm_model,
                    "score": score,
                    "evalue": evalue
                }
            )
df = pd.DataFrame(records)
print(f"Total HMM hits read: {len(df)}")
if df.empty:
    raise ValueError(f"No HMM hits were found in {input_dir}")

score_matrix = df.pivot_table(
    index=[
        "sequence",
        "true_class"
    ],
    columns="hmm_model",
    values="score",
    aggfunc="max"
)
comparison_matrix = score_matrix.fillna(
    float("-inf")
)

hmm_models = list(
    score_matrix.columns
)

print(f"Number of HMM models: {len(hmm_models)}")

classification_results = []
for (
    sequence,
    true_class
), row in comparison_matrix.iterrows():
    sorted_scores = row.sort_values(
        ascending=False
    )
    best_model = sorted_scores.index[0]
    best_score = sorted_scores.iloc[0]
    if len(sorted_scores) > 1:
        second_model = sorted_scores.index[1]
        second_score = sorted_scores.iloc[1]
    else:
        second_model = None
        second_score = np.nan
    if np.isneginf(best_score):
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
        else:
            delta_score = (
                best_score
                -
                second_score
            )
        correct = (
            best_model == true_class
        )
    classification_results.append(
        {
            "sequence": sequence,
            "true_class": true_class,
            "best_model": best_model,
            "best_score": best_score,
            "second_model": second_model,
            "second_score": second_score,
            "delta_score": delta_score,
            "correct_classification": correct
        }
    )
classification_df = pd.DataFrame(
    classification_results
)

score_output = score_matrix.reset_index()
final_df = classification_df.merge(
    score_output,
    on=[
        "sequence",
        "true_class"
    ],
    how="left"
)

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
        best_other_score = (
            comparison_matrix[
                other_models
            ]
            .max(axis=1)
        )
        relative_score = (model_score - best_other_score)
        target_no_hit = np.isneginf(
            model_score
        )
        relative_score[
            target_no_hit
        ] = np.nan
        relative_score = (
            relative_score.replace(
                [
                    np.inf,
                    -np.inf
                ],
                np.nan
            )
        )

    relative_df = (
        relative_score
        .rename(relative_col)
        .reset_index()
    )
    final_df = final_df.merge(
        relative_df,
        on=[
            "sequence",
            "true_class"
        ],
        how="left"
    )

def calculate_metrics(
    actual_positive,
    predicted_positive
):
    actual_positive = np.asarray(
        actual_positive,
        dtype=bool
    )
    predicted_positive = np.asarray(
        predicted_positive,
        dtype=bool
    )
    TP = np.sum(
        actual_positive
        &
        predicted_positive
    )
    FP = np.sum(
        (~actual_positive)
        &
        predicted_positive
    )
    FN = np.sum(
        actual_positive
        &
        (~predicted_positive)
    )
    TN = np.sum(
        (~actual_positive)
        &
        (~predicted_positive)
    )
    TPR = (
        TP / (TP + FN)
        if (TP + FN) > 0
        else np.nan
    )
    FPR = (
        FP / (FP + TN)
        if (FP + TN) > 0
        else np.nan
    )
    FNR = (
        FN / (FN + TP)
        if (FN + TP) > 0
        else np.nan
    )
    precision = (
        TP / (TP + FP)
        if (TP + FP) > 0
        else np.nan
    )
    specificity = (
        TN / (TN + FP)
        if (TN + FP) > 0
        else np.nan
    )
    return {
        "TP": int(TP),
        "FP": int(FP),
        "TN": int(TN),
        "FN": int(FN),
        "TPR": TPR,
        "FPR": FPR,
        "FNR": FNR,
        "precision": precision,
        "specificity": specificity
    }

def scan_cutoffs(
    dataframe,
    model,
    score_column,
    score_type
):
    actual_positive = (
        dataframe["true_class"] == model
    )
    raw_scores = dataframe[
        score_column
    ]
    scores = raw_scores.fillna(
        float("-inf")
    )
    finite_scores = raw_scores[
        np.isfinite(raw_scores)
    ]
    unique_scores = np.sort(
        finite_scores.unique()
    )
    if len(unique_scores) == 0:
        return pd.DataFrame()
    cutoff_results = []
    for cutoff in unique_scores:
        predicted_positive = (
            scores >= cutoff
        )
        metrics = calculate_metrics(
            actual_positive,
            predicted_positive
        )
        if (
            pd.notna(metrics["FPR"])
            and
            pd.notna(metrics["FNR"])
        ):
            balance_difference = abs(
                metrics["FPR"] - metrics["FNR"]
            )
        else:
            balance_difference = np.nan
        if (
            pd.notna(metrics["TPR"])
            and
            pd.notna(metrics["FPR"])
        ):
            youden_J = (
                metrics["TPR"] - metrics["FPR"]
            )
        else:
            youden_J = np.nan
        cutoff_results.append(
            {
                "model": model,
                "score_type": score_type,
                "cutoff": cutoff,
                "TP": metrics["TP"],
                "FP": metrics["FP"],
                "TN": metrics["TN"],
                "FN": metrics["FN"],
                "TPR": metrics["TPR"],
                "FPR": metrics["FPR"],
                "FNR": metrics["FNR"],
                "precision": metrics["precision"],
                "specificity": metrics["specificity"],
                "abs_FPR_minus_FNR":
                    balance_difference,
                "youden_J":
                    youden_J
            }
        )
    return pd.DataFrame(
        cutoff_results
    )
absolute_scan_list = []
print()
print(
    "Scanning absolute-score cutoffs..."
)

for model in hmm_models:
    if model not in final_df.columns:
        continue
    model_scan = scan_cutoffs(
        dataframe=final_df,
        model=model,
        score_column=model,
        score_type="absolute"
    )
    if not model_scan.empty:
        absolute_scan_list.append(
            model_scan
        )
if absolute_scan_list:
    absolute_cutoff_scan = pd.concat(
        absolute_scan_list,
        ignore_index=True
    )
else:
    absolute_cutoff_scan = pd.DataFrame()
relative_scan_list = []
print(
    "Scanning relative-score cutoffs..."
)

for model in hmm_models:
    relative_col = (
        f"{model}_relative_score"
    )
    if relative_col not in final_df.columns:
        continue
    model_scan = scan_cutoffs(
        dataframe=final_df,
        model=model,
        score_column=relative_col,
        score_type="relative"
    )
    if not model_scan.empty:
        relative_scan_list.append(
            model_scan
        )
if relative_scan_list:
    relative_cutoff_scan = pd.concat(
        relative_scan_list,
        ignore_index=True
    )
else:
    relative_cutoff_scan = pd.DataFrame()
def make_summary(
    cutoff_df
):
    if cutoff_df.empty:
        return pd.DataFrame()
    summary_rows = []
    for model, group in cutoff_df.groupby(
        "model"
    ):
        group = group.copy()
        group["total_error_rate"] = (
            group["FPR"]
            +
            group["FNR"]
        )
        min_balance = (
            group[
                "abs_FPR_minus_FNR"
            ].min()
        )
        candidates = group[
            group[
                "abs_FPR_minus_FNR"
            ] == min_balance
        ].copy()
        min_total_error = (
            candidates[
                "total_error_rate"
            ].min()
        )
        candidates = candidates[
            candidates[
                "total_error_rate"
            ] == min_total_error
        ].copy()
        candidates = candidates.sort_values(
            by=[
                "precision",
                "cutoff"
            ],
            ascending=[
                False,
                False
            ]
        )
        best = candidates.iloc[0]
        N_positive = (
            final_df[
                "true_class"
            ] == model
        ).sum()
        N_negative = (
            final_df[
                "true_class"
            ] != model
        ).sum()
        summary_rows.append(
            {
                "model":
                    model,
                "score_type":
                    best["score_type"],
                "N_positive":
                    int(N_positive),
                "N_negative":
                    int(N_negative),
                "selected_cutoff":
                    best["cutoff"],
                "TP":
                    int(best["TP"]),
                "FP":
                    int(best["FP"]),
                "TN":
                    int(best["TN"]),
                "FN":
                    int(best["FN"]),
                "TPR":
                    best["TPR"],
                "FPR":
                    best["FPR"],
                "FNR":
                    best["FNR"],
                "precision":
                    best["precision"],
                "specificity":
                    best["specificity"],
                "abs_FPR_minus_FNR":
                    best[
                        "abs_FPR_minus_FNR"
                    ],
                "FPR_plus_FNR":
                    best[
                        "total_error_rate"
                    ],
                "youden_J":
                    best["youden_J"]
            }
        )

    return pd.DataFrame(
        summary_rows
    )
absolute_summary = make_summary(
    absolute_cutoff_scan
)
relative_summary = make_summary(
    relative_cutoff_scan
)


best_model_stats = []
for model in hmm_models:
    actual_positive = (
        classification_df[
            "true_class"
        ] == model
    )
    predicted_positive = (
        classification_df[
            "best_model"
        ] == model
    )
    metrics = calculate_metrics(
        actual_positive,
        predicted_positive
    )
    best_model_stats.append(
        {
            "model":
                model,
            "TP":
                metrics["TP"],
            "FP":
                metrics["FP"],
            "TN":
                metrics["TN"],
            "FN":
                metrics["FN"],
            "TPR":
                metrics["TPR"],
            "FPR":
                metrics["FPR"],
            "FNR":
                metrics["FNR"],
            "precision":
                metrics["precision"],
            "specificity":
                metrics["specificity"]
        }
    )

best_model_stats_df = pd.DataFrame(
    best_model_stats
)

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
final_df = final_df[
    basic_columns
    +
    model_columns
]
accuracy = (
    classification_df[
        "correct_classification"
    ].mean()

)

if not absolute_summary.empty:
    print(
        absolute_summary[
            [
                "model",
                "N_positive",
                "N_negative",
                "selected_cutoff",
                "TPR",
                "FPR",
                "FNR",
                "precision"
            ]
        ].to_string(
            index=False
        )
    )

if not relative_summary.empty:
    print(
        relative_summary[
            [
                "model",
                "N_positive",
                "N_negative",
                "selected_cutoff",
                "TPR",
                "FPR",
                "FNR",
                "precision"
            ]
        ].to_string(
            index=False
        )

    )
with pd.ExcelWriter(
    output_file,
    engine="openpyxl"
) as writer:
    final_df.to_excel(
        writer,
        sheet_name="score_matrix",
        index=False
    )
    absolute_cutoff_scan.to_excel(
        writer,
        sheet_name="absolute_cutoff_scan",
        index=False
    )
    relative_cutoff_scan.to_excel(
        writer,
        sheet_name="relative_cutoff_scan",
        index=False
    )
    absolute_summary.to_excel(
        writer,
        sheet_name="absolute_summary",
        index=False
    )
    relative_summary.to_excel(
        writer,
        sheet_name="relative_summary",
        index=False
    )
    best_model_stats_df.to_excel(
        writer,
        sheet_name="best_model_stats",
        index=False
    )
    df.to_excel(
        writer,
        sheet_name="raw_scores",
        index=False
    )
print(
    f"Results saved to: {output_file}"
)
