#!/usr/bin/env python

import pandas as pd
import os

# Define Purkinje range
PURKINJE_RANGE = range(547680, 548790)

def analyze_activation_after(time_threshold, myocardial_df):
    """
    Analyze activation times after a given threshold for myocardial nodes.
    
    Args:
        time_threshold (float): Time after which to analyze activations.
    
    Returns:
        dict: Contains earliest and latest activations and total span.
    """
    post_time_df = myocardial_df[myocardial_df["activation_time"] > time_threshold]

    if post_time_df.empty:
        return {"message": "No myocardial activations found after the given time."}

    earliest_row = post_time_df.loc[post_time_df["activation_time"].idxmin()]
    latest_row = post_time_df.loc[post_time_df["activation_time"].idxmax()]

    total_activation_time = latest_row["activation_time"] - earliest_row["activation_time"]

    return {
        "earliest_activation_time": earliest_row["activation_time"],
        "earliest_node": int(earliest_row["node"]),
        "latest_activation_time": latest_row["activation_time"],
        "latest_node": int(latest_row["node"]),
        "total_activation_time_span": total_activation_time
    }

def run():
    parent_folder = "/media/anthonyowusu-mensah/TonySSD/tryNew/noSTARTstate/lumpMatrix/testSimul2/sinusCV"

    file_names = [
        "Ctrl-His-1000-GM-1.0-S1S2-0-Stim-50-Rpmj-45e3-mGNa-1.0-pGNa-1.0-vtx-noEctopy-savF-1-Mlump-0-parabSol-0-fit-1.0",
        "Ctrl-His-1000-GM-0.3-S1S2-0-Stim-50-Rpmj-45e3-mGNa-0.7-pGNa-0.95-vtx-noEctopy-savF-1-Mlump-1-parabSol-0-fit-1.25",
        "LQT2-His-1000-GM-1.0-S1S2-0-Stim-50-Rpmj-45e3-mGNa-1.0-pGNa-1.0-vtx-noEctopy-savF-1-Mlump-0-parabSol-0-fit-1.0",
        "LQT2-His-1000-GM-0.3-S1S2-0-Stim-50-Rpmj-45e3-mGNa-0.7-pGNa-0.95-vtx-noEctopy-savF-1-Mlump-1-parabSol-0-fit-1.25"
    ]

    file_types = ["Ctrl_Sinus", "Ctrl_Drug", "LQT2_Sinus", "LQT2_Drug"]

    time_threshold = 4000

    for idx, file in enumerate(file_names):
        full_file_path = os.path.join(parent_folder, file, "secondAct-thresh.dat")
        print(full_file_path)

        try:
            df = pd.read_csv(full_file_path, delim_whitespace=True, header=None, names=["node", "activation_time"])
        except FileNotFoundError:
            print(f"File not found: {full_file_path}")
            continue

        # Filter to only myocardial nodes
        myocardial_df = df[~df["node"].isin(PURKINJE_RANGE)]

        print(f"\n--- {file_types[idx]} ---")
        result = analyze_activation_after(time_threshold, myocardial_df)
        for key, value in result.items():
            print(f"{key}: {value}")

# Call the run function
run()

