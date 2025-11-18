import numpy as np
import pandas as pd
import os, sys, glob
import yaml
import argparse
import subprocess

# Command : ./apps/EmpiricalPseudoPotentialMain -m Si -b 12 -N 500 -n 10 -r output_dir -j 1 -p LGXK -P -d ../parameter_files/materials-kim-fischetti.yaml

BIN_EPM = "./apps/EmpiricalPseudoPotentialMain"

def run_epm(config_path):
    cmd = [
        BIN_EPM,
        "-m", "Si",
        "-b", "12",
        "-N", "500",
        "-n", "10",
        "-r", "output_dir",
        "-j", "1",
        "-p", "LGXK",
        "-P",
        "-d", config_path
    ]
    subprocess.run(cmd)
    print("EPM run completed.")
    return

def get_results(output_dir):
    data_files = glob.glob(os.path.join(output_dir, "EPM_*.csv"))
    all_data = []
    for file in data_files:
        df = pd.read_csv(file)
        all_data.append(df)
    combined_data = pd.concat(all_data, ignore_index=True)
    return combined_data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some CSV files.")
    parser.add_argument('--config', type=str, default="config.yaml", help="Path to the config file")
    args = parser.parse_args()
    print("Loading config from:", args.config)
    # Copy config file to current directory
    os.system(f"cp {args.config} ./tmp_config.yaml")
    exit(0)