import json
import argparse
import pandas as pd
from pathlib import Path


def get_best_params(lisi_log_path):
    lisi_log_path = Path(lisi_log_path)
    output_dir = Path(lisi_log_path).parent

    if not lisi_log_path.exists():
        raise FileNotFoundError(f"No LISI log found at {lisi_log_path}")

    results_df = pd.read_csv(lisi_log_path)
    best_row = results_df.loc[results_df["lisi_score"].idxmax()]
    best_params = best_row.drop("lisi_score").to_dict()

    best_params_path = output_dir / "best_params.json"
    with open(best_params_path, "w") as f:
        json.dump(best_params, f, indent=2)

    print(f"Best params saved to {best_params_path}")
    print(best_params)

    return best_params


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract best hyperparameters from LISI log")
    parser.add_argument("--log", type=str, required=True, help="Path to LISI log CSV")
    args = parser.parse_args()

    get_best_params(args.log)