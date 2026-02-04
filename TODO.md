# TODO: Data Usage and Scoring Improvements

## 1. Enable Offline Data Usage
- Download Monarch, Orphanet, and gnomAD data dumps.
- Refactor code to load and query local files instead of making API calls.
- Ensure all data dependencies are satisfied for full offline scoring.

## 2. Improve Scoring Methodology
- Move from rule-based multipliers to data-driven calibration.
- Integrate machine learning or statistical models for risk estimation.
- Add uncertainty quantification and validate on real patient data.
