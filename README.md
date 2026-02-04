# Age & Ethnicity Disease Scoring System

## Quick Reference: Data Sources & Scoring Logic

- **Data Sources:** Monarch Initiative (disease/phenotype/inheritance), Orphanet (age of onset, rare disease), gnomAD (population allele frequencies)
- **Key Functions:** See `scoring_engine.py`, `age_demographic_modifier.py`, and `monarch_client.py` for how data is fetched and used.
- **Scoring Logic:**
  - Raw HPO score (patient vs. disease phenotypes)
  - Age, demographic, and inheritance modifiers (from external data)
  - **Final Score = Raw HPO Score × Age Modifier × Demographic Modifier × Inheritance Modifier**

For a detailed breakdown, see [SCORING_LOGIC.md](SCORING_LOGIC.md).

This project provides a scalable, API-driven system for scoring the probability of genetic diseases based on patient phenotypes (HPO terms), age, ethnicity, and inheritance. It integrates live data from Monarch Initiative, gnomAD, and Orphanet, with no hardcoded disease information.

---

## Project Structure

```
scoring_system_live/
├── age_demographic_modifier.py   # Main modifier logic (age, ethnicity, inheritance)
├── scoring_engine.py             # HPO→Disease scoring engine (raw scores)
├── monarch_client.py             # Monarch API client (legacy, used by scoring_engine)
├── test_modifier.py              # Demo/test script for modifier
├── demo.py                       # Example usage/demo
├── cache/                        # Local cache for API responses
│   └── *.json                    # Cached Monarch/gnomAD/Orphanet results
└── README.md                     # Project documentation
```

---

## Scoring Pipeline

### 1. HPO → Disease Scoring (`scoring_engine.py`)
- **Input**: Patient HPO terms, age, population, sex
- **Process**:
  - For each disease, fetch HPO phenotype profile from Monarch
  - Compare patient HPOs to disease HPOs (semantic similarity, Jaccard, etc.)
  - Calculate a raw score (probability or similarity) for each disease
- **Output**: Raw HPO→Disease score for each candidate disease

### 2. Modifier Calculation (`age_demographic_modifier.py`)
- **Input**: Raw HPO→Disease score, patient age, population, sex, disease ID
- **Process**:
  - **Age Modifier**: Uses Monarch and Orphanet to determine how well patient age matches disease onset (e.g., perfect, adjacent, poor)
  - **Demographic Modifier**: Uses gnomAD population allele frequencies for causal genes to enrich/deplete score based on ethnicity
  - **Inheritance Modifier**: Uses Monarch inheritance mode and patient sex to adjust score (e.g., X-linked, autosomal)
- **Output**: Final score = Raw HPO Score × Age Modifier × Demographic Modifier × Inheritance Modifier

### 3. API Integration
- **Monarch Initiative**: Disease phenotypes, inheritance, age of onset
- **gnomAD**: Population allele frequencies for causal genes
- **Orphanet**: Age of onset and natural history for rare diseases

---

## Main Modules

### `scoring_engine.py`
- Implements the initial HPO→Disease scoring logic
- Handles Monarch API queries for disease phenotypes, inheritance, age of onset
- Provides `score_disease()` and `search_and_score()` functions

### `age_demographic_modifier.py`
- Provides the `apply_modifier()` function to adjust raw scores
- Integrates Monarch, gnomAD, and Orphanet APIs
- Contains logic for age, demographic, and inheritance modifiers
- Returns detailed reasoning and raw data for transparency

### `test_modifier.py`
- Demonstrates usage of the modifier system
- Shows how OMIM→MONDO mapping, Orphanet fallback, and gnomAD enrichment work

---

## Example Usage

```python
from scoring_engine import PatientProfile, score
from age_demographic_modifier import apply_modifier

# Step 1: Raw HPO→Disease scoring
patient = PatientProfile(age=5, population='asj', sex='male')
raw_result = score(patient, 'OMIM:272800')  # Tay-Sachs

# Step 2: Apply age/ethnicity/inheritance modifier
hpo_score = 0.75  # Example raw score
final_score, details = apply_modifier(
    hpo_score=hpo_score,
    disease_id='OMIM:272800',
    patient_age=5,
    patient_population='asj',
    patient_sex='male'
)
print(f"Final score: {final_score:.3f}")
print(details)
```

---

## Extending/Customizing

- Plug in your own HPO→Disease scoring function (semantic similarity, ML, etc.)
- Extend modifier logic for other demographic or clinical factors
- Add new API sources or caching strategies

---

## Running Tests

```bash
python3 test_modifier.py
```

---

## License

MIT

---

## Contact

For questions or contributions, open an issue or contact the maintainer.
