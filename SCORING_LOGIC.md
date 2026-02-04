# Scoring System: Data Sources, Functions, and Logic

## Data Sources

- **Monarch Initiative**: Provides disease phenotypes (HPO terms), age of onset, and inheritance mode.
- **Orphanet**: Supplies additional age of onset and rare disease data (used as fallback).
- **gnomAD**: Delivers population allele frequencies for disease genes, enabling ethnicity-based risk adjustment.

## Key Functions

### Monarch Data
- `get_disease(disease_id)`: Fetches disease info (name, phenotypes, inheritance).
- `get_age_of_onset(disease_id)`: Extracts HPO codes for age of onset.
- `get_inheritance(disease_id)`: Extracts inheritance mode.

### gnomAD Data
- `get_gene_variant_summary(gene_symbol)`: Gets carrier rates for a gene in each population.
- `get_population_enrichment(gene_symbol, target_population)`: Calculates enrichment ratio for a population.

### Orphanet Data
- `get_age_of_onset(disease_id)`: Fallback for age of onset if Monarch lacks data.

### Scoring Functions
- `score_disease(patient, disease_id)`: Main scoring function (raw + modifiers).
- `search_and_score(patient, query)`: Searches and scores multiple diseases.
- `AgeDemographicModifier.calculate(...)`: Advanced scoring with all modifiers and live data.

## Scoring Logic

1. **Raw HPO Score**: Measures similarity between patient HPO terms and disease HPO profile.
2. **Modifiers**:
   - **Age Modifier**: Multiplies score based on how well patient age matches disease onset (Monarch/Orphanet).
   - **Demographic Modifier**: Multiplies score based on population enrichment (gnomAD or founder effect table).
   - **Inheritance Modifier**: Multiplies score based on inheritance mode and patient sex (Monarch).

**Formula:**

    Final Score = Raw HPO Score × Age Modifier × Demographic Modifier × Inheritance Modifier

- If data is missing, modifiers default to 1.0 (neutral).
- All data is fetched live and cached for efficiency.

---

For more details, see the code in `scoring_engine.py`, `age_demographic_modifier.py`, and `monarch_client.py`.
