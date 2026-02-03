#!/usr/bin/env python3
"""Test the age demographic modifier with live API data."""

from age_demographic_modifier import AgeDemographicModifier

def main():
    modifier = AgeDemographicModifier()
    
    # Test 1: Gene extraction
    print("=== Test 1: Gene Extraction ===")
    disease_data = modifier.monarch.get_disease('MONDO:0010100')
    genes = modifier.monarch.get_disease_genes('MONDO:0010100', disease_data)
    print(f"Genes for Tay-Sachs: {genes}")
    
    # Test 2: Inheritance extraction  
    print("\n=== Test 2: Inheritance ===")
    inheritance_mode, raw = modifier._extract_inheritance(disease_data)
    print(f"Inheritance mode: {inheritance_mode}")
    
    # Test 3: gnomAD population data
    print("\n=== Test 3: gnomAD Population Data ===")
    if genes:
        enrichment, pop_data = modifier.gnomad.get_population_enrichment(genes[0], 'asj')
        print(f"ASJ enrichment for HEXA: {enrichment:.2f}")
        pop_rates = pop_data.get("all_populations", {})
        for pop, rate in pop_rates.items():
            print(f"  {pop}: {rate:.6f}")
    
    # Test 4: OMIM to MONDO mapping
    print("\n=== Test 4: OMIM→MONDO Mapping ===")
    mondo_id = modifier.monarch.omim_to_mondo('OMIM:219700')
    print(f"OMIM:219700 (CF) → {mondo_id}")
    
    cf_data = modifier.monarch.get_disease('OMIM:219700')
    if cf_data:
        print(f"Disease name: {cf_data.get('name')}")
        resolved = cf_data.get('_resolved_mondo', 'N/A')
        print(f"Resolved MONDO ID: {resolved}")
        
        # Check inheritance
        inh_mode, _ = modifier._extract_inheritance(cf_data)
        print(f"Inheritance: {inh_mode}")
    
    # Test 5: Orphanet age of onset
    print("\n=== Test 5: Orphanet Age of Onset ===")
    orphanet_data = modifier.orphanet.get_age_of_onset('OMIM:219700')
    if orphanet_data:
        print(f"Orphanet age data: {orphanet_data}")
    else:
        print("No Orphanet data (API may require different endpoint)")
    
    # Test 6: Full modifier with OMIM ID
    print("\n=== Test 6: Full Modifier (CF via OMIM ID) ===")
    result = modifier.calculate(
        disease_id='OMIM:219700',
        patient_age=8,
        patient_population='nfe',
        patient_sex='male'
    )
    print(f"Disease: {result.disease_name}")
    print(f"Age modifier: {result.age_modifier:.2f} ({result.age_match_level})")
    print(f"Age source: {result.raw_age_data.get('source', 'N/A')}")
    print(f"Demo modifier: {result.demographic_modifier:.2f} ({result.demographic_category})")
    print(f"Inheritance: {result.inheritance_mode}")
    print(f"Combined: {result.combined_modifier:.2f}")
    print(f"Data sources: {result.data_sources}")
    
    # Test 7: Apply to HPO score
    print("\n=== Test 7: Apply to HPO Score ===")
    from age_demographic_modifier import apply_modifier
    
    hpo_score = 0.75
    final_score, details = apply_modifier(
        hpo_score=hpo_score,
        disease_id='MONDO:0010100',
        patient_age=1.5,
        patient_population='asj',
        patient_sex='male'
    )
    print(f"Original HPO score: {hpo_score}")
    print(f"After age/demo modifiers: {final_score:.3f}")
    
    # Test 8: Sickle Cell with AFR population (should show enrichment)
    print("\n=== Test 8: Sickle Cell (AFR population) ===")
    result = modifier.calculate(
        disease_id='sickle cell disease',
        patient_age=5,
        patient_population='afr',
        patient_sex='female'
    )
    print(f"Disease: {result.disease_name} ({result.disease_id})")
    print(f"Demo modifier: {result.demographic_modifier:.2f} ({result.demographic_category})")
    print(f"Genes used: {result.raw_demographic_data.get('genes_checked', [])}")

if __name__ == "__main__":
    main()
