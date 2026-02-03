#!/usr/bin/env python3
"""
Disease Probability Scoring - Live Demo
========================================

Uses Monarch Initiative API (no downloads needed).
Data is fetched on-demand and cached locally.
"""

from scoring_engine import ScoringEngine, PatientProfile, POPULATIONS


def print_line(char="=", width=60):
    print(char * width)


def demo_single_disease():
    """Score a specific disease."""
    print_line()
    print(" SCORE A SPECIFIC DISEASE")
    print_line()
    
    engine = ScoringEngine()
    
    # Patient: 2-year-old Ashkenazi Jewish child
    patient = PatientProfile(age=2, population="asj", sex="male")
    
    print(f"Patient: {patient.age} year old {patient.sex}")
    print(f"Population: {POPULATIONS[patient.population]}")
    print()
    
    # Search for Tay-Sachs first (Monarch uses MONDO IDs, not OMIM)
    print("Searching for: Tay-Sachs Disease")
    print("-" * 40)
    
    results = engine.search_and_score(patient, "tay-sachs", limit=1)
    
    if results:
        score = results[0]
        print(f"Disease: {score.name}")
        print(f"ID: {score.disease_id}")
        print(f"Risk Category: {score.risk_category.upper()}")
        print(f"Probability: {score.probability}")
        print()
        print("Factors:")
        print(f"  Age Match: {score.age_factor:.0%}")
        print(f"  Ethnicity: {score.ethnicity_factor:.1f}x (founder effect)")
        print(f"  Inheritance: {score.inheritance_factor:.1f}x")
        print()
        print("Explanation:")
        for key, value in score.explanation.items():
            print(f"  {key}: {value}")
    else:
        print("Disease not found")


def demo_search():
    """Search and score diseases."""
    print()
    print_line()
    print(" SEARCH DISEASES")
    print_line()
    
    engine = ScoringEngine()
    
    # Different patients
    patients = [
        PatientProfile(age=5, population="nfe", sex="male"),
        PatientProfile(age=35, population="asj", sex="female"),
        PatientProfile(age=50, population="afr", sex="male"),
    ]
    
    queries = ["cystic fibrosis", "huntington disease", "sickle cell"]
    
    for query in queries:
        print(f"\nSearching: '{query}'")
        print("-" * 40)
        
        # Use first patient for demo
        patient = patients[0]
        results = engine.search_and_score(patient, query, limit=3)
        
        if not results:
            print("  No results found")
            continue
        
        for i, score in enumerate(results, 1):
            print(f"\n  {i}. {score.name}")
            print(f"     ID: {score.disease_id}")
            print(f"     Risk: {score.risk_category} ({score.probability})")


def demo_population_comparison():
    """Compare same disease across populations."""
    print()
    print_line()
    print(" POPULATION COMPARISON: TAY-SACHS")
    print_line()
    
    engine = ScoringEngine()
    
    # Search for Tay-Sachs to get the MONDO ID
    results = engine.client.search_diseases("tay-sachs", limit=1)
    if not results:
        print("Tay-Sachs not found")
        return
    
    disease_id = results[0].get('id')
    disease_name = results[0].get('name')
    age = 2  # Infantile onset
    
    print(f"Disease: {disease_name}")
    print(f"Disease ID: {disease_id}")
    print(f"Patient Age: {age} years")
    print()
    
    populations = ["asj", "nfe", "afr", "eas"]
    
    for pop in populations:
        patient = PatientProfile(age=age, population=pop, sex="unknown")
        score = engine.score_disease(patient, disease_id)
        
        if score:
            print(f"  {POPULATIONS[pop]}:")
            print(f"    Ethnicity Factor: {score.ethnicity_factor:.1f}x")
            print(f"    Adjusted Prevalence: {score.adjusted_prevalence:.2f}/100k")
            print(f"    Risk: {score.risk_category}")
            print()


def demo_age_impact():
    """Show how age affects scoring."""
    print()
    print_line()
    print(" AGE IMPACT: HUNTINGTON DISEASE")
    print_line()
    
    engine = ScoringEngine()
    
    # Search for Huntington
    results = engine.client.search_diseases("huntington disease", limit=1)
    if not results:
        print("Huntington not found")
        return
    
    disease_id = results[0].get('id')
    print(f"Disease ID: {disease_id}")
    print()
    
    ages = [5, 15, 30, 45, 60, 75]
    
    for age in ages:
        patient = PatientProfile(age=age, population="nfe", sex="unknown")
        score = engine.score_disease(patient, disease_id)
        
        if score:
            age_cat = engine._get_age_category(age)[1]
            print(f"  Age {age} ({age_cat}):")
            print(f"    Age Factor: {score.age_factor:.0%}")
            print(f"    Risk: {score.risk_category}")


def main():
    print()
    print_line("=")
    print(" DISEASE PROBABILITY SCORING SYSTEM")
    print(" Using Monarch Initiative API")
    print_line("=")
    print()
    print("This system fetches data on-demand from:")
    print("  • Monarch Initiative (monarchinitiative.org)")
    print("  • Combines OMIM + Orphanet + HPO data")
    print("  • No bulk downloads - just API calls")
    print()
    
    # Run demos
    demo_single_disease()
    demo_search()
    demo_population_comparison()
    demo_age_impact()
    
    print()
    print_line()
    print(" DEMO COMPLETE")
    print_line()
    print()
    print("How to use in your code:")
    print()
    print("  from scoring_engine import PatientProfile, score, search")
    print()
    print("  # Score a specific disease")
    print("  patient = PatientProfile(age=5, population='asj', sex='male')")
    print("  result = score(patient, 'OMIM:272800')")
    print()
    print("  # Search diseases")
    print("  results = search(patient, 'cystic fibrosis')")
    print()


if __name__ == "__main__":
    main()
