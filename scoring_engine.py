"""
Disease Probability Scoring Engine
===================================

Uses Monarch Initiative API for on-demand disease data.
No bulk downloads - fetches only what's needed.

Scoring Formula:
    Adjusted Prevalence = Base Prevalence × Age Factor × Ethnicity Factor × Inheritance Factor
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from monarch_client import MonarchClient, get_client


# =============================================================================
# CONFIGURATION
# =============================================================================

# gnomAD population codes
POPULATIONS = {
    "afr": "African/African American",
    "amr": "Latino/Admixed American", 
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "Finnish",
    "mid": "Middle Eastern",
    "nfe": "European (non-Finnish)",
    "sas": "South Asian",
}

# HPO age of onset mapping
HPO_AGE_OF_ONSET = {
    "HP:0030674": {"name": "Antenatal", "min": -0.75, "max": 0},
    "HP:0003577": {"name": "Congenital", "min": 0, "max": 0},
    "HP:0003623": {"name": "Neonatal", "min": 0, "max": 0.08},
    "HP:0003593": {"name": "Infantile", "min": 0.08, "max": 1},
    "HP:0011463": {"name": "Childhood", "min": 1, "max": 5},
    "HP:0003621": {"name": "Juvenile", "min": 5, "max": 15},
    "HP:0011462": {"name": "Young adult", "min": 16, "max": 40},
    "HP:0003581": {"name": "Adult", "min": 16, "max": 60},
    "HP:0003584": {"name": "Late onset", "min": 40, "max": 100},
}

# Known founder effect diseases (ethnicity risk multipliers)
FOUNDER_EFFECTS = {
    # (disease_keyword, population): risk_multiplier
    ("tay-sachs", "asj"): 30.0,
    ("gaucher", "asj"): 15.0,
    ("canavan", "asj"): 20.0,
    ("niemann-pick", "asj"): 10.0,
    ("familial dysautonomia", "asj"): 15.0,
    ("cystic fibrosis", "nfe"): 2.0,
    ("cystic fibrosis", "asj"): 2.0,
    ("sickle cell", "afr"): 20.0,
    ("mediterranean fever", "mid"): 20.0,
    ("mediterranean fever", "asj"): 15.0,
    ("beta-thalassemia", "mid"): 10.0,
    ("beta-thalassemia", "sas"): 5.0,
    ("huntington", "nfe"): 2.0,
}

# Default prevalence for rare diseases (per 100,000)
DEFAULT_RARE_PREVALENCE = 1.0


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class PatientProfile:
    """Patient demographics."""
    age: int
    population: str  # gnomAD code: afr, amr, asj, eas, fin, mid, nfe, sas
    sex: str = "unknown"  # male, female, unknown
    
    def __post_init__(self):
        if self.population not in POPULATIONS:
            raise ValueError(f"Unknown population '{self.population}'. Use: {list(POPULATIONS.keys())}")


@dataclass
class DiseaseScore:
    """Scoring result."""
    disease_id: str
    name: str
    base_prevalence: float
    adjusted_prevalence: float
    age_factor: float
    ethnicity_factor: float
    inheritance_factor: float
    risk_category: str
    probability: str  # e.g., "1 in 50,000"
    explanation: Dict[str, str] = field(default_factory=dict)


# =============================================================================
# SCORING ENGINE
# =============================================================================

class ScoringEngine:
    """
    On-demand disease probability scoring using Monarch API.
    """
    
    def __init__(self):
        self.client = get_client()
    
    def _get_age_category(self, age: int) -> Tuple[str, str]:
        """Map patient age to HPO category."""
        for hpo_code, info in HPO_AGE_OF_ONSET.items():
            if info["min"] <= age <= info["max"]:
                return hpo_code, info["name"]
        return "HP:0003581", "Adult"
    
    def _calculate_age_factor(self, patient_age: int, disease_onset_codes: List[str]) -> float:
        """
        Calculate age match factor.
        
        Returns:
            1.0 = perfect match
            0.3 = adjacent category
            0.1 = two categories away
            0.05 = no match
        """
        if not disease_onset_codes:
            return 0.5  # Unknown
        
        patient_hpo, _ = self._get_age_category(patient_age)
        
        if patient_hpo in disease_onset_codes:
            return 1.0  # Direct match
        
        # Calculate distance
        ordered = list(HPO_AGE_OF_ONSET.keys())
        try:
            patient_idx = ordered.index(patient_hpo)
        except ValueError:
            return 0.1
        
        min_distance = float('inf')
        for onset_code in disease_onset_codes:
            try:
                onset_idx = ordered.index(onset_code)
                min_distance = min(min_distance, abs(patient_idx - onset_idx))
            except ValueError:
                continue
        
        if min_distance == 1:
            return 0.3
        elif min_distance == 2:
            return 0.1
        return 0.05
    
    def _calculate_ethnicity_factor(self, population: str, disease_name: str) -> float:
        """
        Calculate ethnicity risk factor based on known founder effects.
        """
        disease_lower = disease_name.lower()
        
        for (keyword, pop), factor in FOUNDER_EFFECTS.items():
            if keyword in disease_lower and population == pop:
                return factor
        
        return 1.0
    
    def _calculate_inheritance_factor(self, inheritance: List[str], sex: str) -> float:
        """
        Calculate sex-based inheritance factor.
        
        X-linked recessive: males much higher risk
        """
        if 'x_linked_recessive' in inheritance:
            if sex == 'male':
                return 2.0
            elif sex == 'female':
                return 0.1
        
        if 'x_linked_dominant' in inheritance:
            if sex == 'female':
                return 1.5
            elif sex == 'male':
                return 0.5
        
        return 1.0
    
    def _categorize_risk(self, prevalence: float) -> str:
        """Categorize adjusted prevalence."""
        if prevalence >= 50:
            return "high"
        elif prevalence >= 10:
            return "elevated"
        elif prevalence >= 1:
            return "moderate"
        elif prevalence >= 0.1:
            return "low"
        return "very_low"
    
    def score_disease(self, patient: PatientProfile, disease_id: str) -> Optional[DiseaseScore]:
        """
        Score a single disease for a patient.
        
        Args:
            patient: Patient demographics
            disease_id: OMIM ID (e.g., "219700" or "OMIM:219700")
            
        Returns:
            DiseaseScore or None if disease not found
        """
        # Fetch disease from Monarch API
        disease = self.client.get_disease(disease_id)
        if not disease:
            return None
        
        name = disease.get('name', f'Disease {disease_id}')
        
        # Get age of onset
        onset_codes = self.client.get_age_of_onset(disease_id)
        
        # Get inheritance
        inheritance = self.client.get_inheritance(disease_id)
        
        # Calculate factors
        age_factor = self._calculate_age_factor(patient.age, onset_codes)
        ethnicity_factor = self._calculate_ethnicity_factor(patient.population, name)
        inheritance_factor = self._calculate_inheritance_factor(inheritance, patient.sex)
        
        # Base prevalence (Monarch doesn't provide this directly, use default for rare)
        base_prevalence = DEFAULT_RARE_PREVALENCE
        
        # Adjusted prevalence
        adjusted = base_prevalence * age_factor * ethnicity_factor * inheritance_factor
        
        # Risk category
        risk = self._categorize_risk(adjusted)
        
        # Probability string
        if adjusted > 0:
            odds = int(100000 / adjusted)
            probability = f"1 in {odds:,}"
        else:
            probability = "N/A"
        
        # Build explanation
        patient_age_cat = self._get_age_category(patient.age)[1]
        onset_names = [HPO_AGE_OF_ONSET.get(c, {}).get('name', c) for c in onset_codes]
        
        explanation = {
            "disease_name": name,
            "patient": f"{patient.age}yo {patient.sex} ({POPULATIONS.get(patient.population)})",
            "age_match": f"{age_factor:.0%} (patient: {patient_age_cat}, disease: {', '.join(onset_names) or 'unknown'})",
            "ethnicity_factor": f"{ethnicity_factor:.1f}x",
            "inheritance": f"{', '.join(inheritance) or 'unknown'} (factor: {inheritance_factor:.1f}x)",
            "base_prevalence": f"{base_prevalence}/100k (default for rare disease)",
            "adjusted_prevalence": f"{adjusted:.4f}/100k",
        }
        
        return DiseaseScore(
            disease_id=disease_id,
            name=name,
            base_prevalence=base_prevalence,
            adjusted_prevalence=adjusted,
            age_factor=age_factor,
            ethnicity_factor=ethnicity_factor,
            inheritance_factor=inheritance_factor,
            risk_category=risk,
            probability=probability,
            explanation=explanation
        )
    
    def search_and_score(
        self, 
        patient: PatientProfile, 
        query: str, 
        limit: int = 5
    ) -> List[DiseaseScore]:
        """
        Search for diseases and score them.
        
        Args:
            patient: Patient demographics
            query: Search string (e.g., "cystic fibrosis")
            limit: Max results
            
        Returns:
            List of scored diseases
        """
        # Search via Monarch API
        results = self.client.search_diseases(query, limit=limit * 2)
        
        scores = []
        for result in results:
            disease_id = result.get('id', '')
            if not disease_id:
                continue
            
            score = self.score_disease(patient, disease_id)
            if score:
                scores.append(score)
            
            if len(scores) >= limit:
                break
        
        # Sort by adjusted prevalence
        scores.sort(key=lambda x: x.adjusted_prevalence, reverse=True)
        return scores


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def score(patient: PatientProfile, disease_id: str) -> Dict:
    """
    Quick function to score a disease.
    
    Example:
        patient = PatientProfile(age=5, population="asj", sex="male")
        result = score(patient, "272800")  # Tay-Sachs
    """
    engine = ScoringEngine()
    result = engine.score_disease(patient, disease_id)
    
    if not result:
        return {"error": f"Disease {disease_id} not found"}
    
    return {
        "disease": result.name,
        "id": result.disease_id,
        "probability": result.probability,
        "risk": result.risk_category,
        "factors": {
            "age": f"{result.age_factor:.0%}",
            "ethnicity": f"{result.ethnicity_factor:.1f}x",
            "inheritance": f"{result.inheritance_factor:.1f}x",
        },
        "explanation": result.explanation
    }


def search(patient: PatientProfile, query: str) -> List[Dict]:
    """
    Search and score diseases.
    
    Example:
        patient = PatientProfile(age=35, population="nfe", sex="female")
        results = search(patient, "huntington")
    """
    engine = ScoringEngine()
    results = engine.search_and_score(patient, query)
    
    return [
        {
            "disease": r.name,
            "id": r.disease_id,
            "probability": r.probability,
            "risk": r.risk_category,
        }
        for r in results
    ]
