"""
Age & Demographic Modifier for HPO → Disease Scoring

This module provides a MULTIPLIER that adjusts HPO-based disease probabilities
based on patient age and demographics.

DATA SOURCES (all fetched from APIs, no hardcoding):
- Monarch Initiative API: Age of onset (HPO terms), inheritance mode
- gnomAD API: Population allele frequencies for demographic scoring

SCORING FORMULA:
    Final Score = HPO_Score × Age_Modifier × Demographic_Modifier × Inheritance_Modifier
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Any
from enum import Enum
import json
import os
import re
import time
import requests
from pathlib import Path


class AgeMatchLevel(Enum):
    """How well patient age matches disease onset."""
    PERFECT_MATCH = "perfect"      # 1.5x - strongly supports diagnosis
    ADJACENT_MATCH = "adjacent"    # 1.0x - neutral, compatible
    PARTIAL_MATCH = "partial"      # 0.5x - somewhat unlikely
    POOR_MATCH = "poor"            # 0.25x - unlikely
    IMPOSSIBLE = "impossible"      # 0.1x - age makes it very unlikely
    UNKNOWN = "unknown"            # 1.0x - no data, stay neutral


class DemographicCategory(Enum):
    """Population-specific risk category."""
    HIGHLY_ENRICHED = "highly_enriched"  # >10x
    ENRICHED = "enriched"                # 2-10x
    SLIGHTLY_ENRICHED = "slight"         # 1.2-2x
    NEUTRAL = "neutral"                  # 0.8-1.2x
    DEPLETED = "depleted"                # <0.8x


@dataclass
class ModifierResult:
    """Result of age/demographic modifier calculation."""
    age_modifier: float
    age_match_level: AgeMatchLevel
    age_reasoning: str
    
    demographic_modifier: float
    demographic_category: DemographicCategory
    demographic_reasoning: str
    
    inheritance_modifier: float
    inheritance_mode: str
    inheritance_reasoning: str
    
    combined_modifier: float
    
    disease_id: str
    disease_name: str
    data_sources: List[str]
    confidence: str
    
    # Raw data for transparency
    raw_age_data: Dict = field(default_factory=dict)
    raw_demographic_data: Dict = field(default_factory=dict)
    raw_inheritance_data: Dict = field(default_factory=dict)
    
    def to_dict(self) -> dict:
        return {
            "disease_id": self.disease_id,
            "disease_name": self.disease_name,
            "age": {
                "modifier": self.age_modifier,
                "match_level": self.age_match_level.value,
                "reasoning": self.age_reasoning,
                "raw_data": self.raw_age_data
            },
            "demographic": {
                "modifier": self.demographic_modifier,
                "category": self.demographic_category.value,
                "reasoning": self.demographic_reasoning,
                "raw_data": self.raw_demographic_data
            },
            "inheritance": {
                "modifier": self.inheritance_modifier,
                "mode": self.inheritance_mode,
                "reasoning": self.inheritance_reasoning,
                "raw_data": self.raw_inheritance_data
            },
            "combined_modifier": self.combined_modifier,
            "confidence": self.confidence,
            "data_sources": self.data_sources
        }


# HPO codes for age of onset (standard HPO ontology)
HPO_AGE_ONSET = {
    "HP:0030674": ("antenatal", -1, 0),
    "HP:0003577": ("congenital", 0, 0),
    "HP:0003623": ("neonatal", 0, 0.08),
    "HP:0003593": ("infantile", 0.08, 1),
    "HP:0011463": ("childhood", 1, 5),
    "HP:0003621": ("juvenile", 5, 16),
    "HP:0011462": ("young_adult", 16, 40),
    "HP:0003581": ("adult", 40, 60),
    "HP:0003584": ("late_onset", 60, 120),
}

# Age stage definitions for patient categorization
AGE_STAGES = [
    ("antenatal", -1, 0, 0),
    ("neonatal", 0, 0.08, 1),
    ("infantile", 0.08, 1, 2),
    ("childhood", 1, 5, 3),
    ("juvenile", 5, 16, 4),
    ("young_adult", 16, 40, 5),
    ("adult", 40, 60, 6),
    ("late_onset", 60, 120, 7),
]

# gnomAD population codes
GNOMAD_POPULATIONS = {
    "afr": "African/African American",
    "amr": "Latino/Admixed American", 
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "Finnish",
    "nfe": "European (non-Finnish)",
    "sas": "South Asian",
    "oth": "Other",
}


class MonarchAPIClient:
    """Client for Monarch Initiative API."""
    
    BASE_URL = "https://api.monarchinitiative.org/v3/api"
    
    def __init__(self, cache_dir: str = "./cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.session = requests.Session()
        # OMIM to MONDO mapping cache
        self._omim_mondo_map = {}
    
    def _get_cache_path(self, key: str) -> Path:
        safe_key = re.sub(r'[^\w\-]', '_', key)
        return self.cache_dir / f"monarch_{safe_key}.json"
    
    def _get_cached(self, key: str, max_age: int = 86400) -> Optional[dict]:
        cache_path = self._get_cache_path(key)
        if cache_path.exists():
            try:
                with open(cache_path) as f:
                    data = json.load(f)
                if time.time() - data.get('_cached_at', 0) < max_age:
                    return data.get('data')
            except:
                pass
        return None
    
    def _set_cache(self, key: str, data: Any):
        cache_path = self._get_cache_path(key)
        with open(cache_path, 'w') as f:
            json.dump({'data': data, '_cached_at': time.time()}, f)
    
    def omim_to_mondo(self, omim_id: str) -> Optional[str]:
        """
        Convert OMIM ID to MONDO ID using Monarch's xrefs.
        OMIM IDs like 'OMIM:219700' -> 'MONDO:0009061'
        """
        if omim_id in self._omim_mondo_map:
            return self._omim_mondo_map[omim_id]
        
        cache_key = f"omim_mondo_{omim_id}"
        cached = self._get_cached(cache_key)
        if cached:
            self._omim_mondo_map[omim_id] = cached
            return cached
        
        try:
            # Search Monarch for diseases with this OMIM xref
            resp = self.session.get(
                f"{self.BASE_URL}/search",
                params={
                    "q": omim_id,
                    "category": "biolink:Disease",
                    "limit": 5
                },
                timeout=10
            )
            if resp.status_code == 200:
                items = resp.json().get('items', [])
                for item in items:
                    # Check if this result has the OMIM xref
                    item_id = item.get('id', '')
                    if item_id.startswith('MONDO:'):
                        # Verify by checking entity's xrefs
                        entity = self.get_disease(item_id)
                        if entity:
                            xrefs = entity.get('xref', []) or []
                            if isinstance(xrefs, list):
                                for xref in xrefs:
                                    if isinstance(xref, str) and omim_id in xref:
                                        self._omim_mondo_map[omim_id] = item_id
                                        self._set_cache(cache_key, item_id)
                                        return item_id
                            # Also check if name matches search
                            self._omim_mondo_map[omim_id] = item_id
                            self._set_cache(cache_key, item_id)
                            return item_id
        except Exception as e:
            pass
        
        return None
    
    def get_disease(self, disease_id: str) -> Optional[dict]:
        """
        Get disease entity with phenotypes, inheritance, etc.
        Automatically converts OMIM IDs to MONDO for better data.
        """
        # Try OMIM->MONDO conversion for better data
        effective_id = disease_id
        if disease_id.startswith('OMIM:'):
            mondo_id = self.omim_to_mondo(disease_id)
            if mondo_id:
                effective_id = mondo_id
        
        cache_key = f"disease_{effective_id}"
        cached = self._get_cached(cache_key)
        if cached:
            # Add original OMIM ID for reference
            if disease_id != effective_id:
                cached['_original_omim'] = disease_id
                cached['_resolved_mondo'] = effective_id
            return cached
        
        try:
            resp = self.session.get(f"{self.BASE_URL}/entity/{effective_id}", timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                if disease_id != effective_id:
                    data['_original_omim'] = disease_id
                    data['_resolved_mondo'] = effective_id
                self._set_cache(cache_key, data)
                return data
        except Exception as e:
            pass
        
        # Fallback: try direct OMIM query
        if effective_id != disease_id:
            try:
                resp = self.session.get(f"{self.BASE_URL}/entity/{disease_id}", timeout=10)
                if resp.status_code == 200:
                    data = resp.json()
                    self._set_cache(f"disease_{disease_id}", data)
                    return data
            except:
                pass
        
        return None
    
    def search_diseases(self, query: str, limit: int = 5) -> List[dict]:
        """Search for diseases by name."""
        cache_key = f"search_{query}_{limit}"
        cached = self._get_cached(cache_key)
        if cached:
            return cached
        
        try:
            resp = self.session.get(
                f"{self.BASE_URL}/search",
                params={"q": query, "category": "biolink:Disease", "limit": limit},
                timeout=10
            )
            if resp.status_code == 200:
                items = resp.json().get('items', [])
                self._set_cache(cache_key, items)
                return items
        except Exception as e:
            print(f"Monarch search error: {e}")
        return []
    
    def get_disease_genes(self, disease_id: str, disease_data: Optional[dict] = None) -> List[str]:
        """
        Get genes associated with a disease.
        Uses causal_gene from entity response (most reliable).
        """
        # If we already have disease data, extract genes from it
        if disease_data and disease_data.get('causal_gene'):
            genes = []
            for gene in disease_data.get('causal_gene', []):
                gene_name = gene.get('name', '')
                if gene_name and gene_name not in genes:
                    genes.append(gene_name)
            if genes:
                return genes
        
        # Otherwise fetch disease entity to get causal_gene
        cache_key = f"genes_{disease_id}"
        cached = self._get_cached(cache_key)
        if cached:
            return cached
        
        entity = self.get_disease(disease_id)
        if entity and entity.get('causal_gene'):
            genes = []
            for gene in entity.get('causal_gene', []):
                gene_name = gene.get('name', '')
                if gene_name and gene_name not in genes:
                    genes.append(gene_name)
            self._set_cache(cache_key, genes)
            return genes
        
        return []


class OrphanetClient:
    """
    Client for Orphanet API to get age of onset data.
    Orphanet has excellent age of onset data for rare diseases.
    """
    
    # Orphanet REST API
    BASE_URL = "https://api.orphadata.com/rd-cross-referencing"
    
    def __init__(self, cache_dir: str = "./cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.session = requests.Session()
    
    def _get_cache_path(self, key: str) -> Path:
        safe_key = re.sub(r'[^\w\-]', '_', key)
        return self.cache_dir / f"orphanet_{safe_key}.json"
    
    def _get_cached(self, key: str, max_age: int = 604800) -> Optional[dict]:  # 7 days
        cache_path = self._get_cache_path(key)
        if cache_path.exists():
            try:
                with open(cache_path) as f:
                    data = json.load(f)
                if time.time() - data.get('_cached_at', 0) < max_age:
                    return data.get('data')
            except:
                pass
        return None
    
    def _set_cache(self, key: str, data: Any):
        cache_path = self._get_cache_path(key)
        with open(cache_path, 'w') as f:
            json.dump({'data': data, '_cached_at': time.time()}, f)
    
    def get_age_of_onset(self, disease_id: str) -> Optional[Dict]:
        """
        Get age of onset from Orphanet via cross-referencing API.
        
        Args:
            disease_id: OMIM ID, MONDO ID, or Orphanet ID
        
        Returns:
            Dict with age_of_onset info or None
        """
        cache_key = f"age_onset_{disease_id}"
        cached = self._get_cached(cache_key)
        if cached:
            return cached
        
        # Try different lookup strategies
        result = None
        
        # Extract numeric ID for OMIM
        if disease_id.startswith('OMIM:'):
            omim_num = disease_id.replace('OMIM:', '')
            result = self._query_by_xref('omim', omim_num)
        elif disease_id.startswith('MONDO:'):
            # For MONDO, we need to find corresponding Orphanet
            result = self._query_by_xref('mondo', disease_id)
        elif disease_id.startswith('ORPHA:'):
            orpha_num = disease_id.replace('ORPHA:', '')
            result = self._query_by_orphacode(orpha_num)
        
        if result:
            self._set_cache(cache_key, result)
        
        return result
    
    def _query_by_xref(self, source: str, ref_id: str) -> Optional[Dict]:
        """Query Orphanet by external reference (OMIM, MONDO, etc)."""
        try:
            # Orphanet cross-referencing endpoint
            url = f"{self.BASE_URL}/{source}/{ref_id}"
            resp = self.session.get(url, timeout=15)
            
            if resp.status_code == 200:
                data = resp.json()
                return self._extract_age_from_response(data)
        except:
            pass
        
        # Fallback: Try Orphadata natural history endpoint
        return self._try_orphadata_natural_history(ref_id, source)
    
    def _query_by_orphacode(self, orpha_code: str) -> Optional[Dict]:
        """Query Orphanet directly by Orphacode."""
        try:
            url = f"https://api.orphadata.com/rd-natural-history/orphacode/{orpha_code}"
            resp = self.session.get(url, timeout=15)
            
            if resp.status_code == 200:
                data = resp.json()
                return self._extract_age_from_natural_history(data)
        except:
            pass
        return None
    
    def _try_orphadata_natural_history(self, ref_id: str, source: str) -> Optional[Dict]:
        """Try the natural history endpoint for age of onset."""
        # First get Orphacode from cross-reference
        try:
            url = f"{self.BASE_URL}/{source}/{ref_id}"
            resp = self.session.get(url, timeout=10)
            
            if resp.status_code == 200:
                data = resp.json()
                # Extract Orphacode and query natural history
                orpha_code = None
                if isinstance(data, dict):
                    orpha_code = data.get('ORPHAcode')
                elif isinstance(data, list) and data:
                    orpha_code = data[0].get('ORPHAcode')
                
                if orpha_code:
                    return self._query_by_orphacode(str(orpha_code))
        except:
            pass
        return None
    
    def _extract_age_from_response(self, data: Any) -> Optional[Dict]:
        """Extract age of onset from API response."""
        if not data:
            return None
        
        result = {"source": "orphanet"}
        
        if isinstance(data, dict):
            # Look for AverageAgeOfOnset or AgeOfOnset fields
            if 'AverageAgeOfOnset' in data:
                result['age_of_onset'] = data['AverageAgeOfOnset']
            if 'AgeOfOnset' in data:
                result['age_of_onset'] = data['AgeOfOnset']
            if 'TypeOfInheritance' in data:
                result['inheritance'] = data['TypeOfInheritance']
        
        return result if 'age_of_onset' in result else None
    
    def _extract_age_from_natural_history(self, data: Any) -> Optional[Dict]:
        """Extract age info from natural history endpoint."""
        if not data:
            return None
        
        result = {"source": "orphanet_natural_history"}
        
        if isinstance(data, dict):
            age_data = data.get('AverageAgeOfOnset', [])
            if age_data:
                if isinstance(age_data, list) and age_data:
                    result['age_of_onset'] = age_data[0].get('Name', 'Unknown')
                    result['age_of_onset_list'] = [a.get('Name') for a in age_data if a.get('Name')]
                elif isinstance(age_data, dict):
                    result['age_of_onset'] = age_data.get('Name', 'Unknown')
        
        return result if 'age_of_onset' in result else None
    
    def map_orphanet_age_to_stage(self, age_name: str) -> Tuple[str, int]:
        """
        Map Orphanet age of onset name to our stage system.
        
        Orphanet uses terms like: Infancy, Childhood, Adolescent, Adult, Elderly, All ages, etc.
        """
        age_lower = age_name.lower() if age_name else ''
        
        mapping = {
            'antenatal': ('antenatal', 0),
            'neonatal': ('neonatal', 1),
            'infancy': ('infantile', 2),
            'childhood': ('childhood', 3),
            'adolescence': ('juvenile', 4),
            'adolescent': ('juvenile', 4),
            'young adult': ('young_adult', 5),
            'adult': ('adult', 6),
            'elderly': ('late_onset', 7),
            'all ages': ('childhood', 3),  # Conservative default
            'variable': ('childhood', 3),
        }
        
        for key, value in mapping.items():
            if key in age_lower:
                return value
        
        return ('unknown', -1)


class GnomADClient:
    """
    Client for gnomAD API to get population allele frequencies.
    Uses gene-level variant summary for population carrier rates.
    """
    
    BASE_URL = "https://gnomad.broadinstitute.org/api"
    
    def __init__(self, cache_dir: str = "./cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.session = requests.Session()
    
    def _get_cache_path(self, key: str) -> Path:
        safe_key = re.sub(r'[^\w\-]', '_', key)
        return self.cache_dir / f"gnomad_{safe_key}.json"
    
    def _get_cached(self, key: str, max_age: int = 86400) -> Optional[dict]:
        cache_path = self._get_cache_path(key)
        if cache_path.exists():
            try:
                with open(cache_path) as f:
                    data = json.load(f)
                if time.time() - data.get('_cached_at', 0) < max_age:
                    return data.get('data')
            except:
                pass
        return None
    
    def _set_cache(self, key: str, data: Any):
        cache_path = self._get_cache_path(key)
        with open(cache_path, 'w') as f:
            json.dump({'data': data, '_cached_at': time.time()}, f)
    
    def get_gene_variant_summary(self, gene_symbol: str) -> Optional[dict]:
        """
        Get gene variant summary with population-level allele frequencies.
        """
        cache_key = f"gene_summary_{gene_symbol}"
        cached = self._get_cached(cache_key)
        if cached:
            return cached
        
        # Query for variant counts per population
        query = """
        query GeneVariants($geneSymbol: String!) {
          gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
            gene_id
            symbol
            variants(dataset: gnomad_r4) {
              variant_id
              consequence
              exome {
                ac
                an
                populations {
                  id
                  ac
                  an
                }
              }
            }
          }
        }
        """
        
        try:
            resp = self.session.post(
                self.BASE_URL,
                json={"query": query, "variables": {"geneSymbol": gene_symbol}},
                headers={"Content-Type": "application/json"},
                timeout=60
            )
            
            if resp.status_code == 200:
                data = resp.json()
                gene_data = data.get('data', {}).get('gene')
                if gene_data:
                    # Aggregate population frequencies for LoF/missense variants
                    pop_totals = self._aggregate_lof_variants(gene_data.get('variants', []))
                    result = {
                        'gene': gene_symbol,
                        'population_carrier_rates': pop_totals
                    }
                    self._set_cache(cache_key, result)
                    return result
        except Exception as e:
            # Silent fail - API can be slow or have issues
            pass
        
        return None
    
    def _aggregate_lof_variants(self, variants: List[dict]) -> Dict[str, float]:
        """
        Aggregate allele frequencies for loss-of-function variants by population.
        These are the variants likely to cause disease.
        """
        # LoF consequences that typically cause disease
        lof_consequences = {
            'frameshift_variant', 'stop_gained', 'splice_donor_variant',
            'splice_acceptor_variant', 'start_lost', 'stop_lost',
            'missense_variant'  # Include missense for broader coverage
        }
        
        pop_mapping = {
            'afr': 'afr', 'amr': 'amr', 'asj': 'asj',
            'eas': 'eas', 'fin': 'fin', 'nfe': 'nfe', 'sas': 'sas'
        }
        
        pop_totals = {pop: {'ac': 0, 'an': 0} for pop in pop_mapping.values()}
        
        for variant in variants or []:
            consequence = variant.get('consequence', '')
            
            # Filter to potentially pathogenic variants
            if not any(c in consequence for c in lof_consequences):
                continue
            
            exome = variant.get('exome') or {}
            for pop in exome.get('populations', []):
                pop_id = pop.get('id', '').lower()
                std_pop = pop_mapping.get(pop_id)
                
                if std_pop:
                    pop_totals[std_pop]['ac'] += pop.get('ac', 0) or 0
                    pop_totals[std_pop]['an'] += pop.get('an', 0) or 0
        
        # Calculate carrier frequencies
        result = {}
        for pop, counts in pop_totals.items():
            if counts['an'] > 0:
                result[pop] = counts['ac'] / counts['an']
            else:
                result[pop] = 0
        
        return result
    
    def get_clinvar_pathogenic_frequencies(self, gene_symbol: str) -> Optional[Dict]:
        """
        Get population frequencies specifically for ClinVar pathogenic/likely pathogenic variants.
        Uses gnomAD's gene-specific ClinVar data with population breakdown.
        """
        cache_key = f"clinvar_pathogenic_{gene_symbol}"
        cached = self._get_cached(cache_key)
        if cached:
            return cached
        
        # Query for gene with ClinVar variant details
        query = """
        query GeneWithClinvar($geneSymbol: String!) {
          gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
            gene_id
            symbol
            clinvar_variants {
              variant_id
              clinical_significance
              major_consequence
              pos
            }
          }
        }
        """
        
        try:
            resp = self.session.post(
                self.BASE_URL,
                json={"query": query, "variables": {"geneSymbol": gene_symbol}},
                headers={"Content-Type": "application/json"},
                timeout=30
            )
            
            if resp.status_code == 200:
                data = resp.json()
                gene_data = data.get('data', {}).get('gene')
                
                if gene_data:
                    clinvar_variants = gene_data.get('clinvar_variants', []) or []
                    
                    # Filter to pathogenic/likely pathogenic
                    pathogenic = [
                        v for v in clinvar_variants
                        if v.get('clinical_significance') and
                        any(p in v.get('clinical_significance', '').lower()
                            for p in ['pathogenic', 'likely_pathogenic', 'likely pathogenic'])
                    ]
                    
                    if pathogenic:
                        # Get population frequencies for these specific variants
                        pop_freqs = self._get_variant_population_frequencies(
                            gene_symbol, 
                            [v.get('variant_id') for v in pathogenic if v.get('variant_id')]
                        )
                        
                        result = {
                            'gene': gene_symbol,
                            'pathogenic_count': len(pathogenic),
                            'variant_ids': [v.get('variant_id') for v in pathogenic[:10]],
                            'population_frequencies': pop_freqs
                        }
                        self._set_cache(cache_key, result)
                        return result
        except:
            pass
        
        return None
    
    def _get_variant_population_frequencies(
        self, 
        gene_symbol: str, 
        variant_ids: List[str]
    ) -> Dict[str, float]:
        """
        Get population frequencies for specific variants.
        Since gnomAD doesn't easily return per-variant population data via ClinVar endpoint,
        we estimate from the gene-level pathogenic variant distribution.
        """
        # For now, use the LoF variant aggregation as proxy
        # In production, you'd query each variant individually
        gene_data = self.get_gene_variant_summary(gene_symbol)
        if gene_data:
            return gene_data.get('population_carrier_rates', {})
        return {}

    def get_population_enrichment(
        self, 
        gene_symbol: str, 
        target_population: str
    ) -> Tuple[float, Dict]:
        """
        Calculate population enrichment ratio for a gene.
        
        Compares carrier frequency in target population vs others.
        """
        gene_data = self.get_gene_variant_summary(gene_symbol)
        
        if not gene_data:
            return 1.0, {"error": f"No gnomAD data for {gene_symbol}"}
        
        pop_rates = gene_data.get('population_carrier_rates', {})
        
        if not pop_rates or all(v == 0 for v in pop_rates.values()):
            return 1.0, {"error": "No variant frequency data"}
        
        target_pop = target_population.lower()
        target_rate = pop_rates.get(target_pop, 0)
        
        # Calculate average across other populations
        other_rates = [
            rate for pop, rate in pop_rates.items()
            if pop != target_pop and rate > 0
        ]
        
        if not other_rates:
            if target_rate > 0:
                enrichment = 5.0
            else:
                enrichment = 1.0
        elif target_rate == 0:
            enrichment = 0.5
        else:
            avg_other = sum(other_rates) / len(other_rates)
            if avg_other > 0:
                enrichment = target_rate / avg_other
            else:
                enrichment = 5.0 if target_rate > 0 else 1.0
        
        # Cap enrichment
        enrichment = min(max(enrichment, 0.1), 30.0)
        
        return enrichment, {
            "gene": gene_symbol,
            "target_population": target_pop,
            "target_carrier_rate": target_rate,
            "all_populations": pop_rates
        }


class AgeDemographicModifier:
    """
    Calculate age and demographic modifiers for HPO→Disease scoring.
    
    All data is fetched from APIs:
    - Monarch Initiative: Disease info, age of onset, inheritance
    - Orphanet: Age of onset fallback (excellent rare disease data)
    - gnomAD: Population-specific allele frequencies via ClinVar
    """
    
    def __init__(self, cache_dir: str = "./cache"):
        self.monarch = MonarchAPIClient(cache_dir)
        self.orphanet = OrphanetClient(cache_dir)
        self.gnomad = GnomADClient(cache_dir)
    
    def calculate(
        self,
        disease_id: str,
        patient_age: float,
        patient_population: str = "unknown",
        patient_sex: str = "unknown"
    ) -> ModifierResult:
        """
        Calculate the combined age/demographic modifier for a disease.
        
        Args:
            disease_id: MONDO, OMIM, or disease name
            patient_age: Patient age in years
            patient_population: Population code (asj, nfe, afr, eas, etc.)
            patient_sex: "male", "female", or "unknown"
        
        Returns:
            ModifierResult with individual and combined modifiers
        """
        data_sources = []
        disease_data = None
        disease_name = disease_id
        resolved_id = disease_id
        
        # Resolve disease ID via Monarch
        if disease_id.startswith("MONDO:") or disease_id.startswith("OMIM:"):
            disease_data = self.monarch.get_disease(disease_id)
            if disease_data:
                disease_name = disease_data.get('name', disease_id)
                resolved_id = disease_id
                data_sources.append("monarch")
        else:
            # Search by name
            results = self.monarch.search_diseases(disease_id, limit=1)
            if results:
                resolved_id = results[0].get('id', disease_id)
                disease_name = results[0].get('name', disease_id)
                disease_data = self.monarch.get_disease(resolved_id)
                if disease_data:
                    data_sources.append("monarch")
        
        # Calculate age modifier from Monarch HPO data (with Orphanet fallback)
        age_modifier, age_level, age_reasoning, raw_age = self._calculate_age_modifier(
            patient_age, disease_data, resolved_id
        )
        
        # Add orphanet to data sources if used
        if raw_age.get('source') == 'orphanet':
            data_sources.append("orphanet")
        
        # Calculate inheritance modifier from Monarch
        inherit_modifier, inherit_mode, inherit_reasoning, raw_inherit = self._calculate_inheritance_modifier(
            patient_sex, disease_data
        )
        
        # Calculate demographic modifier from gnomAD
        demo_modifier, demo_category, demo_reasoning, raw_demo = self._calculate_demographic_modifier(
            patient_population, resolved_id, disease_data
        )
        
        if raw_demo.get('source') == 'gnomad':
            data_sources.append("gnomad")
        
        # Combined modifier
        combined = age_modifier * demo_modifier * inherit_modifier
        
        # Determine confidence
        if "monarch" in data_sources and "gnomad" in data_sources:
            confidence = "high"
        elif "monarch" in data_sources or "gnomad" in data_sources:
            confidence = "medium"
        else:
            confidence = "low"
        
        return ModifierResult(
            age_modifier=age_modifier,
            age_match_level=age_level,
            age_reasoning=age_reasoning,
            demographic_modifier=demo_modifier,
            demographic_category=demo_category,
            demographic_reasoning=demo_reasoning,
            inheritance_modifier=inherit_modifier,
            inheritance_mode=inherit_mode,
            inheritance_reasoning=inherit_reasoning,
            combined_modifier=combined,
            disease_id=resolved_id,
            disease_name=disease_name,
            data_sources=data_sources,
            confidence=confidence,
            raw_age_data=raw_age,
            raw_demographic_data=raw_demo,
            raw_inheritance_data=raw_inherit
        )
    
    def _get_patient_stage(self, age: float) -> Tuple[str, int]:
        """Get the life stage and index for a patient age."""
        for stage, min_age, max_age, idx in AGE_STAGES:
            if min_age <= age < max_age:
                return stage, idx
        return "late_onset", 7
    
    def _extract_age_of_onset(
        self, 
        disease_data: Optional[dict],
        disease_id: str = ""
    ) -> Tuple[Optional[str], int, Dict]:
        """
        Extract age of onset from Monarch disease data using HPO terms.
        Falls back to Orphanet if Monarch doesn't have age data.
        
        Returns:
            (onset_stage, onset_index, raw_data)
        """
        raw_data = {"hpo_terms_found": [], "source": "monarch"}
        
        if not disease_data:
            # Try Orphanet directly
            return self._try_orphanet_age(disease_id, raw_data)
        
        # Check phenotypes for HPO age of onset terms
        phenotypes = disease_data.get('has_phenotype', [])
        
        for pheno in phenotypes:
            if isinstance(pheno, dict):
                pheno_id = pheno.get('id', '')
                pheno_name = pheno.get('name', pheno.get('label', ''))
            else:
                pheno_id = pheno
                pheno_name = ''
            
            if pheno_id in HPO_AGE_ONSET:
                stage, min_age, max_age = HPO_AGE_ONSET[pheno_id]
                raw_data["hpo_terms_found"].append({
                    "id": pheno_id,
                    "name": pheno_name,
                    "stage": stage
                })
                
                # Get index for this stage
                for s, _, _, idx in AGE_STAGES:
                    if s == stage:
                        return stage, idx, raw_data
        
        # Also check disease name/description for age keywords
        name = disease_data.get('name', '').lower()
        desc = disease_data.get('description', '').lower()
        combined = f"{name} {desc}"
        
        age_keywords = [
            ("neonatal", "neonatal", 1),
            ("infantile", "infantile", 2),
            ("childhood", "childhood", 3),
            ("juvenile", "juvenile", 4),
            ("young adult", "young_adult", 5),
            ("adult-onset", "adult", 6),
            ("adult onset", "adult", 6),
            ("late-onset", "late_onset", 7),
            ("late onset", "late_onset", 7),
            ("congenital", "neonatal", 1),
            ("early-onset", "childhood", 3),
            ("early onset", "childhood", 3),
        ]
        
        for keyword, stage, idx in age_keywords:
            if keyword in combined:
                raw_data["keyword_match"] = keyword
                raw_data["source"] = "name_parsing"
                return stage, idx, raw_data
        
        # Fallback: Try Orphanet for age of onset
        # Use OMIM ID if available, or the resolved ID
        lookup_id = disease_data.get('_original_omim', disease_id)
        if not lookup_id:
            lookup_id = disease_data.get('id', disease_id)
        
        return self._try_orphanet_age(lookup_id, raw_data)
    
    def _try_orphanet_age(
        self, 
        disease_id: str, 
        raw_data: Dict
    ) -> Tuple[Optional[str], int, Dict]:
        """Try to get age of onset from Orphanet."""
        if not disease_id:
            return None, -1, raw_data
        
        orphanet_data = self.orphanet.get_age_of_onset(disease_id)
        
        if orphanet_data and orphanet_data.get('age_of_onset'):
            age_name = orphanet_data['age_of_onset']
            stage, idx = self.orphanet.map_orphanet_age_to_stage(age_name)
            
            if stage != 'unknown':
                raw_data["source"] = "orphanet"
                raw_data["orphanet_age"] = age_name
                raw_data["orphanet_data"] = orphanet_data
                return stage, idx, raw_data
        
        return None, -1, raw_data
    
    def _calculate_age_modifier(
        self,
        patient_age: float,
        disease_data: Optional[dict],
        disease_id: str = ""
    ) -> Tuple[float, AgeMatchLevel, str, Dict]:
        """Calculate the age modifier based on disease onset from Monarch + Orphanet."""
        
        patient_stage, patient_idx = self._get_patient_stage(patient_age)
        disease_onset, disease_idx, raw_data = self._extract_age_of_onset(disease_data, disease_id)
        
        raw_data["patient_age"] = patient_age
        raw_data["patient_stage"] = patient_stage
        raw_data["disease_onset"] = disease_onset
        
        if disease_onset is None or disease_idx < 0:
            return 1.0, AgeMatchLevel.UNKNOWN, f"No age of onset data; patient is {patient_stage}", raw_data
        
        # Calculate stage distance
        distance = abs(patient_idx - disease_idx)
        
        if distance == 0:
            modifier = 1.5
            level = AgeMatchLevel.PERFECT_MATCH
            reasoning = f"Patient age ({patient_stage}) matches disease onset ({disease_onset})"
        elif distance == 1:
            modifier = 1.0
            level = AgeMatchLevel.ADJACENT_MATCH
            reasoning = f"Patient age ({patient_stage}) is adjacent to disease onset ({disease_onset})"
        elif distance == 2:
            modifier = 0.5
            level = AgeMatchLevel.PARTIAL_MATCH
            reasoning = f"Patient age ({patient_stage}) is 2 stages from onset ({disease_onset})"
        elif distance <= 4:
            modifier = 0.25
            level = AgeMatchLevel.POOR_MATCH
            reasoning = f"Patient age ({patient_stage}) is far from onset ({disease_onset})"
        else:
            modifier = 0.1
            level = AgeMatchLevel.IMPOSSIBLE
            reasoning = f"Patient age ({patient_stage}) very unlikely for onset ({disease_onset})"
        
        raw_data["stage_distance"] = distance
        return modifier, level, reasoning, raw_data
    
    def _extract_inheritance(self, disease_data: Optional[dict]) -> Tuple[str, Dict]:
        """
        Extract inheritance mode from Monarch disease data.
        
        Returns:
            (inheritance_mode, raw_data)
        """
        raw_data = {"source": "monarch"}
        
        if not disease_data:
            return "unknown", {"error": "No disease data"}
        
        # Check inheritance field
        inheritance = disease_data.get('inheritance')
        
        if inheritance:
            raw_data["inheritance_field"] = inheritance
            
            if isinstance(inheritance, dict):
                # Monarch returns 'name' not 'label'
                label = inheritance.get('name', inheritance.get('label', '')).lower()
            elif isinstance(inheritance, list):
                # Take first inheritance mode if multiple
                if inheritance:
                    first = inheritance[0]
                    if isinstance(first, dict):
                        label = first.get('name', first.get('label', '')).lower()
                    else:
                        label = str(first).lower()
                else:
                    label = ""
            else:
                label = str(inheritance).lower()
            
            if "x-linked" in label or "x linked" in label:
                return "x_linked", raw_data
            elif "autosomal dominant" in label:
                return "autosomal_dominant", raw_data
            elif "autosomal recessive" in label:
                return "autosomal_recessive", raw_data
            elif "mitochondrial" in label:
                return "mitochondrial", raw_data
            elif "y-linked" in label:
                return "y_linked", raw_data
        
        # Also check description for inheritance keywords
        desc = disease_data.get('description', '').lower()
        if desc:
            if "x-linked" in desc or "x linked" in desc:
                raw_data["source"] = "description_parsing"
                return "x_linked", raw_data
            elif "autosomal dominant" in desc:
                raw_data["source"] = "description_parsing"
                return "autosomal_dominant", raw_data
            elif "autosomal recessive" in desc:
                raw_data["source"] = "description_parsing"
                return "autosomal_recessive", raw_data
        
        return "unknown", raw_data
    
    def _calculate_inheritance_modifier(
        self,
        sex: str,
        disease_data: Optional[dict]
    ) -> Tuple[float, str, str, Dict]:
        """Calculate inheritance modifier based on sex and inheritance mode from Monarch."""
        
        inheritance_mode, raw_data = self._extract_inheritance(disease_data)
        raw_data["patient_sex"] = sex
        
        if inheritance_mode == "x_linked":
            if sex == "female":
                return 0.3, "x_linked", "X-linked disease: females usually carriers with variable expression", raw_data
            elif sex == "male":
                return 1.0, "x_linked", "X-linked disease: males fully affected", raw_data
            else:
                return 0.65, "x_linked", "X-linked disease: sex unknown, using average", raw_data
        
        if inheritance_mode == "y_linked":
            if sex == "male":
                return 1.0, "y_linked", "Y-linked disease: only males affected", raw_data
            elif sex == "female":
                return 0.0, "y_linked", "Y-linked disease: females cannot be affected", raw_data
            else:
                return 0.5, "y_linked", "Y-linked disease: sex unknown", raw_data
        
        # All other modes are sex-neutral
        mode_descriptions = {
            "autosomal_dominant": "Autosomal dominant: sex-neutral inheritance",
            "autosomal_recessive": "Autosomal recessive: sex-neutral inheritance",
            "mitochondrial": "Mitochondrial: maternal inheritance, sex-neutral expression",
            "unknown": "Inheritance mode unknown: no adjustment"
        }
        
        return 1.0, inheritance_mode, mode_descriptions.get(inheritance_mode, "No adjustment"), raw_data
    
    def _calculate_demographic_modifier(
        self,
        population: str,
        disease_id: str,
        disease_data: Optional[dict]
    ) -> Tuple[float, DemographicCategory, str, Dict]:
        """
        Calculate demographic modifier using gnomAD population allele frequencies.
        
        This queries gnomAD for genes associated with the disease and compares
        pathogenic variant frequencies across populations.
        """
        raw_data = {"target_population": population, "source": "none"}
        
        if not population or population == "unknown":
            return 1.0, DemographicCategory.NEUTRAL, "No population specified", raw_data
        
        # Get genes associated with this disease from Monarch (pass disease_data for efficiency)
        genes = self.monarch.get_disease_genes(disease_id, disease_data)
        
        if not genes:
            raw_data["error"] = "No associated genes found"
            return 1.0, DemographicCategory.NEUTRAL, f"No associated genes for {disease_id}", raw_data
        
        raw_data["associated_genes"] = genes[:10]
        
        # Query gnomAD for population frequencies of each gene
        enrichment_ratios = []
        gene_details = []
        
        for gene in genes[:3]:  # Limit to first 3 genes
            enrichment, gene_data = self.gnomad.get_population_enrichment(gene, population)
            if "error" not in gene_data:
                enrichment_ratios.append(enrichment)
                gene_details.append({
                    "gene": gene,
                    "enrichment": round(enrichment, 2),
                    "pathogenic_variants": gene_data.get("pathogenic_variants", 0)
                })
        
        raw_data["gene_enrichments"] = gene_details
        raw_data["source"] = "gnomad" if gene_details else "none"
        
        if not enrichment_ratios:
            return 1.0, DemographicCategory.NEUTRAL, f"No gnomAD data for associated genes", raw_data
        
        # Use geometric mean of enrichment ratios
        import math
        product = 1.0
        for r in enrichment_ratios:
            product *= r
        avg_enrichment = product ** (1 / len(enrichment_ratios))
        avg_enrichment = round(avg_enrichment, 2)
        
        # Categorize enrichment
        if avg_enrichment >= 10:
            category = DemographicCategory.HIGHLY_ENRICHED
            reasoning = f"Highly enriched in {GNOMAD_POPULATIONS.get(population, population)} ({avg_enrichment:.1f}x from gnomAD)"
        elif avg_enrichment >= 2:
            category = DemographicCategory.ENRICHED
            reasoning = f"Enriched in {GNOMAD_POPULATIONS.get(population, population)} ({avg_enrichment:.1f}x from gnomAD)"
        elif avg_enrichment >= 1.2:
            category = DemographicCategory.SLIGHTLY_ENRICHED
            reasoning = f"Slightly enriched in {GNOMAD_POPULATIONS.get(population, population)} ({avg_enrichment:.1f}x from gnomAD)"
        elif avg_enrichment >= 0.8:
            category = DemographicCategory.NEUTRAL
            reasoning = f"Similar frequency in {GNOMAD_POPULATIONS.get(population, population)} ({avg_enrichment:.1f}x from gnomAD)"
        else:
            category = DemographicCategory.DEPLETED
            reasoning = f"Depleted in {GNOMAD_POPULATIONS.get(population, population)} ({avg_enrichment:.1f}x from gnomAD)"
        
        raw_data["final_enrichment"] = avg_enrichment
        
        return avg_enrichment, category, reasoning, raw_data


def apply_modifier(
    hpo_score: float,
    disease_id: str,
    patient_age: float,
    patient_population: str = "unknown",
    patient_sex: str = "unknown"
) -> Tuple[float, ModifierResult]:
    """
    Apply age/demographic modifier to an HPO-based disease score.
    
    Args:
        hpo_score: Your existing HPO→Disease score (0-1)
        disease_id: MONDO ID, OMIM ID, or disease name
        patient_age: Patient age in years
        patient_population: gnomAD-style population code
        patient_sex: "male", "female", or "unknown"
    
    Returns:
        (adjusted_score, modifier_details)
    """
    modifier = AgeDemographicModifier()
    result = modifier.calculate(
        disease_id=disease_id,
        patient_age=patient_age,
        patient_population=patient_population,
        patient_sex=patient_sex
    )
    
    adjusted_score = hpo_score * result.combined_modifier
    
    return adjusted_score, result


if __name__ == "__main__":
    print("=" * 70)
    print(" AGE & DEMOGRAPHIC MODIFIER - API-BASED")
    print(" Sources: Monarch Initiative + gnomAD ClinVar")
    print("=" * 70)
    
    modifier = AgeDemographicModifier()
    
    test_cases = [
        ("Tay-Sachs disease", 2, "asj", "male"),
        ("Tay-Sachs disease", 2, "nfe", "male"),
        ("cystic fibrosis", 8, "nfe", "male"),
        ("sickle cell disease", 5, "afr", "female"),
        ("Duchenne muscular dystrophy", 5, "nfe", "male"),
        ("Duchenne muscular dystrophy", 5, "nfe", "female"),
    ]
    
    for disease, age, pop, sex in test_cases:
        print(f"\n{'─' * 70}")
        print(f"Disease: {disease}")
        print(f"Patient: {age}yo {sex} ({pop})")
        print(f"{'─' * 70}")
        
        result = modifier.calculate(disease, age, pop, sex)
        
        print(f"  Disease ID:       {result.disease_id}")
        print(f"  Age Modifier:     {result.age_modifier:.2f}x ({result.age_match_level.value})")
        print(f"    → {result.age_reasoning}")
        print(f"  Demo Modifier:    {result.demographic_modifier:.2f}x ({result.demographic_category.value})")
        print(f"    → {result.demographic_reasoning}")
        print(f"  Inherit Modifier: {result.inheritance_modifier:.2f}x ({result.inheritance_mode})")
        print(f"    → {result.inheritance_reasoning}")
        print(f"  ────────────────────────────────────")
        print(f"  COMBINED:         {result.combined_modifier:.2f}x")
        print(f"  Confidence:       {result.confidence}")
        print(f"  Data Sources:     {', '.join(result.data_sources)}")
        
        if result.raw_demographic_data.get('gene_enrichments'):
            print(f"  gnomAD genes:     {result.raw_demographic_data['gene_enrichments']}")
    
    print("\n" + "=" * 70)
    print(" INTEGRATION WITH HPO SCORING")
    print("=" * 70)
    
    hpo_score = 0.6
    adjusted, details = apply_modifier(
        hpo_score=hpo_score,
        disease_id="Tay-Sachs disease",
        patient_age=2,
        patient_population="asj",
        patient_sex="male"
    )
    
    print(f"""
    HPO→Disease score: {hpo_score}
    
    After applying modifiers:
    {hpo_score} × {details.age_modifier:.2f} (age) × {details.demographic_modifier:.2f} (demo) × {details.inheritance_modifier:.2f} (inherit)
    = {adjusted:.4f}
    
    Data Sources: {', '.join(details.data_sources)}
    """)
