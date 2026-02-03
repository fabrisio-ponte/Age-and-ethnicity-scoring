"""
Monarch API Client
==================

Fetches disease data from the Monarch Initiative API.
https://api.monarchinitiative.org

Monarch combines data from:
- OMIM (disease genetics)
- Orphanet (rare disease prevalence)
- HPO (phenotypes and age of onset)
- gnomAD (population genetics)

No downloads needed - all on-demand API calls with local caching.
"""

import json
import requests
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, List, Optional

# Cache directory
CACHE_DIR = Path(__file__).parent / "cache"
CACHE_DIR.mkdir(exist_ok=True)
CACHE_EXPIRY_HOURS = 24

# Monarch API endpoints
MONARCH_BASE = "https://api.monarchinitiative.org/v3/api"


class MonarchClient:
    """
    Client for the Monarch Initiative API.
    """
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            "Accept": "application/json"
        })
    
    def _get_cache_path(self, key: str) -> Path:
        """Get cache file path."""
        safe_key = key.replace(":", "_").replace("/", "_")
        return CACHE_DIR / f"{safe_key}.json"
    
    def _get_cached(self, key: str) -> Optional[dict]:
        """Get cached response if valid."""
        cache_path = self._get_cache_path(key)
        if not cache_path.exists():
            return None
        
        # Check expiry
        modified = datetime.fromtimestamp(cache_path.stat().st_mtime)
        if datetime.now() - modified > timedelta(hours=CACHE_EXPIRY_HOURS):
            return None
        
        with open(cache_path, 'r') as f:
            return json.load(f)
    
    def _set_cached(self, key: str, data: dict):
        """Cache response."""
        cache_path = self._get_cache_path(key)
        with open(cache_path, 'w') as f:
            json.dump(data, f)
    
    def get_disease(self, disease_id: str) -> Optional[dict]:
        """
        Get disease information by ID.
        
        Args:
            disease_id: OMIM ID (e.g., "OMIM:219700") or Orphanet ID
            
        Returns:
            Disease data dict or None
        """
        # Normalize ID format
        if disease_id.isdigit():
            disease_id = f"OMIM:{disease_id}"
        
        # Check cache
        cached = self._get_cached(f"disease_{disease_id}")
        if cached:
            return cached
        
        # API call
        url = f"{MONARCH_BASE}/entity/{disease_id}"
        
        try:
            response = self.session.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            # Cache result
            self._set_cached(f"disease_{disease_id}", data)
            return data
            
        except requests.RequestException as e:
            print(f"API error for {disease_id}: {e}")
            return None
    
    def search_diseases(self, query: str, limit: int = 10) -> List[dict]:
        """
        Search for diseases by name.
        
        Args:
            query: Search string (e.g., "cystic fibrosis")
            limit: Max results
            
        Returns:
            List of disease matches
        """
        # Check cache
        cache_key = f"search_{query.replace(' ', '_')}_{limit}"
        cached = self._get_cached(cache_key)
        if cached:
            return cached.get('results', [])
        
        url = f"{MONARCH_BASE}/search"
        params = {
            "q": query,
            "category": "biolink:Disease",
            "limit": limit
        }
        
        try:
            response = self.session.get(url, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            # Extract results
            results = data.get('items', [])
            
            # Cache
            self._set_cached(cache_key, {'results': results})
            return results
            
        except requests.RequestException as e:
            print(f"Search error: {e}")
            return []
    
    def get_age_of_onset(self, disease_id: str) -> List[str]:
        """
        Extract age of onset HPO terms for a disease.
        
        Returns:
            List of HPO codes for age of onset
        """
        # Get disease entity (includes has_phenotype)
        disease = self.get_disease(disease_id)
        if not disease:
            return []
        
        # HPO age of onset terms
        age_onset_terms = {
            "HP:0030674": "Antenatal onset",
            "HP:0003577": "Congenital onset",
            "HP:0003623": "Neonatal onset",
            "HP:0003593": "Infantile onset",
            "HP:0011463": "Childhood onset",
            "HP:0003621": "Juvenile onset",
            "HP:0011462": "Young adult onset",
            "HP:0003581": "Adult onset",
            "HP:0003584": "Late onset",
        }
        
        found_ages = []
        
        # Check has_phenotype list from disease entity
        phenotype_ids = disease.get('has_phenotype', [])
        if phenotype_ids:
            for pheno_id in phenotype_ids:
                if pheno_id in age_onset_terms:
                    found_ages.append(pheno_id)
        
        return found_ages
    
    def get_inheritance(self, disease_id: str) -> List[str]:
        """
        Get inheritance patterns for a disease.
        
        Returns:
            List like ["autosomal_recessive", "autosomal_dominant"]
        """
        disease = self.get_disease(disease_id)
        if not disease:
            return []
        
        patterns = []
        
        # Check inheritance field from Monarch
        inheritance_info = disease.get('inheritance')
        if inheritance_info:
            inh_name = inheritance_info.get('name', '').lower()
            if 'autosomal recessive' in inh_name:
                patterns.append('autosomal_recessive')
            elif 'autosomal dominant' in inh_name:
                patterns.append('autosomal_dominant')
            elif 'x-linked recessive' in inh_name:
                patterns.append('x_linked_recessive')
            elif 'x-linked dominant' in inh_name:
                patterns.append('x_linked_dominant')
            elif 'mitochondrial' in inh_name:
                patterns.append('mitochondrial')
        
        # Also check description as fallback
        if not patterns:
            description = disease.get('description', '').lower()
            name = disease.get('name', '').lower()
            full_text = description + " " + name
            
            if 'autosomal recessive' in full_text:
                patterns.append('autosomal_recessive')
            if 'autosomal dominant' in full_text:
                patterns.append('autosomal_dominant')
            if 'x-linked recessive' in full_text:
                patterns.append('x_linked_recessive')
            if 'x-linked dominant' in full_text:
                patterns.append('x_linked_dominant')
            if 'mitochondrial' in full_text:
                patterns.append('mitochondrial')
        
        return patterns


# Singleton
_client = None

def get_client() -> MonarchClient:
    """Get singleton Monarch client."""
    global _client
    if _client is None:
        _client = MonarchClient()
    return _client
