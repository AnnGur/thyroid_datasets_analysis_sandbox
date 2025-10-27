from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class AnalysisConfig:
    output_dir: str = 'results'
    comparisons: List[Tuple[str, str]] = None
    
    def __post_init__(self):
        if self.comparisons is None:
            self.comparisons = [
                ('PTC_BRAF_mut', 'PTC_BRAF_wt'),
                ('PTC_BRAF_mut', 'Normal'),
                ('FA', 'FVPTC')
            ]