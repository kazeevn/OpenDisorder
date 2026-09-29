"""opendisorder: is a crystal structure already in the Crystallography Open Database?

Layer 1 (pristine) finds COD entries that are the same structure.
Layer 2 (disorder) finds COD entries with partial occupancies of which the
query is an ordered realization.
"""

from .checker import CheckOptions, NoveltyChecker, NoveltyReport
from .cod import CODClient
from .index import CODIndex
from .matching import (
    DisorderMatcher,
    MatchKind,
    MatchResult,
    PristineMatcher,
    SiblingOrderingMatcher,
)

__all__ = [
    "CheckOptions",
    "CODClient",
    "CODIndex",
    "DisorderMatcher",
    "MatchKind",
    "MatchResult",
    "NoveltyChecker",
    "NoveltyReport",
    "PristineMatcher",
    "SiblingOrderingMatcher",
]

__version__ = "0.1.0"
