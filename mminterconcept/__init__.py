"""
MMInteroperabilityConcept
The concept implementation for the MolSSI MM Interoperability project
"""

# Add imports here
from .mminterconcept import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
