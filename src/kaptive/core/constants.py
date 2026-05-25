"""
Constants shared across modules
"""
from enum import StrEnum, IntEnum, Flag, auto


# String Enums ---------------------------------------------------------------------------------------------------------
class FeatureType(StrEnum):
    """Genomic feature types commonly encountered in genome annotation."""
    CDS = "CDS"
    MOBILE_ELEMENT = "mobile_element"
    REGULATORY = "regulatory"
    REPEAT_REGION = "repeat_region"


class Orientation(StrEnum):
    """The relative strand orientation between two genomic features."""
    SAME = "same strand"
    OPPOSITE = "opposite strand"
    NONE = "-"


class Effect(Flag):
    """
    Biological impact of a structural variant on a genomic feature.
    
    Can be combined using bitwise OR (|) to represent multiple concurrent effects.
    """
    NONE = auto()
    UPREGULATED = auto()
    TRUNCATED = auto()
    DISRUPTED = auto()


# Int Enums ------------------------------------------------------------------------------------------------------------
class Context(IntEnum):
    """
    Spatial relationship between two genomic intervals.

    Used to describe the position of a 'passenger' or 'flank' gene relative
    to a primary target (e.g. a mobile element insertion site).
    """
    UPSTREAM = auto()
    DOWNSTREAM = auto()
    INSIDE = auto()
    OVERLAPPING = auto()
    OVERLAPPING_START = auto()
    OVERLAPPING_END = auto()


class Strand(IntEnum):
    """
    Integer representation of genomic strand orientation.

    Supports conversion from common string formats (+, -, 1, -1).
    """
    FORWARD = 1
    REVERSE = -1
    UNSTRANDED = 0

    @classmethod
    def _missing_(cls, value):
        if isinstance(value, bytes):
            value = value.decode('ascii')
        if isinstance(value, str):
            if value == '+' or value == '1' or value == '+1':
                return Strand.FORWARD
            elif value == '-' or value == '-1':
                return Strand.REVERSE
        return Strand.UNSTRANDED

    def __str__(self):
        if self == Strand.FORWARD: return '+'
        if self == Strand.REVERSE: return '-'
        return '.'