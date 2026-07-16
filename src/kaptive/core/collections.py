from __future__ import annotations
from typing import Protocol, TypeVar, Iterable, Any
import numpy as np
import numpy.typing as npt

T = TypeVar('T', covariant=True)
S = TypeVar('S', bound='BatchedContainer')

class BatchedContainer(Protocol[T, S]):
    """A formal contract for Structure-of-Array (SoA) and batched containers in Kaptive.
    
    This protocol ensures that all high-performance vectorized collections provide a standard
    interface for instantiation, concatenation, and indexing without relying on slow
    dynamic mixins or runtime type introspection.
    """
    
    def __len__(self) -> int:
        """Returns the number of records in the batch."""
        ...
    
    def __getitem__(self, item: int | slice | npt.NDArray | list) -> T | S:
        """Accesses records by index, slice, or boolean/integer array mask.
        
        Args:
            item: An integer index, slice, or mask/indices array.
            
        Returns:
            T | S: A single scalar record (T) if an integer is provided, or a new
                batched collection (S) if a slice or mask is provided.
        """
        ...
    
    @classmethod
    def empty(cls: type[S]) -> S:
        """Creates an empty, 0-length collection with correctly typed arrays."""
        ...
    
    @classmethod
    def concat(cls: type[S], batches: Iterable[S]) -> S:
        """Concatenates multiple collections into a single, larger collection.
        
        Args:
            batches: An iterable of collections of the same type.
            
        Returns:
            S: A new, combined collection.
        """
        ...


class RaggedArrayContainer(BatchedContainer[T, S], Protocol[T, S]):
    """A formal contract for ragged SoA containers.
    
    Ragged containers store variable-length data (like sequences or CIGARs) in a flat,
    contiguous memory layout, and manage boundaries using `offsets` and `lengths`.
    """
    
    offsets: npt.NDArray[np.int32]
    lengths: npt.NDArray[np.int32]
