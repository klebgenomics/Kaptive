r"""Protocol contracts for Structure-of-Arrays (SoA) and batched containers.

This module defines formal protocols ([`BatchedContainer`][kaptive.core.collections.BatchedContainer]
and [`RaggedArrayContainer`][kaptive.core.collections.RaggedArrayContainer]) that govern high-performance,
vectorized containers across Kaptive.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Protocol, Self, TypeVar

import numpy as np
import numpy.typing as npt

T = TypeVar("T", covariant=True)
S = TypeVar("S", bound="BatchedContainer[Any, Any]")  # type: ignore


class BatchedContainer(Protocol[T, S]):
    r"""A formal contract for Structure-of-Arrays (SoA) and batched containers in Kaptive.

    This protocol ensures that all high-performance vectorized collections provide a standard
    interface for instantiation, concatenation, and indexing without relying on slow
    dynamic mixins or runtime type introspection.

    Type Parameters:
        T: The scalar record type returned when accessing a single index.
        S: The batched container type returned when indexing with slices or array masks.
    """

    def __len__(self) -> int:
        r"""Return the number of records in the batch.

        Returns:
            int: The total count of items in the container.
        """
        ...

    def __getitem__(self, item: int | slice | npt.NDArray[Any] | list[Any]) -> T | S:  # type: ignore
        r"""Access records by index, slice, or boolean/integer array mask.

        Args:
            item (int | slice | npt.NDArray | list): An integer index, slice, or mask/indices array.

        Returns:
            T | S: A single scalar record (`T`) if an integer index is provided, or a new
                batched collection (`S`) if a slice or mask is provided.

        Raises:
            IndexError: If an integer index is out of bounds.
        """
        ...

    @classmethod
    def empty(cls) -> Self:
        r"""Create an empty, 0-length collection with correctly typed arrays.

        Returns:
            S: An empty instance of the container.
        """
        ...

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:
        r"""Concatenate multiple collections into a single, larger collection.

        Args:
            batches (Iterable[S]): An iterable of collections of the same type.

        Returns:
            S: A new, combined collection.

        Raises:
            ValueError: If the input iterable is empty or contains incompatible batches.
        """
        ...


class RaggedArrayContainer(BatchedContainer[T, S], Protocol[T, S]):
    r"""A formal contract for ragged Structure-of-Arrays (SoA) containers.

    Ragged containers store variable-length data (such as sequences or CIGAR operations) in a flat,
    contiguous memory layout, managing sequence boundaries using `offsets` and `lengths` arrays.

    Attributes:
        offsets (npt.NDArray[np.int32]): Starting indices into the flat data array for each record.
        lengths (npt.NDArray[np.int32]): Number of elements in the flat array for each record.
    """

    offsets: npt.NDArray[np.int32]
    lengths: npt.NDArray[np.int32]
