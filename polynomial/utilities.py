import math
from typing import *

T = TypeVar('T')

def is_prime(n: int) -> bool:
        if n == 2:
            return True
        if n < 2 or n % 2 == 0:
            return False
        for i in range(3, math.isqrt(n) + 1, 2):
            if n % i == 0:
                return False
        return True

def lsum(a: List[T], when_empty: T | None = None) -> T:
    if len(a) == 0:
        return when_empty
    else:
        s = a[0]
        for x in a[1:]:
            s = s + x

        return s

def lsumprod(a: List[T], b: List[T], when_empty: T | None = None) -> T:
    if len(a) != len(b):
        raise ValueError(f"Lengths of a:[{a}] and b:[{b}] must be the same.")
    elif len(a) == 0:
        return when_empty
    else:
        s = a[0] * b[0]
        for i in range(1, len(a)):
            s = s + a[i] * b[i]

        return s