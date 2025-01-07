from typing import *

from .qring import QRPoly
from .utilities import lsumprod

class QRPolySamples(list[QRPoly]):
    def __init__(self, degree: int, modulus: int, samples: Iterable[QRPoly]):
        self.n = degree
        self.p = modulus
        super().__init__(samples)

    def copy(self, degree: int = None, modulus: int = None) -> 'QRPolySamples':
        return QRPolySamples(degree or self.n, modulus or self.p, [poly.copy() for poly in self])

    def hashing(self, another: 'QRPolySamples') -> QRPoly:
        return lsumprod(self, another, QRPoly(self.n, self.p))
    
    def set_degree(self, degree: int):
        self.n = degree
        for poly in self:
            poly.set_degree(degree)

    def set_module(self, modulus: int):
        self.p = modulus
        for poly in self:
            poly.set_module(modulus)
    
    def __mul__(self, another: QRPoly) -> 'QRPolySamples':
        return QRPolySamples(self.n, self.p, [self[i] * another for i in range(len(self))])
    
    def __rmul__(self, another: QRPoly) -> 'QRPolySamples':
        return self * another
    
    def __floordiv__(self, another: QRPoly) -> 'QRPolySamples':
        return QRPolySamples(self.n, self.p, [self[i] // another for i in range(len(self))])
    
    def __add__(self, another: 'QRPolySamples') -> 'QRPolySamples':
        return QRPolySamples(self.n, self.p, [self[i] + another[i] for i in range(len(self))])
    
    def __sub__(self, another: 'QRPolySamples') -> 'QRPolySamples':
        return QRPolySamples(self.n, self.p, [self[i] - another[i] for i in range(len(self))])
    
    def __lt__(self, another: int) -> list[bool]:
        return [all(self[i] < another) for i in range(len(self))]
    
    def __le__(self, another: int) -> list[bool]:
        return [all(self[i] <= another) for i in range(len(self))]
    
    def __gt__(self, another: int) -> list[bool]:
        return [all(self[i] > another) for i in range(len(self))]
    
    def __ge__(self, another: int) -> list[bool]:
        return [all(self[i] >= another) for i in range(len(self))]

    @staticmethod
    def random(degree: int, modulus: int, num_samples: int) -> 'QRPolySamples':
        return QRPolySamples(degree, modulus, [QRPoly.random(degree, modulus) for _ in range(num_samples)])
    
    def __repr__(self):
        return f"QRPolySamples({super().__repr__()}, deg={self.n}, mod={self.p})"
    
    def __str__(self):
        return f"{super().__str__()}"