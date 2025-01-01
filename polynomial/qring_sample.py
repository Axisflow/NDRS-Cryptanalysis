from typing import *

from .qring import QRPoly
from .utilities import lsumprod

class QRPolySamples(list[QRPoly]):
    def __init__(self, samples: Iterable[QRPoly]):
        super().__init__(samples)

    def hashing(self, another: 'QRPolySamples') -> QRPoly:
        return lsumprod(self, another, QRPoly(1, 1))
    
    def set_module(self, modulus: int):
        for p in self:
            p.set_module(modulus)
    
    def __mul__(self, another: QRPoly) -> 'QRPolySamples':
        return QRPolySamples([self[i] * another for i in range(len(self))])
    
    def __rmul__(self, another: QRPoly) -> 'QRPolySamples':
        return self * another
    
    def __add__(self, another: 'QRPolySamples') -> 'QRPolySamples':
        return QRPolySamples([self[i] + another[i] for i in range(len(self))])
    
    def __sub__(self, another: 'QRPolySamples') -> 'QRPolySamples':
        return QRPolySamples([self[i] - another[i] for i in range(len(self))])
    
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
        return QRPolySamples([QRPoly.random(degree, modulus) for _ in range(num_samples)])
    
    def __repr__(self):
        return f"QRPolySamples({repr(super())})"
    
    def __str__(self):
        return f"{str(super())}"