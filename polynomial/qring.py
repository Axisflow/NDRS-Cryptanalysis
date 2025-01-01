from typing import *
from random import randint
from sympy import symbols, Poly, div, GF, invert

class QRPoly(list[int]):
    def __init__(self, degree: int, modulus: int, coeffs: list[int] = [0]):
        self.n = degree
        self.p = modulus
        self.x = symbols('x')
        self.field = GF(self.p)
        self.modulus = Poly(self.x ** self.n + 1, self.x, modulus=self.p)
        self.poly = Poly(coeffs, self.x, modulus=self.p)
        self.__set_in_ring()

    def __set_in_ring(self):
        _, self.poly = div(self.poly, self.modulus, domain=self.field)
        self.clear()
        self.extend(self.poly.all_coeffs())

    def set_module(self, modulus: int):
        self.p = modulus
        self.field = GF(self.p)
        self.modulus = Poly(self.x ** self.n + 1, self.x, modulus=self.p)
        self.__set_in_ring()

    def set_degree(self, degree: int):
        self.n = degree
        self.modulus = Poly(self.x ** self.n + 1, self.x, modulus=self.p)
        self.__set_in_ring()

    def divmod(self, other: 'QRPoly') -> tuple['QRPoly', 'QRPoly']:
        q, r = div(self.poly, other.poly, domain=self.field)
        return QRPoly(self.n, self.p, q.all_coeffs()), QRPoly(self.n, self.p, r.all_coeffs())

    def __add__(self, other: 'QRPoly') -> 'QRPoly':
        return QRPoly(self.n, self.p, (self.poly + other.poly).all_coeffs())
    
    def __sub__(self, other: 'QRPoly') -> 'QRPoly':
        return QRPoly(self.n, self.p, (self.poly - other.poly).all_coeffs())
    
    def __mul__(self, other: 'QRPoly') -> 'QRPoly':
        return QRPoly(self.n, self.p, (self.poly * other.poly).all_coeffs())
    
    def __floordiv__(self, other: 'QRPoly') -> 'QRPoly':
        q, _ = self.divmod(other)
        return q
    
    def __mod__(self, other: 'QRPoly') -> 'QRPoly':
        _, r = self.divmod(other)
        return r
    
    def invertible(self) -> bool:
        try:
            invert(self.poly, self.modulus, domain=self.field)
            return True
        except:
            return False
    
    def inverse(self) -> 'QRPoly':
        return QRPoly(self.n, self.p, invert(self.poly, self.modulus, domain=self.field).all_coeffs())
    
    def __pow__(self, other: int) -> 'QRPoly':
        if other == 0:
            return QRPoly(self.n, self.p, [1])
        elif other == 1:
            return QRPoly(self.n, self.p, self.poly.all_coeffs())
        elif other < 0:
            base = self.inverse()
            other = -other
        else:
            base = self

        if other % 2 == 0:
            return (base ** (other // 2)) ** 2
        else:
            return base * (base ** (other - 1))
    
    def __eq__(self, other: 'QRPoly') -> bool:
        return self.poly == other.poly
    
    def __ne__(self, other: 'QRPoly') -> bool:
        return self.poly != other.poly
    
    def __lt__(self, other: int) -> list[bool]:
        return [all(self[i] < other) for i in range(len(self))]
    
    def __le__(self, other: int) -> list[bool]:
        return [all(self[i] <= other) for i in range(len(self))]
    
    def __gt__(self, other: int) -> list[bool]:
        return [all(self[i] > other) for i in range(len(self))]
    
    def __ge__(self, other: int) -> list[bool]:
        return [all(self[i] >= other) for i in range(len(self))]
    
    def __repr__(self):
        return f"QuotientRing(degree={self.n}, modulus={self.p}, poly={self.poly})"
    
    def __str__(self):
        return f"{self.poly}"
    
    @staticmethod
    def random(degree: int, modulus: int) -> 'QRPoly':
        return QRPoly(degree, modulus, [randint(0, modulus - 1) for _ in range(degree)])

    @staticmethod
    def one(degree: int, modulus: int) -> 'QRPoly':
        return QRPoly(degree, modulus, [1])

# test
if __name__ == "__main__":
    a = QRPoly(5, 11, [1, 2, 3, 4, 5])
    b = QRPoly(5, 11, [5, 4, 3, 2, 1])
    
    print(f'a = {a}')
    print(f'b = {b}')
    
    print(f'a + b = {a + b}')
    
    print(f'a - b = {a - b}')
    
    print(f'a * b = {a * b}')
    
    print(f'a // b = {a // b}')
    print(f'a % b = {a % b}')
    
    print(f'b // a = {b // a}')
    print(f'b % a = {b % a}')
    
    a.set_module(3)
    print(f'set mod 1: a = {a}')
    
    b.set_module(3)
    print(f'set mod 1: b = {b}')

    while True:
        c = QRPoly.random(5, 11)
        if c.invertible():
            print(f'c = {c}')
            print(f'c^-1 = {c.inverse()}')
            print(f'c * c^-1 = {c * c.inverse()}')
            print(f'c^-1 * c = {c.inverse() * c}')
            break
