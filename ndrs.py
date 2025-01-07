from typing import *
import math
import random as rnd

from polynomial.qring_sample import QRPolySamples
from polynomial.qring import QRPoly
from polynomial.utilities import is_prime, lsum

class Key(QRPolySamples):
    def __init__(self, qrpoly_samples: QRPolySamples):
        super().__init__(qrpoly_samples.n, qrpoly_samples.p, qrpoly_samples)
    def __repr__(self):
        return f"Key({super().__repr__()})"

class SecretKey(Key):
    def __repr__(self):
        return f"SecretKey({super().__repr__()})"

class PublicKey(Key):
    def __repr__(self):
        return f"PublicKey({super().__repr__()})"

class KeyPair(NamedTuple):
    pk: PublicKey
    sk: SecretKey
    def __repr__(self):
        return f"KeyPair(pk={self.pk}, sk={self.sk})"

class Signature(NamedTuple):
    pks: list[PublicKey]
    b_hat: QRPolySamples
    A: QRPoly
    z_hats: list[QRPolySamples]
    vs: list[QRPoly]
    def __repr__(self):
        return f"Signature(pks={self.pks}, b_hat={self.b_hat}, A={self.A}, z_hats={self.z_hats}, vs={self.vs})"

class Evidence(NamedTuple):
    sigmai: QRPoly
    alphai: QRPoly
    betai: QRPoly
    zi_hat: QRPolySamples
    ei: QRPoly
    def __repr__(self):
        return f"Evidence(sigmai={self.sigmai}, alphai={self.alphai}, betai={self.betai}, zi_hat={self.zi_hat}, ei={self.ei})"

class NDRS:
    Ds_max: int = 1
    Ds_mod: int = Ds_max * 2 + 1

    def __init__(self, k: int, c: int = 3):
        self.k = k # security parameter
        self.c = c # A constant for robustness
        self.n = int(2 ** (math.floor(math.log2(k)) + 1))
        self.m = round((3 + 2 * c / 3) * math.log2(self.n))
        self.__sqrtnlogn = int(math.sqrt(self.n) * math.log2(self.n))
        self.Dy_max = int(self.m * self.n**1.5 * math.log2(self.n))
        self.Dh_max = self.Dy_max + self.__sqrtnlogn
        self.Dz_max = self.Dy_max - self.__sqrtnlogn
        self.Dy_mod = self.Dy_max * 2 + 1
        self.Dh_mod = self.Dh_max * 2 + 1
        self.Dz_mod = self.Dz_max * 2 + 1
        self.D_mod = self.find_p()
        self.S = QRPoly(self.n, self.D_mod)
        while self.S == QRPoly(self.n, self.D_mod):
            self.S = QRPoly.random(self.n, self.D_mod)

    def __repr__(self) -> str:
        return f"P(k={self.k}, c={self.c}, n={self.n}, m={self.m}, D_mod={self.D_mod}, Dy_mod={self.Dy_mod}, Dh_mod={self.Dh_mod}, Dz_mod={self.Dz_mod}, Ds_mod={self.Ds_mod}, S={self.S})"

    def find_p(self) -> int:
        p = self.n ** (4 + self.c)
        while True:
            if is_prime(p) and p % 8 == 3:
                return p
            p += 1

    def hash1(self, idx: int, pk: PublicKey) -> QRPoly:
        # concatenate idx and elements of pk
        cancatenated = []
        for i in range(self.m):
            _p = list(pk[i].copy())
            _p.extend([0] * (self.n - len(_p)))
            cancatenated.extend(_p)

        coeffs = []
        for i in range(self.n):
            coeffs.append(hash(tuple([idx] + cancatenated[i*self.m:(i+1)*self.m])))
            
        return QRPoly(self.n, self.Ds_mod, coeffs)

    def hash2(self, sum_alphas: QRPoly, betas: List[QRPoly], A: QRPoly, pks: List[PublicKey], msg: str) -> QRPoly:
        return self.hash1(hash(msg), lsum(pks) * (sum_alphas + lsum(betas) + A))

    def hash3(self, alphai: QRPoly, betai: QRPoly, A: QRPoly, pks: List[PublicKey], msg: str) -> QRPoly:
        return self.hash1(hash(msg), (alphai + betai + A) * lsum(pks))

    def key_gen(self) -> KeyPair:
        st0_idx = -1
        while st0_idx == -1:
            si_hat = QRPolySamples.random(self.n, self.Ds_mod, self.m)
            for i in range(self.m):
                if si_hat[i].invertible():
                    st0_idx = i
                    break
        
        ai_hat = QRPolySamples.random(self.n, self.D_mod, self.m)
        ai_hat[st0_idx] = (self.S - lsum([ai_hat[i] * si_hat[i] for i in range(self.m) if i != st0_idx])) * si_hat[st0_idx].inverse()
        return KeyPair(PublicKey(ai_hat), SecretKey(si_hat))
    
    def sign(self, signer: KeyPair, other_signers: List[PublicKey], msg: str) -> Signature:
        signer_idx = rnd.randint(0, len(other_signers))
        pks = other_signers[0:signer_idx] + [signer.pk] + other_signers[signer_idx:]
        signers_idices = list(range(len(pks)))
        signers_idices.remove(signer_idx)
        
        sigmaj = QRPoly(self.n, self.Ds_mod)
        while sigmaj == QRPoly(self.n, self.Ds_mod):
            b_hat = QRPolySamples.random(self.n, self.D_mod, self.m)
            sigmaj = b_hat.hashing(signer.sk)

        A = sigmaj - self.S * self.hash1(signer_idx, signer.pk) 

        while True:
            yj_hat = QRPolySamples.random(self.n, self.Dy_mod, self.m)
            alphaj = signer.pk.hashing(yj_hat)
            betaj = b_hat.hashing(yj_hat)

            sum_alphas = alphaj.copy()
            betas = list[QRPoly]()
            vs = list[QRPoly]()
            z_hats = list[QRPolySamples]()
            for i in signers_idices:
                z_hats.append(QRPolySamples.random(self.n, self.Dz_mod, self.m))
                vs.append(QRPoly.random(self.n, self.Ds_mod))
                sigmai = self.S * self.hash1(i, pks[i]) + A
                sum_alphas = sum_alphas + pks[i].hashing(z_hats[-1]) - self.S * vs[-1]
                betas.append(b_hat.hashing(z_hats[-1]) - sigmai * vs[-1])

            betas.insert(signer_idx, betaj)
            v = self.hash2(sum_alphas, betas, A, pks, msg)
            vs.insert(signer_idx, v - lsum(vs))
            z_hats.insert(signer_idx, yj_hat + signer.sk * vs[signer_idx])
            
            if all(-self.Dz_max <= z_hats[signer_idx] <= self.Dz_max) and all(-self.Ds_max <= vs[signer_idx] <= self.Ds_max):
                break
        
        return Signature(pks, b_hat, A, z_hats, vs)
    
    def verify(self, msg: str, sig: Signature) -> bool:
        sigma_prime = list[QRPoly]()
        alpha_prime = list[QRPoly]()
        beta_prime = list[QRPoly]()
        for i in range(len(sig.pks)):
            sigma_prime.append(self.S * self.hash1(i, sig.pks[i]) + sig.A)
            alpha_prime.append(sig.pks[i].hashing(sig.z_hats[i]) - self.S * sig.vs[i])
            beta_prime.append(sig.b_hat.hashing(sig.z_hats[i]) - sigma_prime[i] * sig.vs[i])

        v_prime = self.hash2(lsum(alpha_prime), beta_prime, sig.A, sig.pks, msg)
        return v_prime == lsum(sig.vs)
    
    def evidence_gen(self, signer: KeyPair, msg: str, sig: Signature) -> QRPoly:
        if not self.verify(msg, sig):
            raise ValueError("Invalid signature or message.")
        
        sigmai = sig.b_hat.hashing(signer.sk)
        yi_hat = QRPolySamples.random(self.n, self.Dy_mod, self.m)
        alphai = signer.pk.hashing(yi_hat)
        betai = sig.b_hat.hashing(yi_hat)
        ei = self.hash3(alphai, betai, sig.A, sig.pks, msg)
        zi_hat = yi_hat + signer.sk * ei

        return Evidence(sigmai, alphai, betai, zi_hat, ei)
    
    def evidence_check(self, singer: KeyPair, msg: str, sig: Signature, evi: Evidence) -> bool:
        if not self.verify(msg, sig):
            raise ValueError("Invalid signature or message.")
        
        alphai_prime = singer.pk.hashing(evi.zi_hat) - self.S * evi.ei
        betai_prime = sig.b_hat.hashing(evi.zi_hat) - evi.sigmai * evi.ei
        ei_prime = self.hash3(alphai_prime, betai_prime, sig.A, sig.pks, msg)
        if ei_prime != evi.ei:
            raise ValueError("Invalid evidence.")
        
        return evi.sigmai == self.S * self.hash1(sig.pks.index(singer.pk), singer.pk) + sig.A

class Frameable_NDRS(NDRS):
    def __init__(self, k: int, c: int = 3):
        super().__init__(k, c)

    def key_gen(self) -> KeyPair:
        key = super().key_gen()
        # find a element of secret key that is invertible and make it be the first element of the secret key
        st0_idx = -1
        while st0_idx == -1:
            for i in range(self.m):
                if key.sk[i].invertible():
                    st0_idx = i
                    break

        key.sk[0], key.sk[st0_idx] = key.sk[st0_idx], key.sk[0]
        key.pk[0], key.pk[st0_idx] = key.pk[st0_idx], key.pk[0]
        return key
    
    def fake_skey_gen(self, original_key: KeyPair) -> Tuple[KeyPair, SecretKey]:
        while True:
            sk = QRPolySamples.random(self.n, self.Ds_mod, self.m)
            if sk[1].invertible():
                break
        
        sk[0] = QRPoly(self.n, self.Ds_mod)
        pk = original_key.pk.copy()
        pk[1] = (self.S - lsum([pk[i] * sk[i] for i in range(2, self.m)])) * sk[1].inverse()
        pk[0] = (self.S - lsum([pk[i] * sk[i] for i in range(1, self.m)])) * original_key.sk[0].inverse()
        
        return KeyPair(pk, original_key.sk), sk
    
    def frameably_sign(self, signer: KeyPair, other_signers: List[PublicKey], msg: str, framed_idx: int) -> Signature:
        signer_idx = rnd.randint(0, len(other_signers))
        pks = other_signers[0:signer_idx] + [signer.pk] + other_signers[signer_idx:]
        signers_idices = list(range(len(pks)))
        signers_idices.remove(signer_idx)
        
        b_hat = other_signers[framed_idx] * (self.hash1(framed_idx, other_signers[framed_idx]) - self.hash1(signer_idx, signer.pk)) * self.S / (self.S - other_signers[framed_idx].hashing(signer.sk))
        sigmaj = b_hat.hashing(signer.sk)

        A = sigmaj - self.hash1(signer_idx, signer.pk) * self.S

        yj_hat = QRPolySamples.random(self.n, self.Dy_mod, self.m)
        alphaj = signer.pk.hashing(yj_hat)
        betaj = b_hat.hashing(yj_hat)

        sum_alphas = alphaj.copy()
        betas = list[QRPoly]()
        vs = list[QRPoly]()
        z_hats = list[QRPolySamples]()
        for i in signers_idices:
            z_hats.append(QRPolySamples.random(self.n, self.Dz_mod, self.m))
            vs.append(QRPoly.random(self.n, self.Ds_mod))
            sigmai = self.S * self.hash1(i, pks[i]) + A
            sum_alphas += pks[i].hashing(z_hats[-1]) - self.S * vs[-1]
            betas.append(b_hat.hashing(z_hats[-1]) - sigmai * vs[-1])

        betas.insert(signer_idx, betaj)
        v = self.hash2(sum_alphas, betas, A, pks, msg)
        vs.insert(signer_idx, v - lsum(vs))
        z_hats.insert(signer_idx, yj_hat + signer.sk * vs[signer_idx])
        
        return Signature(pks, b_hat, A, z_hats, vs)
