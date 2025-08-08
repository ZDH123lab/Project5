import numpy as np
import random
import hashlib
import time
from typing import Tuple, List

# SM2椭圆曲线参数
SM2_p = 0x8542D69E4C044F18E8B92435BF6FF7DE457283915C45517D722EDB8B08F1DFC3
SM2_a = 0x787968B4FA32C3FD2417842E73BBFEFF2F3C848B6831D7E0EC65228B3937E498
SM2_b = 0x63E4C6D3B23B0C849CF84241484BFE48F61D59A5B16BA06E6E12D1DA27C5249A
SM2_n = 0x8542D69E4C044F18E8B92435BF6FF7DD297720630485628D5AE74EE7C32E79B7
SM2_Gx = 0x421DEBD61B62EAB6746434EBC3CC315E32220B3BADD50BDC4C4E6C147FEDD43D
SM2_Gy = 0x0680512BCBB42C07D47349D2153B70C4E5D7FDFCBFA36EA1A85841B9E46E09A2

# 预计算参数
SM2_h = 1  # 辅因子


class ECPoint:
    """椭圆曲线点类（雅可比坐标）"""

    def __init__(self, x: int, y: int, z: int = 1):
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other) -> bool:
        if not isinstance(other, ECPoint):
            return False
        x1 = (self.x * pow(other.z, 2, SM2_p)) % SM2_p
        x2 = (other.x * pow(self.z, 2, SM2_p)) % SM2_p
        y1 = (self.y * pow(other.z, 3, SM2_p)) % SM2_p
        y2 = (other.y * pow(self.z, 3, SM2_p)) % SM2_p
        return x1 == x2 and y1 == y2

    def is_infinity(self) -> bool:
        """判断是否是无穷远点"""
        return self.z == 0

    def to_affine(self) -> 'ECPoint':
        """转换为仿射坐标"""
        if self.is_infinity():
            return ECPoint(0, 0, 0)
        z_inv = pow(self.z, SM2_p - 2, SM2_p)
        z_inv_sqr = (z_inv * z_inv) % SM2_p
        x = (self.x * z_inv_sqr) % SM2_p
        y = (self.y * z_inv_sqr * z_inv) % SM2_p
        return ECPoint(x, y, 1)

    def __str__(self) -> str:
        if self.is_infinity():
            return "Infinity"
        affine = self.to_affine()
        return f"({hex(affine.x)}, {hex(affine.y)})"


# 生成元点G
G = ECPoint(SM2_Gx, SM2_Gy)


def mod_inv(a: int, p: int) -> int:
    """扩展欧几里得算法求模逆"""
    old_r, r = a, p
    old_s, s = 1, 0
    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
    return old_s % p


def point_double(P: ECPoint) -> ECPoint:
    """椭圆曲线点倍运算 (雅可比坐标)"""
    if P.is_infinity():
        return P

    # 点倍公式 (4M + 4S + 10A)
    X, Y, Z = P.x, P.y, P.z
    XX = (X * X) % SM2_p
    YY = (Y * Y) % SM2_p
    YYYY = (YY * YY) % SM2_p
    ZZ = (Z * Z) % SM2_p

    S = (4 * X * YY) % SM2_p
    M = (3 * XX + SM2_a * ZZ * ZZ) % SM2_p

    X3 = (M * M - 2 * S) % SM2_p
    Y3 = (M * (S - X3) - 8 * YYYY) % SM2_p
    Z3 = (2 * Y * Z) % SM2_p

    return ECPoint(X3, Y3, Z3)


def point_add(P: ECPoint, Q: ECPoint) -> ECPoint:
    """椭圆曲线点加运算 (雅可比坐标)"""
    if P.is_infinity():
        return Q
    if Q.is_infinity():
        return P

    # 点加公式 (12M + 4S + 7A)
    X1, Y1, Z1 = P.x, P.y, P.z
    X2, Y2, Z2 = Q.x, Q.y, Q.z

    Z1Z1 = (Z1 * Z1) % SM2_p
    Z2Z2 = (Z2 * Z2) % SM2_p
    U1 = (X1 * Z2Z2) % SM2_p
    U2 = (X2 * Z1Z1) % SM2_p
    S1 = (Y1 * Z2 * Z2Z2) % SM2_p
    S2 = (Y2 * Z1 * Z1Z1) % SM2_p

    if U1 == U2:
        if S1 != S2:
            return ECPoint(0, 0, 0)  # 无穷远点
        else:
            return point_double(P)

    H = (U2 - U1) % SM2_p
    HH = (H * H) % SM2_p
    HHH = (H * HH) % SM2_p
    r = (S2 - S1) % SM2_p

    X3 = (r * r - HHH - 2 * U1 * HH) % SM2_p
    Y3 = (r * (U1 * HH - X3) - S1 * HHH) % SM2_p
    Z3 = (Z1 * Z2 * H) % SM2_p

    return ECPoint(X3, Y3, Z3)


def point_neg(P: ECPoint) -> ECPoint:
    """椭圆曲线点取负"""
    if P.is_infinity():
        return P
    return ECPoint(P.x, (-P.y) % SM2_p, P.z)


def point_mul_standard(k: int, P: ECPoint) -> ECPoint:
    """标准实现的标量乘法"""
    k = k % SM2_n
    result = ECPoint(0, 0, 0)  # 无穷远点
    current = P

    while k > 0:
        if k & 1:
            result = point_add(result, current)
        current = point_double(current)
        k >>= 1

    return result


class SIMDOptimizer:
    """SIMD优化实现 """

    def __init__(self):
        self.p = SM2_p
        self.zero = 0

    def mod_add_simd(self, a: int, b: int) -> int:
        """模加法"""
        return (a + b) % self.p

    def mod_sub_simd(self, a: int, b: int) -> int:
        """模减法"""
        return (a - b) % self.p

    def mod_mul_simd(self, a: int, b: int) -> int:
        """模乘法"""
        return (a * b) % self.p

    def point_double_simd(self, P: ECPoint) -> ECPoint:
        """使用优化的点倍运算"""
        if P.is_infinity():
            return P

        X, Y, Z = P.x, P.y, P.z
        XX = self.mod_mul_simd(X, X)
        YY = self.mod_mul_simd(Y, Y)
        YYYY = self.mod_mul_simd(YY, YY)
        ZZ = self.mod_mul_simd(Z, Z)

        S = self.mod_mul_simd(4, self.mod_mul_simd(X, YY))
        M = self.mod_add_simd(
            self.mod_mul_simd(3, XX),
            self.mod_mul_simd(
                self.mod_mul_simd(SM2_a, ZZ),
                ZZ
            )
        )

        X3 = self.mod_sub_simd(
            self.mod_mul_simd(M, M),
            self.mod_mul_simd(2, S)
        )
        Y3 = self.mod_sub_simd(
            self.mod_mul_simd(M, self.mod_sub_simd(S, X3)),
            self.mod_mul_simd(8, YYYY)
        )
        Z3 = self.mod_mul_simd(
            self.mod_mul_simd(2, Y),
            Z
        )

        return ECPoint(X3, Y3, Z3)

    def point_add_simd(self, P: ECPoint, Q: ECPoint) -> ECPoint:
        """使用优化的点加运算"""
        if P.is_infinity():
            return Q
        if Q.is_infinity():
            return P

        X1, Y1, Z1 = P.x, P.y, P.z
        X2, Y2, Z2 = Q.x, Q.y, Q.z

        Z1Z1 = self.mod_mul_simd(Z1, Z1)
        Z2Z2 = self.mod_mul_simd(Z2, Z2)
        U1 = self.mod_mul_simd(X1, Z2Z2)
        U2 = self.mod_mul_simd(X2, Z1Z1)
        S1 = self.mod_mul_simd(Y1, self.mod_mul_simd(Z2, Z2Z2))
        S2 = self.mod_mul_simd(Y2, self.mod_mul_simd(Z1, Z1Z1))

        if U1 == U2:
            if S1 != S2:
                return ECPoint(0, 0, 0)  # 无穷远点
            else:
                return self.point_double_simd(P)

        H = self.mod_sub_simd(U2, U1)
        HH = self.mod_mul_simd(H, H)
        HHH = self.mod_mul_simd(H, HH)
        r = self.mod_sub_simd(S2, S1)

        X3 = self.mod_sub_simd(
            self.mod_sub_simd(self.mod_mul_simd(r, r), HHH),
            self.mod_mul_simd(
                self.mod_mul_simd(2, U1),
                HH
            )
        )
        Y3 = self.mod_sub_simd(
            self.mod_mul_simd(
                r,
                self.mod_sub_simd(
                    self.mod_mul_simd(U1, HH),
                    X3
                )
            ),
            self.mod_mul_simd(S1, HHH)
        )
        Z3 = self.mod_mul_simd(
            self.mod_mul_simd(Z1, Z2),
            H
        )

        return ECPoint(X3, Y3, Z3)

    def scalar_mul_simd(self, k: int, P: ECPoint) -> ECPoint:
        """使用优化的标量乘法"""
        k = k % SM2_n
        result = ECPoint(0, 0, 0)  # 无穷远点
        current = P

        while k > 0:
            if k & 1:
                result = self.point_add_simd(result, current)
            current = self.point_double_simd(current)
            k >>= 1

        return result


def bytes_to_int(b: bytes) -> int:
    """将字节串转换为整数"""
    return int.from_bytes(b, byteorder='big')


def int_to_bytes(x: int, length: int = None) -> bytes:
    """将整数转换为字节串"""
    if length is None:
        length = (x.bit_length() + 7) // 8 or 1
    return x.to_bytes(length, byteorder='big')


def hash_msg(msg: bytes) -> int:
    """计算消息的哈希值"""
    h = hashlib.sha256(msg).digest()
    return bytes_to_int(h)


def kdf(z: bytes, klen: int) -> bytes:
    """密钥派生函数"""
    v = 256 // 8
    ct = 0x00000001
    ha = b''
    for i in range((klen + v - 1) // v):
        ct_bytes = ct.to_bytes(4, byteorder='big')
        ha += hashlib.sha256(z + ct_bytes).digest()
        ct += 1
    return ha[:klen]


def encrypt_sm2(public_key: ECPoint, msg: str) -> bytes:
    """使用SM2加密消息"""
    # 1. 将消息转换为字节
    msg_bytes = msg.encode('utf-8')
    msg_len = len(msg_bytes)

    # 2. 生成随机数k
    while True:
        k = random.randint(1, SM2_n - 1)

        # 3. 计算C1 = [k]G
        simd = SIMDOptimizer()
        C1 = simd.scalar_mul_simd(k, G)
        C1_bytes = int_to_bytes(C1.to_affine().x) + int_to_bytes(C1.to_affine().y)

        # 4. 计算椭圆曲线点[k]P = (x2, y2)
        kP = simd.scalar_mul_simd(k, public_key)
        x2 = kP.to_affine().x
        y2 = kP.to_affine().y
        x2_bytes = int_to_bytes(x2)
        y2_bytes = int_to_bytes(y2)

        # 5. 计算t = KDF(x2 || y2, msg_len)
        t = kdf(x2_bytes + y2_bytes, msg_len)
        if any(t):  # t不全为0
            break

    # 6. 计算C2 = M ⊕ t
    C2 = bytes(a ^ b for a, b in zip(msg_bytes, t))

    # 7. 计算C3 = Hash(x2 || M || y2)
    C3 = hash_msg(x2_bytes + msg_bytes + y2_bytes)
    C3_bytes = int_to_bytes(C3, 32)  # 假设哈希输出256位

    # 8. 输出密文C = C1 || C3 || C2
    return C1_bytes + C3_bytes + C2


def decrypt_sm2(private_key: int, ciphertext: bytes) -> str:
    """使用SM2解密消息"""
    # 1. 从密文中提取C1, C3, C2
    point_len = (SM2_p.bit_length() + 7) // 8
    C1_bytes = ciphertext[:2 * point_len]
    C3_bytes = ciphertext[2 * point_len:2 * point_len + 32]  # 假设哈希输出256位
    C2_bytes = ciphertext[2 * point_len + 32:]
    msg_len = len(C2_bytes)

    # 2. 从C1中恢复椭圆曲线点
    x1 = bytes_to_int(C1_bytes[:point_len])
    y1 = bytes_to_int(C1_bytes[point_len:])
    C1 = ECPoint(x1, y1)

    # 3. 计算椭圆曲线点[dB]C1 = (x2, y2)
    simd = SIMDOptimizer()
    x2y2 = simd.scalar_mul_simd(private_key, C1)
    x2 = x2y2.to_affine().x
    y2 = x2y2.to_affine().y
    x2_bytes = int_to_bytes(x2)
    y2_bytes = int_to_bytes(y2)

    # 4. 计算t = KDF(x2 || y2, msg_len)
    t = kdf(x2_bytes + y2_bytes, msg_len)
    if not any(t):  # t全为0
        raise ValueError("KDF failed - t is all zero")

    # 5. 计算M = C2 ⊕ t
    msg_bytes = bytes(a ^ b for a, b in zip(C2_bytes, t))

    # 6. 计算u = Hash(x2 || M || y2)
    u = hash_msg(x2_bytes + msg_bytes + y2_bytes)
    u_bytes = int_to_bytes(u, 32)

    # 7. 验证u == C3
    if u_bytes != C3_bytes:
        raise ValueError("Hash verification failed")

    # 8. 返回明文
    return msg_bytes.decode('utf-8')


def test_sm2_encryption():
    print("\n=== SM2加密测试 ===")

    # 生成密钥对
    private_key = random.randint(1, SM2_n - 1)
    simd = SIMDOptimizer()
    public_key = simd.scalar_mul_simd(private_key, G)

    print(f"私钥: {hex(private_key)}")
    print(f"公钥: {public_key}")

    # 要加密的消息
    msg = "zdh202200460110"
    print(f"\n原始消息: {msg}")

    # 加密
    start = time.time()
    ciphertext = encrypt_sm2(public_key, msg)
    encrypt_time = time.time() - start
    print(f"\n加密结果 (十六进制): {ciphertext.hex()}")
    print(f"加密耗时: {encrypt_time:.6f}s")

    # 解密
    start = time.time()
    decrypted_msg = decrypt_sm2(private_key, ciphertext)
    decrypt_time = time.time() - start
    print(f"\n解密结果: {decrypted_msg}")
    print(f"解密耗时: {decrypt_time:.6f}s")

    # 验证
    print(f"\n解密是否成功: {decrypted_msg == msg}")


if __name__ == "__main__":
    test_sm2_encryption()