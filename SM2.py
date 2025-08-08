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
    """椭圆曲线点倍运算"""
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
    """椭圆曲线点加运算 """
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


def scalar_mult(k: int, P: ECPoint) -> ECPoint:
    """标量乘法 (使用双倍-加法算法)"""
    result = ECPoint(0, 0, 0)  # 无穷远点
    current = P

    while k > 0:
        if k & 1:
            result = point_add(result, current)
        current = point_double(current)
        k >>= 1

    return result


def bytes_to_int(b: bytes) -> int:
    """将字节串转换为整数"""
    return int.from_bytes(b, 'big')


def int_to_bytes(x: int, length: int = None) -> bytes:
    """将整数转换为字节串"""
    if length is None:
        length = (x.bit_length() + 7) // 8
    return x.to_bytes(length, 'big')


def hash_msg(msg: bytes) -> bytes:
    """SHA-256"""
    return hashlib.sha256(msg).digest()


def kdf(z: bytes, klen: int) -> bytes:
    """密钥派生函数 (KDF)"""
    v = 256 // 8  # SM3输出长度(这里用SHA-256代替)
    ct = 0x00000001
    ha = b''

    for i in range((klen + v - 1) // v):
        ct_bytes = ct.to_bytes(4, 'big')
        ha += hash_msg(z + ct_bytes)
        ct += 1

    return ha[:klen]


def sm2_encrypt(public_key: ECPoint, msg: str) -> Tuple[bytes, float]:
    """SM2加密算法，返回加密结果和耗时(秒)"""
    start_time = time.time()

    # 将消息转换为字节
    msg_bytes = msg.encode('utf-8')
    msg_len = len(msg_bytes)

    # 步骤1: 产生随机数k
    while True:
        k = random.randint(1, SM2_n - 1)
        # 步骤2: 计算椭圆曲线点C1 = [k]G
        C1 = scalar_mult(k, G).to_affine()

        # 步骤3: 计算椭圆曲线点S = [h]Pb
        S = scalar_mult(SM2_h, public_key)
        if S.is_infinity():
            continue  # 如果S是无穷远点，重新选择k

        # 步骤4: 计算椭圆曲线点[k]Pb = (x2, y2)
        kPb = scalar_mult(k, public_key).to_affine()
        x2 = kPb.x
        y2 = kPb.y

        # 步骤5: 计算t = KDF(x2 || y2, msg_len)
        x2_bytes = int_to_bytes(x2, (x2.bit_length() + 7) // 8)
        y2_bytes = int_to_bytes(y2, (y2.bit_length() + 7) // 8)
        t = kdf(x2_bytes + y2_bytes, msg_len)

        # 检查t是否为全0
        if all(b == 0 for b in t):
            continue  # 如果t全0，重新选择k

        # 步骤6: 计算C2 = M ⊕ t
        C2 = bytes(m ^ t[i] for i, m in enumerate(msg_bytes))

        # 步骤7: 计算C3 = Hash(x2 || M || y2)
        C3 = hash_msg(x2_bytes + msg_bytes + y2_bytes)

        # 步骤8: 输出密文C = C1 || C3 || C2
        C1_bytes = int_to_bytes(C1.x) + int_to_bytes(C1.y)
        end_time = time.time()
        return C1_bytes + C3 + C2, end_time - start_time


def sm2_decrypt(private_key: int, ciphertext: bytes) -> Tuple[str, float]:
    """SM2解密算法，返回解密结果和耗时(秒)"""
    start_time = time.time()

    # 解析密文
    point_len = (SM2_p.bit_length() + 7) // 8
    C1_bytes = ciphertext[:2 * point_len]
    C3_len = 32  # SM3输出长度(这里用SHA-256代替)
    C3 = ciphertext[2 * point_len:2 * point_len + C3_len]
    C2 = ciphertext[2 * point_len + C3_len:]
    msg_len = len(C2)

    # 步骤1: 从C1中取出椭圆曲线点
    x1 = bytes_to_int(C1_bytes[:point_len])
    y1 = bytes_to_int(C1_bytes[point_len:])
    C1 = ECPoint(x1, y1)

    # 步骤2: 验证C1是否满足椭圆曲线方程
    # (这里省略验证步骤)

    # 步骤3: 计算椭圆曲线点S = [h]C1
    S = scalar_mult(SM2_h, C1)
    if S.is_infinity():
        raise ValueError("S is infinity point")

    # 步骤4: 计算[dB]C1 = (x2, y2)
    dB_C1 = scalar_mult(private_key, C1).to_affine()
    x2 = dB_C1.x
    y2 = dB_C1.y

    # 步骤5: 计算t = KDF(x2 || y2, msg_len)
    x2_bytes = int_to_bytes(x2, (x2.bit_length() + 7) // 8)
    y2_bytes = int_to_bytes(y2, (y2.bit_length() + 7) // 8)
    t = kdf(x2_bytes + y2_bytes, msg_len)

    # 检查t是否为全0
    if all(b == 0 for b in t):
        raise ValueError("t is all zero")

    # 步骤6: 计算M' = C2 ⊕ t
    msg_bytes = bytes(c ^ t[i] for i, c in enumerate(C2))

    # 步骤7: 计算u = Hash(x2 || M' || y2)
    u = hash_msg(x2_bytes + msg_bytes + y2_bytes)

    # 步骤8: 验证u == C3
    if u != C3:
        raise ValueError("Hash verification failed")

    # 步骤9: 返回明文M'
    end_time = time.time()
    return msg_bytes.decode('utf-8'), end_time - start_time


# 示例使用
if __name__ == "__main__":
    # 生成密钥对
    private_key = random.randint(1, SM2_n - 1)
    public_key = scalar_mult(private_key, G)

    # 要加密的消息
    message = "zdh202200460110"
    print(f"原始消息: {message}")

    # 加密并计算耗时
    ciphertext, encrypt_time = sm2_encrypt(public_key, message)
    print(f"加密结果 (十六进制): {ciphertext.hex()}")
    print(f"加密耗时: {encrypt_time * 1000:.4f} 毫秒")

    # 解密并计算耗时
    decrypted, decrypt_time = sm2_decrypt(private_key, ciphertext)
    print(f"解密结果: {decrypted}")
    print(f"解密耗时: {decrypt_time * 1000:.4f} 毫秒")

    # 验证
    assert decrypted == message
    print("加密解密验证成功!")

    # 性能测试 - 多次运行取平均
    test_runs = 10
    total_encrypt_time = 0
    total_decrypt_time = 0

    for _ in range(test_runs):
        # 加密
        _, encrypt_time = sm2_encrypt(public_key, message)
        total_encrypt_time += encrypt_time

        # 解密
        _, decrypt_time = sm2_decrypt(private_key, ciphertext)
        total_decrypt_time += decrypt_time

    avg_encrypt_time = total_encrypt_time / test_runs
    avg_decrypt_time = total_decrypt_time / test_runs

    print(f"\n性能测试 (平均{test_runs}次运行):")
    print(f"平均加密耗时: {avg_encrypt_time * 1000:.4f} 毫秒")
    print(f"平均解密耗时: {avg_decrypt_time * 1000:.4f} 毫秒")