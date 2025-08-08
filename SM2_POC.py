import hashlib
import random
from dataclasses import dataclass

# SM2曲线参数
p = int('0x8542D69E4C044F18E8B92435BF6FF7DE457283915C45517D722EDB8B08F1DFC3', 16)
a = int('0x787968B4FA32C3FD2417842E73BBFEFF2F3C848B6831D7E0EC65228B3937E498', 16)
b = int('0x63E4C6D3B23B0C849CF84241484BFE48F61D59A5B16BA06E6E12D1DA27C5249A', 16)
n = int('0x8542D69E4C044F18E8B92435BF6FF7DD297720630485628D5AE74EE7C32E79B7', 16)
Gx = int('0x421DEBD61B62EAB6746434EBC3CC315E32220B3BADD50BDC4C4E6C147FEDD43D', 16)
Gy = int('0x0680512BCBB42C07D47349D2153B70C4E5D7FDFCBFA36EA1A85841B9E46E09A2', 16)


@dataclass
class Point:
    x: int = -1
    y: int = -1
    is_infinite: bool = False

    def __post_init__(self):
        if self.x == -1 and self.y == -1:
            self.is_infinite = True
        else:
            self.is_infinite = False

    def __eq__(self, other):
        if self.is_infinite and other.is_infinite:
            return True
        return self.x == other.x and self.y == other.y


def mod_inv(a, modulus):
    """扩展欧几里得算法求模逆"""
    if modulus == 1:
        return 0
    lm, hm = 1, 0
    low, high = a % modulus, modulus
    while low > 1:
        ratio = high // low
        nm, new = hm - lm * ratio, high - low * ratio
        lm, hm, low, high = nm, lm, new, low
    return lm % modulus


def point_add(P, Q):
    """椭圆曲线点加法"""
    if P.is_infinite:
        return Q
    if Q.is_infinite:
        return P
    if P.x == Q.x:
        if P.y != Q.y:
            return Point(is_infinite=True)
        lam = (3 * P.x ** 2 + a) * mod_inv(2 * P.y, p) % p
    else:
        lam = (Q.y - P.y) * mod_inv(Q.x - P.x, p) % p
    x3 = (lam ** 2 - P.x - Q.x) % p
    y3 = (lam * (P.x - x3) - P.y) % p
    return Point(x3, y3)


def scalar_mult(k, P):
    """椭圆曲线标量乘法（double-and-add算法）"""
    if k == 0 or P.is_infinite:
        return Point(is_infinite=True)
    result = Point(is_infinite=True)
    current = P
    while k:
        if k & 1:
            result = point_add(result, current)
        current = point_add(current, current)
        k >>= 1
    return result


def generate_keypair():
    """生成SM2密钥对"""
    d = random.randint(1, n - 1)
    Q = scalar_mult(d, Point(Gx, Gy))
    return d, Q


def sm2_sign(private_key, msg):
    """SM2签名算法"""
    k = random.randint(1, n - 1)
    R = scalar_mult(k, Point(Gx, Gy))
    if R.is_infinite:
        return None
    e = int(hashlib.sha256(msg.encode()).hexdigest(), 16) % n
    r = (e + R.x) % n
    s = mod_inv(1 + private_key, n) * (k - r * private_key) % n
    return (r, s, k, R)  # 返回k和R用于POC验证


def ecdsa_sign(private_key, msg):
    """ECDSA签名算法"""
    k = random.randint(1, n - 1)
    R = scalar_mult(k, Point(Gx, Gy))
    if R.is_infinite:
        return None
    e = int(hashlib.sha256(msg.encode()).hexdigest(), 16) % n
    r = R.x % n
    s = mod_inv(k, n) * (e + r * private_key) % n
    return (r, s, k, R)


# 三种签名误用情况的POC验证

def case1_same_user_same_k():
    """情况1: 同一用户对两个消息使用相同k"""
    dA, _ = generate_keypair()
    msg1, msg2 = "Hello", "World"

    # 使用相同k签名两个消息
    _, s1, k, R = sm2_sign(dA, msg1)
    _, s2, _, _ = sm2_sign(dA, msg2)

    # 推导私钥（根据文档公式）
    numerator = (s2 - s1) % n
    denominator = (s1 - s2 + R.x - R.x) % n  # 实际应为(s1 - s2) + (r1 - r2)但r相同
    dA_derived = numerator * mod_inv(denominator, n) % n

    return dA == dA_derived


def case2_different_users_same_k():
    """情况2: 两个用户使用相同k"""
    dA, _ = generate_keypair()
    dB, _ = generate_keypair()
    msg1, msg2 = "Alice", "Bob"
    k = random.randint(1, n - 1)  # 共享相同k

    # Alice签名
    R = scalar_mult(k, Point(Gx, Gy))
    e1 = int(hashlib.sha256(msg1.encode()).hexdigest(), 16) % n
    r1 = (e1 + R.x) % n
    s1 = mod_inv(1 + dA, n) * (k - r1 * dA) % n

    # Bob签名（使用相同k）
    e2 = int(hashlib.sha256(msg2.encode()).hexdigest(), 16) % n
    r2 = (e2 + R.x) % n
    s2 = mod_inv(1 + dB, n) * (k - r2 * dB) % n

    # Alice推导Bob私钥
    dB_derived = (k - s2) * mod_inv(s2 + r2, n) % n
    # Bob推导Alice私钥
    dA_derived = (k - s1) * mod_inv(s1 + r1, n) % n

    return (dA == dA_derived, dB == dB_derived)


def case3_same_d_and_k_different_algorithms():
    """情况3: 相同d和k用于ECDSA和SM2"""
    d, _ = generate_keypair()
    msg = "CriticalMessage"
    k = random.randint(1, n - 1)  # 共享k

    # ECDSA签名
    R = scalar_mult(k, Point(Gx, Gy))
    e1 = int(hashlib.sha256(msg.encode()).hexdigest(), 16) % n
    r1 = R.x % n
    s1 = mod_inv(k, n) * (e1 + r1 * d) % n

    # SM2签名（相同d和k）
    e2 = int(hashlib.sha256(msg.encode()).hexdigest(), 16) % n
    r2 = (e2 + R.x) % n
    s2 = mod_inv(1 + d, n) * (k - r2 * d) % n

    # 推导私钥
    term1 = e1 * mod_inv(s1, n) % n
    term2 = s2 * (1 + d) % n
    term3 = r2 * d % n
    d_derived = (term1 - s2) * mod_inv(s2 + r2 - r1 * mod_inv(s1, n), n) % n

    return d == d_derived

# 执行POC验证

if __name__ == "__main__":
    print("情况1验证结果:", "成功" if case1_same_user_same_k() else "失败")

    resultA, resultB = case2_different_users_same_k()
    print(f"情况2验证: Alice私钥推导{'成功' if resultA else '失败'}, Bob私钥推导{'成功' if resultB else '失败'}")

    print("情况3验证结果:", "成功" if case3_same_d_and_k_different_algorithms() else "失败")