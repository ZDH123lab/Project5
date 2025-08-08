import hashlib
import random
import os
from dataclasses import dataclass

# 椭圆曲线参数保持不变
SM2_p = 0x8542D69E4C044F18E8B92435BF6FF7DE457283915C45517D722EDB8B08F1DFC3
SM2_a = 0x787968B4FA32C3FD2417842E73BBFEFF2F3C848B6831D7E0EC65228B3937E498
SM2_b = 0x63E4C6D3B23B0C849CF84241484BFE48F61D59A5B16BA06E6E12D1DA27C5249A
SM2_n = 0x8542D69E4C044F18E8B92435BF6FF7DD297720630485628D5AE74EE7C32E79B7
SM2_Gx = 0x421DEBD61B62EAB6746434EBC3CC315E32220B3BADD50BDC4C4E6C147FEDD43D
SM2_Gy = 0x0680512BCBB42C07D47349D2153B70C4E5D7FDFCBFA36EA1A85841B9E46E09A2


@dataclass
class Point:
    """椭圆曲线点类"""
    x: int = -1
    y: int = -1
    is_infinite: bool = False

    def __post_init__(self):
        """初始化点状态"""
        self.is_infinite = (self.x == -1 and self.y == -1)


def mod_inv(a, m):
    """计算模逆"""
    try:
        g, x, y = extended_gcd(a, m)
        if g != 1:
            return None
        return x % m
    except Exception as e:
        print(f"模逆计算失败: a={a}, m={m}, 错误={str(e)}")
        return None


def extended_gcd(a, b):
    """扩展欧几里得算法"""
    if a == 0:
        return (b, 0, 1)
    g, y, x = extended_gcd(b % a, a)
    return (g, x - (b // a) * y, y)


def point_add(P, Q):
    """椭圆曲线点加法"""
    if P.is_infinite:
        return Q
    if Q.is_infinite:
        return P
    if P.x == Q.x and P.y != Q.y:
        return Point()

    # 计算斜率
    if P.x == Q.x:
        numerator = (3 * P.x ** 2 + SM2_a) % SM2_p
        denominator = (2 * P.y) % SM2_p
    else:
        numerator = (Q.y - P.y) % SM2_p
        denominator = (Q.x - P.x) % SM2_p

    if denominator == 0:
        raise ValueError("斜率计算失败（分母为零）")

    inv_denominator = mod_inv(denominator, SM2_p)
    if inv_denominator is None:
        raise ValueError("斜率计算失败（分母不可逆）")

    lam = (numerator * inv_denominator) % SM2_p
    x = (lam ** 2 - P.x - Q.x) % SM2_p
    y = (lam * (P.x - x) - P.y) % SM2_p
    return Point(x, y)


def scalar_mult(k, P):
    """椭圆曲线标量乘法"""
    result = Point()
    addend = P

    while k:
        if k & 1:
            result = point_add(result, addend)
        addend = point_add(addend, addend)
        k >>= 1

    return result


def generate_keypair():
    """生成密钥对"""
    for _ in range(100):  # 最大尝试次数
        d = random.randint(1, SM2_n - 1)
        Q = scalar_mult(d, Point(SM2_Gx, SM2_Gy))
        if not Q.is_infinite:
            return d, Q
    raise RuntimeError("密钥对生成失败")


def hash_to_int(message):
    """哈希消息并转换为整数"""
    return int.from_bytes(hashlib.sha256(message).digest(), 'big') % SM2_n


def sign_ecdsa(private_key, message):
    """ECDSA签名"""
    z = hash_to_int(message)

    while True:
        k = random.randint(1, SM2_n - 1)
        R = scalar_mult(k, Point(SM2_Gx, SM2_Gy))

        if R.is_infinite:
            continue

        r = R.x % SM2_n
        k_inv = mod_inv(k, SM2_n)

        if k_inv is None:
            continue

        s = (z + r * private_key) * k_inv % SM2_n

        if 0 < r < SM2_n and 0 < s < SM2_n:
            return r, s


def verify_ecdsa(public_key, message, signature):
    """ECDSA签名验证"""
    r, s = signature

    if not (0 < r < SM2_n and 0 < s < SM2_n):
        return False

    z = hash_to_int(message)
    s_inv = mod_inv(s, SM2_n)

    if s_inv is None:
        return False

    u1 = z * s_inv % SM2_n
    u2 = r * s_inv % SM2_n

    P1 = scalar_mult(u1, Point(SM2_Gx, SM2_Gy))
    P2 = scalar_mult(u2, public_key)
    X = point_add(P1, P2)

    return not X.is_infinite and X.x % SM2_n == r


def main():
    """主函数"""
    try:
        private_key, public_key = generate_keypair()
        print("密钥对生成成功:")
        print(f"私钥: {private_key}")
        print(f"公钥: ({hex(public_key.x)[2:]}, {hex(public_key.y)[2:]})")

        message = "伪造中本聪的数字签名".encode('utf-8')
        print("\n===== 签名流程开始 =====")
        print(f"消息: {message.decode('utf-8')}")

        r, s = sign_ecdsa(private_key, message)
        print(f"\n签名生成成功!")
        print(f"r: {r}")
        print(f"s: {s}")

        print("\n验证签名...")
        is_valid = verify_ecdsa(public_key, message, (r, s))
        print(f"验证结果: {'✓ 通过' if is_valid else '✗ 失败'}")
    except Exception as e:
        print(f"错误: {str(e)}")


if __name__ == "__main__":
    main()