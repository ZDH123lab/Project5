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


class ECPoint:
    """椭圆曲线点类"""

    def __init__(self, x: int, y: int, z: int = 1):
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other) -> bool:
        if not isinstance(other, ECPoint):
            return False
        return self.x == other.x and self.y == other.y and self.z == other.z

    def is_infinity(self) -> bool:
        return self.z == 0

    def to_affine(self) -> 'ECPoint':
        if self.is_infinity():
            return ECPoint(0, 0, 0)
        z_inv = pow(self.z, SM2_p - 2, SM2_p)
        x = (self.x * z_inv * z_inv) % SM2_p
        y = (self.y * z_inv * z_inv * z_inv) % SM2_p
        return ECPoint(x, y, 1)


G = ECPoint(SM2_Gx, SM2_Gy)


def point_add(P: ECPoint, Q: ECPoint) -> ECPoint:
    if P.is_infinity():
        return Q
    if Q.is_infinity():
        return P
    if P == Q:
        return point_double(P)

    # 点加公式
    Z1Z1 = (P.z * P.z) % SM2_p
    Z2Z2 = (Q.z * Q.z) % SM2_p
    U1 = (P.x * Z2Z2) % SM2_p
    U2 = (Q.x * Z1Z1) % SM2_p
    S1 = (P.y * Q.z * Z2Z2) % SM2_p
    S2 = (Q.y * P.z * Z1Z1) % SM2_p

    H = (U2 - U1) % SM2_p
    I = (4 * H * H) % SM2_p
    J = (H * I) % SM2_p
    r = (2 * (S2 - S1)) % SM2_p
    V = (U1 * I) % SM2_p

    X3 = (r * r - J - 2 * V) % SM2_p
    Y3 = (r * (V - X3) - 2 * S1 * J) % SM2_p
    Z3 = ((P.z + Q.z) * (P.z + Q.z) - Z1Z1 - Z2Z2) % SM2_p
    Z3 = (Z3 * H) % SM2_p

    return ECPoint(X3, Y3, Z3)


def point_double(P: ECPoint) -> ECPoint:
    if P.is_infinity():
        return P

    # 点倍公式
    A = (P.x * P.x) % SM2_p
    B = (P.y * P.y) % SM2_p
    C = (B * B) % SM2_p
    D = (2 * ((P.x + B) * (P.x + B) - A - C)) % SM2_p
    E = (3 * A) % SM2_p
    F = (E * E) % SM2_p

    X3 = (F - 2 * D) % SM2_p
    Y3 = (E * (D - X3) - 8 * C) % SM2_p
    Z3 = (2 * P.y * P.z) % SM2_p

    return ECPoint(X3, Y3, Z3)


def point_mul(k: int, P: ECPoint) -> ECPoint:
    result = ECPoint(0, 0, 0)
    current = P

    while k > 0:
        if k & 1:
            result = point_add(result, current)
        current = point_double(current)
        k >>= 1

    return result


def bytes_to_int(b: bytes) -> int:
    return int.from_bytes(b, 'big')


def int_to_bytes(x: int, length: int = None) -> bytes:
    if length is None:
        length = (x.bit_length() + 7) // 8 or 1
    return x.to_bytes(length, 'big')


def hash_msg(msg: bytes) -> bytes:
    return hashlib.sha256(msg).digest()


def kdf(z: bytes, klen: int) -> bytes:
    hash_len = 32
    reps = (klen + hash_len - 1) // hash_len
    derived = b''
    for i in range(reps):
        derived += hash_msg(z + int_to_bytes(i + 1, 4))
    return derived[:klen]


def sm2_encrypt(public_key: ECPoint, msg: str) -> bytes:
    msg_bytes = msg.encode('utf-8')
    msg_len = len(msg_bytes)

    while True:
        k = random.randint(1, SM2_n - 1)
        C1 = point_mul(k, G).to_affine()
        S = point_mul(k, public_key).to_affine()

        if S.is_infinity():
            continue

        x2 = S.x
        y2 = S.y

        t = kdf(int_to_bytes(x2) + int_to_bytes(y2), msg_len)
        if all(b == 0 for b in t):
            continue

        C2 = bytes(m ^ t[i] for i, m in enumerate(msg_bytes))
        C3 = hash_msg(int_to_bytes(x2) + msg_bytes + int_to_bytes(y2))

        C1_bytes = int_to_bytes(C1.x) + int_to_bytes(C1.y)
        return C1_bytes + C2 + C3


def sm2_decrypt(private_key: int, ciphertext: bytes) -> str:
    point_len = (SM2_p.bit_length() + 7) // 8
    C1_x = bytes_to_int(ciphertext[:point_len])
    C1_y = bytes_to_int(ciphertext[point_len:2 * point_len])
    C1 = ECPoint(C1_x, C1_y)

    hash_len = 32
    C3 = ciphertext[-hash_len:]
    C2 = ciphertext[2 * point_len:-hash_len]
    msg_len = len(C2)

    S = point_mul(private_key, C1).to_affine()
    if S.is_infinity():
        raise ValueError("解密失败: S是无穷远点")

    x2 = S.x
    y2 = S.y

    t = kdf(int_to_bytes(x2) + int_to_bytes(y2), msg_len)
    if all(b == 0 for b in t):
        raise ValueError("解密失败: t是全0")

    msg_bytes = bytes(c ^ t[i] for i, c in enumerate(C2))
    u = hash_msg(int_to_bytes(x2) + msg_bytes + int_to_bytes(y2))

    if u != C3:
        raise ValueError("解密失败: 哈希验证不通过")

    return msg_bytes.decode('utf-8')


def generate_keypair() -> Tuple[int, ECPoint]:
    private_key = random.randint(1, SM2_n - 1)
    public_key = point_mul(private_key, G)
    return private_key, public_key


def main():
    # 生成密钥对
    private_key, public_key = generate_keypair()
    print(f"私钥: {hex(private_key)}")
    print(f"公钥: ({hex(public_key.x)}, {hex(public_key.y)})")

    # 要加密的消息
    message = "zdh202200460110"
    print(f"\n原始消息: {message}")

    # 加密
    start_enc = time.time()
    ciphertext = sm2_encrypt(public_key, message)
    enc_time = time.time() - start_enc
    print(f"\n加密耗时: {enc_time:.6f}秒")
    print(f"加密结果 (十六进制): {ciphertext.hex()}")
    print(f"密文长度: {len(ciphertext)} 字节")

    # 解密
    start_dec = time.time()
    try:
        decrypted_msg = sm2_decrypt(private_key, ciphertext)
        dec_time = time.time() - start_dec
        print(f"\n解密耗时: {dec_time:.6f}秒")
        print(f"解密结果: {decrypted_msg}")

        if decrypted_msg == message:
            print("解密成功: 解密结果与原始消息一致")
        else:
            print("解密失败: 解密结果与原始消息不一致")

        print(f"\n总耗时: {enc_time + dec_time:.6f}秒")
        print(f"加密/解密速度比: {enc_time / dec_time:.2f}")

    except ValueError as e:
        print(f"\n解密失败: {str(e)}")


if __name__ == "__main__":
    main()