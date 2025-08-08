import hashlib as hs
from random import SystemRandom as SR
from dataclasses import make_dataclass

# SM2椭圆曲线参数
SM2_p = 0x8542D69E4C044F18E8B92435BF6FF7DE457283915C45517D722EDB8B08F1DFC3
SM2_a = 0x787968B4FA32C3FD2417842E73BBFEFF2F3C848B6831D7E0EC65228B3937E498
SM2_b = 0x63E4C6D3B23B0C849CF84241484BFE48F61D59A5B16BA06E6E12D1DA27C5249A
SM2_n = 0x8542D69E4C044F18E8B92435BF6FF7DD297720630485628D5AE74EE7C32E79B7
SM2_Gx = 0x421DEBD61B62EAB6746434EBC3CC315E32220B3BADD50BDC4C4E6C147FEDD43D
SM2_Gy = 0x0680512BCBB42C07D47349D2153B70C4E5D7FDFCBFA36EA1A85841B9E46E09A2

ECPoint = make_dataclass('ECPoint', [('x_coord', int, -1), ('y_coord', int, -1), ('inf_flag', bool, False)])

def _initialize_point(self):
    if self.x_coord == -1 and self.y_coord == -1:
        self.inf_flag = True
    else:
        self.inf_flag = False

ECPoint.__post_init__ = _initialize_point

def modular_inverse(val, mod):
    prev, curr = 1, 0
    a, b = val % mod, mod
    while a > 1:
        quotient = b // a
        new_prev = curr - prev * quotient
        new_a = b - a * quotient
        prev, curr = new_prev, prev
        a, b = new_a, a
    return prev % mod

def elliptic_add(P, Q):
    if P.inf_flag:
        return Q
    if Q.inf_flag:
        return P
    if P.x_coord == Q.x_coord:
        if P.y_coord != Q.y_coord:
            return ECPoint(-1, -1, True)
        slope = ((3 * P.x_coord ** 2 + SM2_a) * modular_inverse(2 * P.y_coord, SM2_p)) % SM2_p
    else:
        slope = ((Q.y_coord - P.y_coord) * modular_inverse(Q.x_coord - P.x_coord, SM2_p)) % SM2_p
    new_x = (slope ** 2 - P.x_coord - Q.x_coord) % SM2_p
    new_y = (slope * (P.x_coord - new_x) - P.y_coord) % SM2_p
    return ECPoint(new_x, new_y)

def point_multiplication(scalar, point):
    if point.inf_flag or scalar == 0:
        return ECPoint(-1, -1, True)
    res = ECPoint(-1, -1, True)
    current = point
    while scalar:
        if scalar & 1:
            res = elliptic_add(res, current)
        current = elliptic_add(current, current)
        scalar >>= 1
    return res

INF_POINT = ECPoint()

def create_key_pair():
    rng = SR()
    while True:
        secret = rng.randint(1, SM2_n - 1)
        pub_point = point_multiplication(secret, ECPoint(SM2_Gx, SM2_Gy))
        if not pub_point.inf_flag:
            break
    return secret, pub_point

def create_signature(priv_key, msg):
    rng = SR()
    k_val = rng.randint(1, SM2_n - 1)
    R_point = point_multiplication(k_val, ECPoint(SM2_Gx, SM2_Gy))
    if R_point.inf_flag:
        return None
    hash_val = int(hs.sha256(msg.encode()).hexdigest()[:64], 16) % SM2_n
    sig = (modular_inverse((k_val % SM2_n), SM2_n) * (hash_val + priv_key % SM2_n)) % SM2_n
    return (R_point, sig)

def check_signature(pub_key, msg, sig):
    R, s_val = sig
    if R.inf_flag or s_val < 1 or s_val > SM2_n - 1:
        return False
    msg_hash = int(hs.sha256(msg.encode()).hexdigest()[:64], 16) % SM2_n
    inv_s = modular_inverse(s_val % SM2_n, SM2_n)
    u1 = (msg_hash * inv_s) % SM2_n
    u2 = (s_val * inv_s) % SM2_n
    X = (u1 * SM2_Gx + u2 * pub_key.x_coord) % SM2_p
    Y = (u1 * SM2_Gy + u2 * pub_key.y_coord) % SM2_p
    return R == ECPoint(X, Y)

# 用户密钥生成和签名过程
priv_A, pub_A = create_key_pair()
priv_B, pub_B = create_key_pair()

msg_A = "Message from Alice"
msg_B = "Message from Bob"

sig_A = create_signature(priv_A, msg_A)
sig_B = create_signature(priv_B, msg_B)

print("Alice签名数据:", sig_A)
print("Bob签名数据:", sig_B)

# 私钥推导过程
hash_A = int(hs.sha256(msg_A.encode()).hexdigest()[:64], 16) % SM2_n
R_A, s_A = sig_A
k_A = modular_inverse((s_A % SM2_n), SM2_n) * (hash_A - R_A.x_coord) % SM2_n
derived_priv_A = (hash_A - k_A * R_A.x_coord) % SM2_n
print("推导Alice私钥:", derived_priv_A)

hash_B = int(hs.sha256(msg_B.encode()).hexdigest()[:64], 16) % SM2_n
R_B, s_B = sig_B
k_B = modular_inverse((s_B % SM2_n), SM2_n) * (hash_B - R_B.x_coord) % SM2_n
derived_priv_B = (hash_B - k_B * R_B.x_coord) % SM2_n
print("推导Bob私钥:", derived_priv_B)