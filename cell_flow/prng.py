"""确定性随机数：自带的 32 位 MT19937，保证与运行环境无关。

标准库 ``random.Random`` 虽然也确定，但 Cell Flow 要求“相同输入、参数、版本
跨运行逐字节一致”，这里固定一个简单、明确的伪随机实现，避免对外部状态
（如 ``PYTHONHASHSEED``、库版本）的任何依赖。
"""


def _uint32(x: int) -> int:
    return x & 0xFFFFFFFF


class MT19937:
    """经典 32 位 Mersenne Twister（参数与参考实现一致）。"""

    _W = 32
    _N = 624
    _M = 397
    _A = 0x9908B0DF
    _F = 1812433253
    _U = 11
    _D = 0xFFFFFFFF
    _S = 7
    _B = 0x9D2C5680
    _T = 15
    _C = 0xEFC60000
    _L = 18

    def __init__(self, seed: int):
        if not isinstance(seed, int):
            raise TypeError("seed 必须是整数")
        state = [0] * self._N
        state[0] = _uint32(seed)
        for i in range(1, self._N):
            prev = state[i - 1]
            state[i] = _uint32(
                self._F * (prev ^ (prev >> (self._W - 2))) + i
            )
        self._state = state
        self._index = self._N

    def _twist(self) -> None:
        n, m, a = self._N, self._M, self._A
        lower_mask = _uint32(0x7FFFFFFF)
        upper_mask = _uint32(0x80000000)
        state = self._state
        for i in range(n):
            x = _uint32((state[i] & upper_mask) | (state[(i + 1) % n] & lower_mask))
            x_a = x >> 1
            if x & 1:
                x_a ^= a
            state[i] = state[(i + m) % n] ^ x_a
        self._index = 0

    def next_u32(self) -> int:
        if self._index >= self._N:
            self._twist()
        y = self._state[self._index]
        y ^= (y >> self._U) & self._D
        y ^= (y << self._S) & self._B
        y ^= (y << self._T) & self._C
        y ^= y >> self._L
        self._index += 1
        return _uint32(y)

    def random(self) -> float:
        """[0, 1) 上的均匀浮点（53 位，两个 32 位随机数拼成）。"""
        high = self.next_u32() >> 5  # 27 bits
        low = self.next_u32() >> 6  # 26 bits
        return (high * (1 << 26) + low) / float(1 << 53)

    def randbelow(self, n: int) -> int:
        """[0, n) 内的均匀整数，n >= 1。"""
        if n <= 0:
            raise ValueError("randbelow 需要正整数上界")
        if n == 1:
            return 0
        # 拒绝采样，严格均匀
        limit = 1 << 32
        bound = limit - (limit % n)
        while True:
            r = self.next_u32()
            if r < bound:
                return r % n

    def gaussian(self) -> float:
        """标准正态样本（Box-Muller，成对缓存另一个）。"""
        if self._spare is not None:
            value = self._spare
            self._spare = None
            return value
        # u1 取 (0,1]，避免 log(0)
        u1 = 1.0 - self.random()
        u2 = self.random()
        import math

        radius = math.sqrt(-2.0 * math.log(u1))
        angle = 2.0 * math.pi * u2
        z0 = radius * math.cos(angle)
        self._spare = radius * math.sin(angle)
        return z0

    _spare: float | None = None
