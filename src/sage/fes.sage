import math
from random import randint

def run_fes(f_sys, n):
  for f in f_sys: # Should be changed to the "clamped" version
    s = init(f, n)
    while s["i"] < 2^n:
      next(s)
      if s["y"] == 0:
        return s["y"]

def next(s):
  s["i"] += 1
  k1 = BIT1(s["i"])
  k2 = BIT2(s["i"])
  if math.log2(s["i"]).is_integer():
    s["deriv1"][k1] = (s["deriv1"][k1] + s["deriv2"][k1][k2])
  s["y"] = (s["y"] + s["deriv1"][k1])

def init(f, n):
  s = dict()
  s["i"] = 0
  s["x"] = 0
  s["y"] = f.constant_coefficient()
  s["deriv2"] = dict()
  for k in range(1, n):
    s["deriv2"][k] = dict()
    for j in range(k):
      s["deriv2"][k][j] = f.derivative(f.variables()[k], f.variables()[j])
  s["deriv1"] = dict()
  s["deriv1"][0] = f.coefficient({v: (1 if i == 0 else 0) for i, v in enumerate(f.variables())})
  for k in range(1, n):
    s["deriv1"][k] = s["deriv2"][k][k-1] + f.coefficient({v: (1 if i == k else 0) for i, v in enumerate(f.variables())})
  return s


def BIT1(i):
  for idx, c in enumerate(reversed(bin(i ^ (i >> 1))[2:])):
    if c == 1: return idx

def BIT2(i):
  return len(bin(i ^ (i >> 1))[2:]) - 1

def main():
  for _ in range(10):
    m = randint(5, 10)
    n = randint(2, 5)
    R = GF(2)[", ".join(["x" + str(i) for i in range(n)])]
    f_sys = []
    for i in range(m):
      f_sys.append(R(GF(2)[R.gens()].random_element(degree=2)))
    run_fes(f_sys, n)

if __name__ == "__main__":
  main() 
