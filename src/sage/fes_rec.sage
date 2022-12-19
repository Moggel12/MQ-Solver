from itertools import combinations

# Get index position of first bit set (if any).
def bit1(x):
  x = x&-x

  x = int(x).bit_length()

  return None if x == 0 else x-1

# Get indices of set bits as array.
def bits(x):
  if x == 0:
   return None

  ret = []
  while x > 0: 
    ret.append(bit1(x))

    x = x ^^ (x&-x)

  return ret


def fes_eval(f, part = None):
  n = len(f.parent().gens())

  if part == None:
    part = n

  res = [None] * 2^n

  X = R.gens()

  d = {}

  d[0] = f.constant_coefficient()

  for deg in reversed(range(1, f.degree()+1)):
    for idx in combinations(range(n), deg):
      x = [1 if i in idx else 0 for i in range(n)]

      ## evaluate partial derivative in corresponding Gray-code order
      #val = derivative(f, *[X[i] for i in idx])(*[a ^^ b for a,b in zip(x, x[1:] + [0])])

      val = f.monomial_coefficient(prod([X[i] for i in idx]))

      for ddeg in range(1, deg+1): #f.degree()-1):
        for didx in combinations([i-1 for i in idx if i-1 >= 0 and i-1 not in idx], ddeg):
          if len(didx) + len(idx) <= f.degree():
            val += f.monomial_coefficient(prod([X[i] for i in didx] + [X[i] for i in idx]))
  
      d[sum([2^i for i in idx])] = val

  res[0] = d[0]

  for si in range(1, 2^n):
    k = bits(si)[:f.degree()]
  
    for j in reversed(range(0, len(k))):
      d[sum([2^i for i in k[:j]])] = d[sum([2^i for i in k[:j]])] + d[sum([2^i for i in k[:j+1]])]
  
    if len(bits(si)) <= part:
      res[si ^^ (si >> 1)] = d[0]

  return res

def fes_recover(evall, n, degree, f = None):  # parameter f for debugging polynomial interpolation
  #n = len(f.parent().gens())
  #degree = f.degree()

  res = [None] * 2^n

  d = {}

  res[0] = evall[0]
  d[0] = evall[0]

  for si in range(1, 2^n):
    if len(bits(si)) > degree:
      # We have the required derivatives; compute the missing evluation value.

      k = bits(si)[:degree]

      for j in reversed(range(0, len(k))):
        d[sum([2^i for i in k[:j]])] = d[sum([2^i for i in k[:j]])] ^^ d[sum([2^i for i in k[:j+1]])]
    
      res[si ^^ (si >> 1)] = d[0]

    else:
      # We need to interpolate s.d1 and s.d2.

      k = bits(si)[:degree]
  
      i = si ^^ (si >> 1)
    
      prev = d[0]
      d[0] = evall[i]
    
      for j in range(1, len(k)+1):
        if j < len(k):
          tmp = d[sum([2^i for i in k[:j]])]
  
        d[sum([2^i for i in k[:j]])] = d[sum([2^i for i in k[:j-1]])] ^^ prev
  
        if j < len(k):
          prev = tmp
    
      res[si ^^ (si >> 1)] = d[0]

##  # For quadratic polynomails:
##  # Recover linear coefficients by adding up to two second derivatives:
##  for i in range(0, n-1):
##    #s.d[1][i] += s.d[2][n-1][i]
##    d[2^i] += d[2^i + 2^(n-1)]
##  
##  for i in range(1, n):
##    #s.d[1][i] += s.d[2][i][i-1]
##    d[2^i] += d[2^i + 2^(i-1)]
##
##  g = evall[0] + \
##      sum([X[i] * d[2^i] for i in range(n)]) + \
##      sum([X[i]*X[j] * d[2^i + 2^j] for i in range(1, n) for j in range(0,i)])
##
##  print()
##  print(g)
##  print(f == g)
##  print()

  return res

def test():
    n = 8
    degree = 3

    R = GF(2)[", ".join(["x" + str(i) for i in range(n)])]

    print(R)
    print()


    X = R.gens()

    f = GF(2).random_element()

    for d in range(1, degree+1):
      for mon in combinations(range(n), d):
        f += prod([X[i] for i in mon]) * GF(2).random_element()
      

    print("Input polynomial:", f)
    print()

    print("Performing full evaluation...")
    eval_full = fes_eval(f)

    print(eval_full)

    check = [f(*reversed([int(i) for i in f"{x:0{n}b}"])) for x in range(2^n)]

    print(check)

    print(eval_full == check)


    print("\nPerforming partial evaluation and recovery...")
    eval_part = fes_eval(f, f.degree())

    print("[" + ", ".join(['-' if v == None else str(v) for v in eval_part]) + "]")

    recovered = fes_recover(eval_part, len(f.parent().gens()), f.degree(), f)

    print(recovered)

    print(recovered == eval_full)

if __name__ == "__main__":
    test()
