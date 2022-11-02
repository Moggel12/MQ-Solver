def rotate_right(x, n):
    return x[-n:] + x[:-n]

def pi(n):
    if n <= 1:
        return (0,)
    x = pi(n - 1) + (n - 1,)
    return rotate_right(tuple(x[k] for k in x), 1)

def p(n, j, reverse=False):
    if n == 1 and j == 0:
        if not reverse:
            yield (0,)
            yield (1,)
        else:
            yield (1,)
            yield (0,)
    elif j >= 0 and j < n:
        perm = pi(n - 1)
        if not reverse:
            for x in p(n - 1, j - 1):
                yield (1,) + tuple(x[k] for k in perm)
            for x in p(n - 1, j):
                yield (0,) + x
        else:
            for x in p(n - 1, j, reverse=True):
                yield (0,) + x
            for x in p(n - 1, j - 1, reverse=True):
                yield (1,) + tuple(x[k] for k in perm)

def monotonic(n):
    for i in range(n):
        for x in (p(n, i) if i % 2 == 0 else p(n, i, reverse=True)):
            yield tuple(reversed(x))

def monotonic_bounded(n, hw):
    prev = 0
    cur = 0 
    for j, i in enumerate(monotonic(n)):
        cur = sum(i)
        if (prev > hw) and (cur > hw):
            prev = cur
            continue
        if (cur <= hw):
            yield j, i
        prev = cur

def main():
    prev = 0
    cur  = 0

    show = False

    weight = 2

    for j, i in enumerate(monotonic(4)):
        cur = sum(i)

        if (prev != weight) and (cur != weight):
            show = False

        if (cur == weight):
            show = True

        #if sum(i) <= 2:
        if show:
            print(j, i)

        prev = cur

if __name__ == '__main__':
    main()