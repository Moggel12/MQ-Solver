
def convert(v, n):
    v = bin(v)[2:]
    return list(map(lambda i : (int(v[-i]) if i <= len(v) else 0), range(1, n + 1)))

def index_of(y_list):
    return sum(b << i for i, b in enumerate(y_list))
