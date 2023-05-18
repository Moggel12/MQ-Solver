#ifndef MQ_UNI_H
#define MQ_UNI_H

//////////////// MACROS FOR OPERATIONS ON INTEGERS ////////////////
#define GF2_ADD(a, b) (a ^ b)
#define GF2_MUL(a, b) (a & b)
#define GRAY(i) (i ^ (i >> 1))

#define INT_LSHIFT(i, w) i << w
#define INT_RSHIFT(i, w) i >> w
#define INT_IDX(p, i) ((p >> i) & 1)
#define INT_SETBIT(p, i, b) (p ^ (b << i))
#define INT_FF (-1)
#define INT_0 (0)
#define INT_1 1
#define INT_LSB(i) (i & -i)
#define INT_IS_ZERO(p) (p == 0)
#define INT_MASK(b) ((1ull << b) - 1)
#define INT_EQ(a, b) a == b

#endif  // !MQ_UNI_H
