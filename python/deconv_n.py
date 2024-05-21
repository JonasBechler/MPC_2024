import numpy as np
import control as ct

def deconv_n(A, B, n):
    len_A = len(A)
    len_B = len(B)
    max_len = max(len_A, len_B, n)

    if len_A < max_len:
        A = np.concatenate((A, np.zeros(max_len - len_A)))

    x, r = np.polydiv(A[:len_B], B)

    A = np.concatenate((r, A[len_B:]))

    if n == 1:
        C = x
        D = np.concatenate(([0], A[1:]))
    else:
        C_next, D_next = deconv_n(A[1:], B, n - 1)
        C = np.concatenate((x, C_next))
        D = np.concatenate(([0], D_next))

    return C, D