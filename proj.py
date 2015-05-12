#!/usr/bin/env python
import ctypes
import matplotlib.pyplot as plt
import numpy as np

def mkdirs(path):
    import os
    try:
        os.makedirs(path)
    except OSError:
        pass

def signature(*args, **kwargs):
    import base64, hashlib, pickle
    return base64.b32encode(hashlib.sha1(
        pickle.dumps((args, kwargs))).digest()).lower()

def cached(name_template, load, dump):
    def go(func):
        def go(*args, **kwargs):
            fn = name_template.format(signature(*args, **kwargs)
                                      .decode("latin-1"))
            try:
                return load(fn)
            except OSError:
                pass
            import os
            result = func(*args, **kwargs)
            mkdirs(os.path.dirname(fn))
            dump(fn, result)
            return result
        return go
    return go

import scipy.special as sps
def bessel_j(l, z):
    '''Spherical Bessel function of the first kind'''
    z = complex(z.real, z.imag)
    return sps.sph_jn(l, z)[0][-1]

def map_bessel_j(l, zs):
    import numpy as np
    return np.array([bessel_j(l, z) for z in zs])

def simple_figure():
    import matplotlib.pyplot as plt
    fig  = plt.figure()
    axes = fig.add_subplot(111)
    return fig, axes

# ----------------------------------------------------------------------------
# Integration
# ----------------------------------------------------------------------------

def integrate(f, a, b, **kwargs):
    '''Similar to spi.quad but works on complex numbers.
    Only returns the value of the integral, however.'''
    from scipy.integrate import quad
    re_integral, _ = quad(lambda z: f(z).real, a, b, **kwargs)
    im_integral, _ = quad(lambda z: f(z).imag, a, b, **kwargs)
    return complex(re_integral, im_integral)

def trapezoidal_rule(count):
    import numpy as np
    xs = np.linspace(-1, 1, count)
    ws = np.full_like(xs, 2 / (count - 1))
    ws[0]  /= 2
    ws[-1] /= 2
    return xs, ws

def make_node_generator(rule):
    def generate_nodes(x1, x2, count):
        import numpy as np
        xs, ws = np.polynomial.legendre.leggauss(count)
        c = (x2 - x1) / 2.
        return c * xs + (x2 + x1) / 2., c * ws
    return generate_nodes

def piecewise_curve(start, segments, node_generator):
    import numpy as np
    xs = []
    ws = []
    for point, count in segments:
        x, w = node_generator(start, point, count)
        start = point
        xs.append(x)
        ws.append(w)
    return np.concatenate(xs), np.concatenate(ws)

def triangle_curve(k_start, k_max, k_peak, k_inflection,
                   count, node_generator):
    segment_count = count // 3
    return piecewise_curve(k_start, (
        (k_peak, segment_count),
        (k_inflection, segment_count),
        (k_max, count - 2 * segment_count),
    ), node_generator)

gauss_legendre_nodes = make_node_generator(np.polynomial.legendre.leggauss)
trapezoidal_nodes    = make_node_generator(trapezoidal_rule)

# ----------------------------------------------------------------------------
# C utility
# ----------------------------------------------------------------------------

def run_interruptibly(func):
    '''Run a given function in a background thread, attempting to join every
    100ms.  Useful for C calls that may block for a long time, allowing the
    user to interrupt the call.  Returns the result of the call.'''
    import threading
    ret = []
    thread = threading.Thread(target=lambda: ret.append(func()))
    thread.daemon = True
    thread.start()
    while thread.is_alive():
        thread.join(.1)
    return ret[0]

def ctypes_vector(x):
    '''Marshal a vector into C.'''
    import ctypes
    assert x.dtype in (float, complex)
    assert len(x.shape) == 1
    return (
        x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_size_t(x.shape[0]),
    )

def ctypes_matrix(x):
    '''Marshal a matrix into C.'''
    import ctypes
    assert x.dtype in (float, complex)
    assert len(x.shape) == 2
    assert x.flags["C_CONTIGUOUS"]
    return (
        x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_size_t(x.shape[0]),
        ctypes.c_size_t(x.shape[1]),
    )

# ----------------------------------------------------------------------------

class WoodsSaxonParams(ctypes.Structure):
    _fields_ = (
        ("V0",  ctypes.c_double),
        ("Rws", ctypes.c_double),
        ("aws", ctypes.c_double),
    )

def get_woods_saxon(ws):
    from numpy import exp
    return lambda R: ws.V0 / (1 + exp((R - ws.Rws) / ws.aws))

def V_eff(V, l, mu2):
    return lambda R: l * (l + 1) / R ** 2 / mu2 + V(R)

libproj = ctypes.cdll.LoadLibrary("./libproj.so")
@cached(__file__ + ".cache/v-matrix/{0}.npy", np.load, np.save)
def calc_V_matrix(k, l, ws, R_max, abs_err, rel_err, limit):
    import numpy as np
    print("generating V matrix...")
    Vm = np.empty((len(k), len(k)), dtype=complex)
    run_interruptibly(lambda: libproj.generate_v_matrix(*(
        ctypes_matrix(Vm)[:1] +
        ctypes_vector(k) +
        (ctypes.c_uint(l),
         ws,
         ctypes.c_double(R_max),
         ctypes.c_double(abs_err),
         ctypes.c_double(rel_err),
         ctypes.c_uint(limit))
    )))
    print("done.")
    return Vm

@cached(__file__ + ".cache/bessel-matrix/{0}.npy", np.load, np.save)
def calc_bessel_matrix(k, l, R):
    import numpy as np
    print("generating Bessel function matrix...")
    len_R = len(R)
    len_k = len(k)
    Jm = np.empty((len_R, len_k), dtype=complex)
    for i in range(len_R):
        for j in range(len_k):
            Jm[i, j] = bessel_j(l, k[j] * R[i])
    print("done.")
    return Jm

def plot_bessel():
    import numpy as np
    from numpy import abs
    R = np.linspace(0, 20, 150)
    fig, axes = simple_figure()
    axes.plot(R, abs(R * map_bessel_j(l, R)) ** 2)
    axes.plot(R, abs(R * map_bessel_j(l, (.8 + 0.1j) * R)) ** 2)
    axes.set_ylim(0, 10)
    return fig

def plot_interactions(V, l_range, mu2):
    import numpy as np
    R = np.linspace(0.1, 10)
    fig, axes = simple_figure()
    for l in l_range:
        axes.plot(R, V_eff(V, l, mu2)(R), label=str(l))
    axes.set_ylim((-70, 70))
    axes.legend()
    return fig

def solve(mu2, ws, l, R_max, count, k_start, ka, kb, k_max,
          abs_err=1e-8, rel_err=1e-8, gk_limit=50,
          node_generator=gauss_legendre_nodes):
    import numpy as np
    from numpy import sqrt
    k, w = triangle_curve(k_start, k_max, ka, kb, count, node_generator)
    Vm = calc_V_matrix(k, l, ws, R_max, abs_err, rel_err, gk_limit)
    Hm = np.empty((count, count), dtype=complex)
    for i in range(count):
        for j in range(count):
            z = sqrt(w[i] * w[j]) * k[i] * k[j] * Vm[i, j]
            if i == j:
                z += k[i] ** 2 / mu2
            Hm[i, j] = z
    Es, phis = np.linalg.eig(Hm)
    return k, w, Es, phis

def find_closest(target, array):
    from numpy import abs
    return array[abs(array - target).argmin()]

# ----------------------------------------------------------------------------

def simple_run():
    import numpy as np
    from numpy import abs, pi, sqrt
    k, w, Es, phis = solve(
        mu2      = mu2,
        ws       = ws,
        l        = l,
        R_max    = 40., # fm
        count    = 80,
        k_start  =  0.01,
        ka       =  1. - .3j,
        kb       =  2.,
        k_max    = 10.,
    )

    print("energies:")
    print(Es)

    eig_k = sqrt(mu2 * Es + 0j)
    fig1, axes = simple_figure()
    axes.scatter(k.real, k.imag, 50, "blue", linewidth=0)
    axes.scatter(eig_k.real, eig_k.imag, 50, "red", linewidth=0)

    R = np.linspace(0.001, 200, 2000)
    Jm = calc_bessel_matrix(k, l, R)

    fig2, axes = simple_figure()
    for E, phi_ in zip(Es, phis.T):
        if not (abs(E - (1.346-0.1118j)) < 0.1 or E < 0):
           continue
        phi = sqrt(w) * k * phi_
        u = R * 1j ** l * sqrt(2 / pi) * np.dot(Jm, phi)

        axes.plot(R, V_eff(V, l, mu2)(R) * .01, label="V")
        axes.plot(R, abs(u)**2, label=repr(E))
        axes.set_ylim((-.6, .6))
        axes.legend()

    return fig1, fig2

def convergence_run1():
    import matplotlib.cm as cm
    import numpy as np
    from numpy import sqrt
    fig, axes = simple_figure()
    counts = tuple(range(30, 60, 1)) + tuple(range(60, 100, 2))
    colors = cm.rainbow(np.linspace(0, 1, len(counts)))
    for (count, c) in zip(counts, colors):
        _, _, Es, _ = solve(
            mu2      = mu2,
            ws       = ws,
            l        = l,
            R_max    = 40.,
            count    = count,
            k_start  =  0.01,
            ka       =  1. - .3j,
            kb       =  2.,
            k_max    = 10.,
        )
        ek = sqrt(mu2 * Es + 0j)
        axes.scatter(ek.real, ek.imag, 50, color=c,
                     linewidth=0, label=repr(count))
    axes.legend()
    return fig

def convergence_run2():
    import numpy as np
    from numpy import abs, sqrt
    fig, axes = simple_figure()
    xs = np.linspace(.02, .1, 9)
    k_last = None
    changes = []
    for x in xs:
        _, _, Es, _ = solve(
            mu2      = mu2,
            ws       = ws,
            l        = l,
            R_max    = 40.,
            count    = 120,
            k_start  =  0.01,
            ka       =  .25 - x * 1j,
            kb       =  2.,
            k_max    =  6.,
        )
        ek = sqrt(mu2 * Es + 0j)
        k = find_closest(0.254003 - 0.0105314j, ek)
        if k_last is not None:
            changes.append(abs(k - k_last))
        k_last = k
    axes.plot(xs[1:], changes)
    axes.set_yscale("log")
    return fig

A  = 10                        # (mass of core)
ws = WoodsSaxonParams(
    V0  = -61.1,               # MeV
    aws =   0.65,              # fm
    Rws =   1.2 * A ** (1/3.), # fm
)
l   = 2
mu2 = 0.0478450                # 1/(MeV*fm^2)
V   = get_woods_saxon(ws)

simple_figs = simple_run()
convergence_fig1 = convergence_run1()
convergence_fig2 = convergence_run2()
plt.show()
