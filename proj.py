#!/usr/bin/env python
import ctypes, itertools, multiprocessing
import matplotlib.pyplot as plt
import numpy as np

libproj = ctypes.cdll.LoadLibrary("./libproj.so")

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
            except IOError:
                pass
            import os
            result = func(*args, **kwargs)
            mkdirs(os.path.dirname(fn))
            dump(fn, result)
            return result
        return go
    return go

def bessel_j(l, z):
    '''Spherical Bessel function of the first kind'''
    import numpy
    j = numpy.asanyarray(z, dtype=complex).flatten()
    run_interruptibly(lambda: libproj.spherical_bessel_jl_many(*(
        (ctypes_vector(j)[0], ctypes.c_double(l)) +
        ctypes_vector(j)
    )))
    if isinstance(j, numpy.ndarray):
        j = j.reshape(z.shape)
    return j

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
    if not ret:
        raise Exception("background thread failed")
    return ret[0]

def ctypes_vector(x):
    '''Marshal a vector into C.'''
    import ctypes, numpy
    if x is None:
        return (ctypes.POINTER(ctypes.c_double)(), ctypes.c_size_t(0))
    assert x.dtype in (float, numpy.float64, complex, numpy.complex128)
    assert len(x.shape) == 1
    return (
        x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_size_t(x.shape[0]),
    )

def ctypes_matrix(x, with_stride=False):
    '''Marshal a matrix into C.'''
    import ctypes, numpy
    if x is None:
        return (
            (ctypes.POINTER(ctypes.c_double)(),) +
            ((ctypes.c_size_t(0),) if with_stride else ()) +
            (ctypes.c_size_t(0),
             ctypes.c_size_t(0))
        )
    assert x.dtype in (float, numpy.float64, complex, numpy.complex128)
    assert len(x.shape) == 2
    # if stride isn't specified then it must have a natural stride
    assert with_stride or x.strides[0] == x.shape[1] * x.itemsize
    assert x.strides[1] == x.itemsize
    assert x.flags["C_CONTIGUOUS"]
    return (
        (x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),) +
        ((ctypes.c_size_t(x.strides[0] // x.itemsize),)
         if with_stride else ()) +
        (ctypes.c_size_t(x.shape[0]),
         ctypes.c_size_t(x.shape[1]))
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

def calc_V_matrix_chunk(args):
    import numpy
    k1, k2, l, ws, R_max, abs_err, rel_err, limit = args
    len_k2 = len(k2) if k2 is not None else len(k1)
    V = numpy.empty((len(k1), len_k2), dtype=complex)
    c_generate = libproj.generate_v_matrix
    c_generate.restype = ctypes.c_uint
    num_eval = run_interruptibly(lambda: c_generate(*(
        ctypes_matrix(V, with_stride=True)[:2] +
        ctypes_vector(k1) +
        ctypes_vector(k2) +
        (ctypes.c_uint(l),
         ws,
         ctypes.c_double(R_max),
         ctypes.c_double(abs_err),
         ctypes.c_double(rel_err),
         ctypes.c_uint(limit))
    )))
    return (V, num_eval)

def calc_V_matrix(pool=None, chunk_size=4):
    '''Calculate the V matrix by subdividing it into chunks and running the
    calculation in parallel.  The `chunk_size` argument determines the minimum
    number of rows or columns per chunk.'''
    if pool:
        pmap = pool.imap
    else:
        pmap = map
    @cached(__file__ + ".cache/v-matrix/{0}.npy", np.load, np.save)
    def self(k, l, ws, R_max, abs_err, rel_err, limit):
        import itertools, timeit, sys
        import numpy

        # start
        sys.stdout.write("generating V matrix... ")
        sys.stdout.flush()
        t0 = timeit.default_timer()

        # build the chunk indices
        len_k = len(k)
        if not len_k:  # we don't support zero sizes so handle them here first
            return numpy.empty((0, 0), dtype=complex)
        indices = numpy.array_split(numpy.arange(len_k),
                                    max(1, len_k // chunk_size))
        biindices = tuple(itertools.combinations_with_replacement(indices, 2))

        # run calculations in process pool
        result = pmap(calc_V_matrix_chunk,
                      ((k[i1], k[i2] # include k[i2] only on off-diagonal chunks
                        if i1[0] != i2[0] else None,
                        l, ws, R_max, abs_err, rel_err, limit)
                       for (i1, i2) in biindices))

        # collect results, making use of symmetry
        V = numpy.empty((len_k, len_k), dtype=complex)
        total_num_eval = 0
        for (i1, i2), (subV, num_eval) in zip(biindices, result):
            total_num_eval += num_eval
            V[i1[0]:(i1[-1] + 1),
              i2[0]:(i2[-1] + 1)] = subV
            if not i1[0] == i2[0]:      # off-diagonal
                V[i2[0]:(i2[-1] + 1),
                  i1[0]:(i1[-1] + 1)] = subV.T

        # finish
        print("done ({0:.4} s, {1} V(r) evaluations)."
              .format(timeit.default_timer() - t0, total_num_eval))
        return V

    return self

def plot_bessel():
    import numpy as np
    from numpy import abs
    R = np.linspace(0, 20, 150)
    fig, axes = simple_figure()
    axes.plot(R, abs(R * bessel_j(l, R)) ** 2)
    axes.plot(R, abs(R * bessel_j(l, (.8 + 0.1j) * R)) ** 2)
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
          node_generator=gauss_legendre_nodes, pool=None):
    import numpy as np
    from numpy import sqrt
    k, w = triangle_curve(k_start, k_max, ka, kb, count, node_generator)
    Vm = calc_V_matrix(pool=pool)(k, l, ws, R_max, abs_err, rel_err, gk_limit)
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

def simple_run(pool=None):
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
        pool     = pool,
    )

    print("energies:")
    print(Es)

    eig_k = sqrt(mu2 * Es + 0j)
    fig1, axes = simple_figure()
    axes.scatter(k.real, k.imag, 50, "blue", linewidth=0)
    axes.scatter(eig_k.real, eig_k.imag, 50, "red", linewidth=0)

    R = np.linspace(0.001, 200, 2000)
    Jm = bessel_j(l, np.outer(R, k))

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

def convergence_run1(pool=None):
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
            pool     = pool,
        )
        ek = sqrt(mu2 * Es + 0j)
        axes.scatter(ek.real, ek.imag, 50, color=c,
                     linewidth=0, label=repr(count))
    axes.legend()
    return fig

def convergence_run2(pool=None):
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
            pool     = pool,
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

pool = multiprocessing.Pool()
print("note: using {0} threads.".format(multiprocessing.cpu_count()))

simple_figs = simple_run(pool=pool)
convergence_fig1 = convergence_run1(pool=pool)
convergence_fig2 = convergence_run2(pool=pool)
plt.show()
