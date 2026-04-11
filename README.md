# Colleague Rootfinders

A Fortran library for finding all roots of functions in a specified region. Includes rootfinders for real functions on intervals, complex analytic functions near real intervals, and complex analytic functions inside a square in the complex plane.

## Algorithms

The library implements two complementary approaches:

**Chebyshev expansion on intervals** ([arXiv:2102.12186](https://arxiv.org/abs/2102.12186)): Expands the function in a Chebyshev basis of order n on [a,b]. The three-term recurrence of the Chebyshev polynomials gives rise to the colleague matrix — a Hermitian tridiagonal matrix plus a rank-one correction — whose eigenvalues are the roots of the approximating polynomial. By exploiting the special structure of the matrix, these are computed using an O(n²) QR eigensolver with componentwise backward stability. Works for real or complex-valued functions on real intervals. Adaptive binary subdivision for functions requiring high resolution. Modified from code written by Kirill Serkh.

**Polynomial expansion on squares** ([arXiv:2307.14494](https://arxiv.org/abs/2307.14494)): Expands a complex analytic function on the boundary of a square in a basis of order n consisting of special polynomials constructed from an unconjugated inner product with random weights, which satisfy a three-term recurrence while remaining well-conditioned. The three-term recurrence gives rise to a generalized colleague matrix — a complex symmetric tridiagonal matrix plus a rank-one correction — whose eigenvalues are the roots of the approximating polynomial. As in the Chebyshev case above, these are computed using an O(n²) QR algorithm, here with complex orthogonal rotations. Adaptive quadtree subdivision for large regions with many roots. Optional residue-based root counting for verification.

The Chebyshev method is provably componentwise backward stable. The square method achieves machine precision experimentally. Both support optional Newton refinement.

## Main subroutines

### Real roots on [a,b]

```fortran
call droots_cheb(ifnewton, fun, par1, par2, a, b, n, eps,
     1    roots, nroots, errest, ier)

call droots_cheb_adap(ifnewton, fun, par1, par2, a, b, n, eps, nexdp,
     1    roots, nroots, errest, ier)
```

- `fun(x, par1, par2, val, dval)` — real function; derivative `dval` is only used when Newton refinement is enabled (`ifnewton=1`)
- `n` — Chebyshev expansion order (recommended 40-200)
- No LAPACK dependency

### Complex roots near [a,b]

```fortran
call zroots_cheb(ifnewton, fun, par1, par2, a, b, n, eps, delta,
     1    croots, nroots, errest, ier)

call zroots_cheb_adap(ifnewton, fun, par1, par2,
     1    a, b, n, eps, delta, nexdp, croots, nroots, errest, ier)
```

- `fun(z, par1, par2, val, dval)` — complex analytic function; derivative `dval` is only used when Newton refinement is enabled (`ifnewton=1`)
- `delta` — semi-minor axis of Bernstein ellipse for root filtering
- Finds complex roots in a strip near [a,b]
- No LAPACK dependency

### Complex roots in a square

```fortran
call zrootsq(ifprint, ifnewton, ifres, fun, par1, par2,
     1    norder, eps, center, sqw, nrtot, errest, roots, ier)

call zrootsq_adap(ifprint, ifnewton, ifres, fun, par1, par2,
     1    norder, eps, nexdp, center, sqw,
     2    nrtot, errest, roots, centers, nc, ier)
```

- `fun(z, par1, par2, val, dval)` — complex analytic function; derivative `dval` is only used when Newton refinement (`ifnewton=1`) or residue-based root counting (`ifres=1`) is enabled
- `norder` — expansion order (max 100, recommended 40 for adaptive version)
- `center` — center of search square (complex)
- `sqw` — side length of search square
- Residue-based root counting for verification (`ifres=1`)
- Requires LAPACK (`zgeqrf`, `zunmqr`, `ztrtrs`)

## Library files

| File | Description | Dependencies |
|------|-------------|-------------|
| `src/droots_cheb.f` | Real rootfinder on [a,b] | `herm_p_rank1.f` |
| `src/zroots_cheb.f` | Complex rootfinder near [a,b] | `herm_p_rank1.f` |
| `src/zrootsq.f` | Complex rootfinder on squares via generalized colleague matrix | `csym_p_rank1.f`, LAPACK |
| `src/herm_p_rank1.f` | O(n²) eigensolver for Hermitian tridiagonal + rank-1 (colleague matrix) | none |
| `src/csym_p_rank1.f` | O(n²) eigensolver for complex symmetric tridiagonal + rank-1 (generalized colleague matrix) | none |

## Build

```bash
cd test
bash test_droots_cheb       # real rootfinder
bash test_droots_cheb_adap  # adaptive real rootfinder
bash test_zroots_cheb       # complex rootfinder on interval
bash test_zroots_cheb_adap  # adaptive complex rootfinder on interval
bash test_zrootsq           # complex rootfinder on square
bash test_zrootsq_adap      # adaptive complex rootfinder on square
bash test_herm_p_rank1      # Hermitian eigensolver test
```

Requires `gfortran` and LAPACK (via Accelerate on macOS or `-llapack` on Linux). The Chebyshev rootfinders (`droots_cheb`, `zroots_cheb`) do not require LAPACK.

## References

- K. Serkh and V. Rokhlin, "A provably componentwise backward stable O(n²) QR algorithm for the diagonalization of colleague matrices." [arXiv:2102.12186](https://arxiv.org/abs/2102.12186)
- H. Zhang and V. Rokhlin, "Finding roots of complex analytic functions via generalized colleague matrices." [arXiv:2307.14494](https://arxiv.org/abs/2307.14494)
