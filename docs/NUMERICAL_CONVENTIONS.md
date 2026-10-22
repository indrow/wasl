# Numerical Conventions

## Conservation law

Wasl uses

```text
partial_t U + div(F(U)) = S(U).
```

The numerical flux is `F*(U_L, U_R, n_LR)`.

## Reference line

The reference line is `xi in [-1, 1]`, oriented from local node 0 to local node
1. Use an orthonormal Legendre modal basis. Mode `j` has polynomial degree `j`,
and a Pp line has `p + 1` modes. Mode zero represents the cell average.

Interior points have two adjacent elements. Boundary points have one adjacent
element and an invalid right element. The canonical normal points from the left
element to the right element; boundary normals point outward.

## Reference triangle (after 1D)

Vertices are ordered:

```text
v0 = (-1, -1)
v1 = ( 1, -1)
v2 = (-1,  1)
```

Element connectivity must be normalized to positive physical orientation.

## Faces

- A face key is the sorted pair of its global vertex IDs.
- Interior faces have exactly two adjacent elements.
- Boundary faces have one adjacent element and an invalid right element.
- The canonical normal points from left element to right element.
- Boundary normals point outward from the left element.
- Left/right ownership is deterministic and independent of MSH record order.
- Face permutations map neighbor traces into the canonical face coordinate.

## Modal basis

- Use an orthonormal Dubiner-type basis on the reference triangle.
- Triangle mode count is `(p + 1)(p + 2)/2`.
- Mode zero is the cell-average mode.
- Modes are ordered by total degree, then by the first reference-coordinate
  exponent.
- Use one polynomial order and element type per launch.

## State layout

Persisted coefficients are structure-of-arrays:

```text
U[variable][mode][element]
index = ((variable * mode_count) + mode) * element_count + element
```

Euler 2D conserved variables are:

```text
U[0] = density
U[1] = x-momentum
U[2] = y-momentum
U[3] = total energy
```

## Precision

- `double` is the reference precision.
- Every regression tolerance belongs with its benchmark, not in a global
  catch-all constant.

## Invariants

The relevant tests must preserve constant states, shared-face conservation,
free streams, cell average under WENO limiting, positive supported Euler
states, expected smooth convergence, and serial/partitioned agreement.
