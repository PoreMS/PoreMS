################################################################################
# Geometry                                                                     #
#                                                                              #
"""Here basic geometric functions are noted."""
################################################################################


import math
import numpy as np


def dot_product(vec_a, vec_b):
    r"""Calculate the dot product of two vectors
    :math:`\boldsymbol{a},\boldsymbol{b}\in\mathbb{R}^n`

    .. math::

        \text{dot}(\boldsymbol{a},\boldsymbol{b})=
        \begin{pmatrix}a_1\\\vdots\\a_n\end{pmatrix}\cdot
        \begin{pmatrix}b_1\\\vdots\\b_n\end{pmatrix}=
        a_1\cdot b_1+a_2\cdot b_2+\dots+a_n\cdot b_n.

    Parameters
    ----------
    vec_a : list
        First vector :math:`\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\boldsymbol{b}`

    Returns
    -------
    dot : float
        Dot product value
    """
    return sum(ai*bi for ai, bi in zip(vec_a, vec_b))


def length(vec):
    r"""Calculate the length of a vector
    :math:`\boldsymbol{a}\in\mathbb{R}^n`

    .. math::

        \text{length}(\boldsymbol{a})=|\boldsymbol{a}|
        =\sqrt{\boldsymbol{a}\cdot\boldsymbol{a}}

    Parameters
    ----------
    vec : list
        Vector a

    Returns
    -------
    length : float
        Vector length
    """
    return math.sqrt(sum(x*x for x in vec))


def vector(pos_a, pos_b):
    r"""Calculate the vector between to two positions
    :math:`\boldsymbol{a},\boldsymbol{b}\in\mathbb{R}^n`

    .. math::

        \text{vec}(\boldsymbol{a},\boldsymbol{b})
        =\begin{pmatrix}b_1-a_1\\\vdots\\b_n-a_n\end{pmatrix}

    Parameters
    ----------
    pos_a : list
        First position :math:`\boldsymbol{a}`
    pos_b : list
        Second position :math:`\boldsymbol{b}`

    Returns
    -------
    vector : list
        Bond vector
    """
    if len(pos_a) != len(pos_b):
        print("Vector: Wrong dimensions...")
        return None
    return [b - a for a, b in zip(pos_a, pos_b)]


def unit(vec):
    r"""Transform a vector :math:`\boldsymbol{a}\in\mathbb{R}^n` into a
    unit vector

    .. math::

        \text{unit}(\boldsymbol{a})
        =\frac{\boldsymbol{a}}{|\boldsymbol{a}|}

    Parameters
    ----------
    vec : list
        Vector a

    Returns
    -------
    vec : list
        Unit vector
    """
    n = math.sqrt(sum(x*x for x in vec))
    if n == 0:
        return [0.0] * len(vec)
    return [x/n for x in vec]


def cross_product(vec_a, vec_b):
    r"""Calculate the cross product of two three-dimensional vectors
    :math:`\boldsymbol{a},\boldsymbol{b}\in\mathbb{R}^3`

    .. math::

        \text{cross}(\boldsymbol{a},\boldsymbol{b})=\begin{pmatrix}
        a_2\cdot b_3-a_3\cdot b_2\\
        a_3\cdot b_1-a_1\cdot b_3\\
        a_1\cdot b_2-a_2\cdot b_1
        \end{pmatrix}

    Parameters
    ----------
    vec_a : list
        First vector :math:`\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\boldsymbol{b}`

    Returns
    -------
    vec : list
        Cross product vector
    """
    a, b = vec_a, vec_b
    return [a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]]


def angle(vec_a, vec_b, is_deg=True):
    r"""Calculate the angle between two vectors
    :math:`\boldsymbol{a},\boldsymbol{b}\in\mathbb{R}^n`

    .. math::

        \text{angle}=\cos^{-1}\frac{\boldsymbol{a}\cdot\boldsymbol{b}}
        {|\boldsymbol{a}||\boldsymbol{a}|}

    Parameters
    ----------
    vec_a : list
        First vector :math:`\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\boldsymbol{b}`
    is_deg : bool, optional
        True if the output should be in degree

    Returns
    -------
    angle : float
        Angle
    """
    dot = sum(ai*bi for ai, bi in zip(vec_a, vec_b))
    la = math.sqrt(sum(x*x for x in vec_a))
    lb = math.sqrt(sum(x*x for x in vec_b))
    cos_val = max(-1.0, min(1.0, dot / (la * lb)))
    a_rad = math.acos(cos_val)
    return math.degrees(a_rad) if is_deg else a_rad


def angle_polar(pos, is_deg=False):
    r"""Calculate the polar angle of a position vector
    :math:`\boldsymbol{a}\in\mathbb{R}^3`, which is the angle of the
    x-axis towards the reflected position vector on the x-y-plane

    .. math::

        \text{polar}(\boldsymbol{a})=\arctan2(x,y)\left\{
        \begin{array}{ll}
        \tan^{-1}\left(\frac{y}{x}\right)&x>0\\
        \tan^{-1}\left(\frac{y}{x}\right)+\pi&x<0,y>0\\
        \pm\pi&x<0,y=0\\
        \tan^{-1}\left(\frac{y}{x}\right)-\pi&x<0,y<0\\
        +\frac{\pi}{2}&x=0,y>0\\
        -\frac{\pi}{2}&x=0,y<0
        \end{array}
        \right.

    with :math:`x` as the first vector entry and :math:`y` as the second.

    Parameters
    ----------
    pos : list
        Position vector :math:`\boldsymbol{a}`
    is_deg : bool, optional
        True if the output should be in degree

    Returns
    -------
    angle : float
        Polar angle
    """
    a = math.atan2(pos[1], pos[0])
    return math.degrees(a) if is_deg else a


def angle_azi(pos, is_deg=False):
    r"""Calculate the azimuthal angle of a position vector
    :math:`\boldsymbol{a}\in\mathbb{R}^3`, which is the angle of the
    position vector towards the x-y-plane

    .. math::

        \text{azimut}(\boldsymbol{a})
        =\cos^{-1}\frac{y}{|\boldsymbol{a}|}

    with :math:`y` as the second vector entry.

    Parameters
    ----------
    pos : list
        Position vector :math:`\boldsymbol{a}`
    is_deg : bool, optional
        True if the output should be in degree

    Returns
    -------
    angle : float
        Azimuthal angle
    """
    v = pos
    n = math.sqrt(sum(x*x for x in v))
    a = math.acos(max(-1.0, min(1.0, v[2] / n))) if n != 0 else math.acos(0)
    return math.degrees(a) if is_deg else a


def main_axis(inp, dim=3):
    """Return the three-dimensional unit-vector of the main axes.
    Input is either integer or string

    * 1 or "x" for the x-axis
    * 2 or "y" for the y-axis
    * 3 or "z" for the z-axis

    Parameters
    ----------
    inp : integer, string
        Axis type input
    dim : integer, optional
        Number of dimensions

    Returns
    -------
    vec : list
        Unit vector
    """
    axis_error = "Wrong axis definition..."

    if isinstance(inp, str):
        mapping = {"x": 0, "y": 1, "z": 2}
        if inp not in mapping:
            return axis_error
        idx = mapping[inp]
    elif isinstance(inp, int):
        if inp not in (1, 2, 3):
            return axis_error
        idx = inp - 1
    else:
        return axis_error

    v = [0.0] * dim
    v[idx] = 1.0
    return v


def rotate(data, axis, angle, is_deg, dim=3):
    r"""Rotate a vector :math:`\boldsymbol{a}\in\mathbb{R}^3`
    along an axis :math:`\boldsymbol{b}\in\mathbb{R}^3` with angle
    :math:`\alpha\in\mathbb{R}`. The rotation is performed using the
    rotation-matrix

    .. math::

        \boldsymbol{R}_\boldsymbol{n}(\alpha)=\begin{pmatrix}
        n_1^2 (1-\cos\alpha)+   \cos\alpha&n_1n_2(1-\cos\alpha)-n_3\sin\alpha&n_1n_3(1-\cos\alpha)+n_2\sin\alpha\\
        n_2n_1(1-\cos\alpha)+n_3\sin\alpha&n_2^2 (1-\cos\alpha)+   \cos\alpha&n_2n_3(1-\cos\alpha)-n_1\sin\alpha\\
        n_3n_1(1-\cos\alpha)-n_2\sin\alpha&n_3n_2(1-\cos\alpha)+n_1\sin\alpha&n_3^2 (1-\cos\alpha)+   \cos\alpha
        \end{pmatrix}

    where :math:`n_i` are the entries for the unit vector
    :math:`\boldsymbol{n}` of the axis. The new coordinates
    :math:`\boldsymbol{c}` are then calculated using a matrix-vector
    multiplication

    .. math::

        \boldsymbol{c}=\boldsymbol{R}_\boldsymbol{n}\boldsymbol{a}.

    Parameters
    ----------
    data : list
        Vector :math:`\boldsymbol{a}`
    axis : integer, string, list
        Axis :math:`\boldsymbol{b}`
    angle : float
        Angle
    is_deg : bool
        True if the input is in degree
    dim : integer, optional
        Number of dimensions

    Returns
    -------
    coord : list or numpy.ndarray
        Vector c as the result of the rotation (list for single 3D vectors,
        numpy array for multi-dimensional data)
    """
    angle = math.radians(angle) if is_deg else float(angle)

    if isinstance(axis, np.ndarray):
        nx, ny, nz = float(axis[0]), float(axis[1]), float(axis[2])
        nn = math.sqrt(nx*nx + ny*ny + nz*nz)
        if nn > 0:
            nx, ny, nz = nx/nn, ny/nn, nz/nn
    elif isinstance(axis, list):
        if len(axis) == dim:
            n = axis
            nn = math.sqrt(n[0]**2 + n[1]**2 + n[2]**2)
            nx, ny, nz = (n[0]/nn, n[1]/nn, n[2]/nn) if nn > 0 else (n[0], n[1], n[2])
        elif len(axis) == 2:
            v = vector(axis[0], axis[1])
            nn = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
            nx, ny, nz = (v[0]/nn, v[1]/nn, v[2]/nn) if nn > 0 else (v[0], v[1], v[2])
        else:
            print("Rotate: Wrong vector dimensions.")
            return None
    else:
        n = main_axis(axis)
        if isinstance(n, str):
            print("Rotate: " + n)
            return None
        nx, ny, nz = n[0], n[1], n[2]

    c = math.cos(angle)
    s = math.sin(angle)
    ic = 1.0 - c

    # Inline rotation matrix coefficients
    r00 = nx*nx*ic + c;    r01 = nx*ny*ic - nz*s; r02 = nx*nz*ic + ny*s
    r10 = ny*nx*ic + nz*s; r11 = ny*ny*ic + c;    r12 = ny*nz*ic - nx*s
    r20 = nz*nx*ic - ny*s; r21 = nz*ny*ic + nx*s; r22 = nz*nz*ic + c

    # Fast path: single 3D vector (most common in pore generation)
    if isinstance(data, (list, tuple)) and len(data) == dim:
        d0, d1, d2 = data[0], data[1], data[2]
        return [r00*d0 + r01*d1 + r02*d2,
                r10*d0 + r11*d1 + r12*d2,
                r20*d0 + r21*d1 + r22*d2]

    # Array path: shape plotting data
    R = np.array([[r00, r01, r02],
                  [r10, r11, r12],
                  [r20, r21, r22]])
    try:
        return np.einsum("ij,j...->i...", R, np.asarray(data, dtype=float))
    except (ValueError, TypeError):
        d = data
        return [r00*d[0] + r01*d[1] + r02*d[2],
                r10*d[0] + r11*d[1] + r12*d[2],
                r20*d[0] + r21*d[1] + r22*d[2]]
