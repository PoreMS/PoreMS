################################################################################
# Geometry                                                                     #
#                                                                              #
"""Here basic geometric functions are noted."""
################################################################################


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
    return float(np.dot(vec_a, vec_b))


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
    return float(np.linalg.norm(vec))


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
    vector : numpy.ndarray
        Bond vector
    """
    a = np.asarray(pos_a)
    b = np.asarray(pos_b)

    if a.shape != b.shape:
        print("Vector: Wrong dimensions...")
        return None

    return b - a


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
    vec : numpy.ndarray
        Unit vector
    """
    v = np.asarray(vec, dtype=float)
    n = np.linalg.norm(v)
    return v / n if n != 0 else v


def cross_product(vec_a, vec_b):
    r"""Calculate the cross product of two three-dimensional vectors
    :math:`\boldsymbol{a},\boldsymbol{b}\in\mathbb{R}^3`

    .. math::

        \text{cross}(\boldsymbol{a},\boldsymbol{b})=\begin{pmatrix}
        a_2\cdot b_3-a_3\cdot b_2\\
        a_3\cdot b_1-a_1\cdot b_4\\
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
    vec : numpy.ndarray
        Cross product vector
    """
    return np.cross(vec_a, vec_b)


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
    cos_val = np.dot(vec_a, vec_b) / (np.linalg.norm(vec_a) * np.linalg.norm(vec_b))
    a = np.arccos(np.clip(cos_val, -1.0, 1.0))
    return float(np.degrees(a)) if is_deg else float(a)


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
    a = np.arctan2(pos[1], pos[0])
    return float(np.degrees(a)) if is_deg else float(a)


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
    n = float(np.linalg.norm(pos))
    a = np.arccos(pos[2] / n) if n != 0 else np.arccos(0)
    return float(np.degrees(a)) if is_deg else float(a)


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
    vec : numpy.ndarray
        Unit vector
    """
    axis_error = "Wrong axis definition..."

    if isinstance(inp, str):
        mapping = {"x": 1, "y": 2, "z": 3}
        if inp not in mapping:
            return axis_error
        axis = mapping[inp]
    elif isinstance(inp, int):
        if inp not in (1, 2, 3):
            return axis_error
        axis = inp
    else:
        return axis_error

    v = np.zeros(dim)
    v[axis - 1] = 1.0
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
    coord : numpy.ndarray
        Vector c as the result of the rotation
    """
    angle = np.radians(angle) if is_deg else float(angle)

    if isinstance(axis, np.ndarray):
        # Already a vector — use directly
        n = axis.astype(float)
    elif isinstance(axis, list):
        if len(axis) == dim:
            n = np.asarray(axis, dtype=float)
        elif len(axis) == 2:
            v = vector(axis[0], axis[1])
            if v is None:
                print("Rotate: Wrong vector dimensions.")
                return None
            n = np.asarray(v, dtype=float)
        else:
            print("Rotate: Wrong vector dimensions.")
            return None
    else:
        n = main_axis(axis)
        if isinstance(n, str):
            print("Rotate: " + n)
            return None
        n = np.asarray(n, dtype=float)

    n_norm = np.linalg.norm(n)
    if n_norm > 0:
        n = n / n_norm
    n1, n2, n3 = n[0], n[1], n[2]
    c = np.cos(angle)
    s = np.sin(angle)

    R = np.array([
        [n1*n1*(1-c)+c,    n1*n2*(1-c)-n3*s, n1*n3*(1-c)+n2*s],
        [n2*n1*(1-c)+n3*s, n2*n2*(1-c)+c,    n2*n3*(1-c)-n1*s],
        [n3*n1*(1-c)-n2*s, n3*n2*(1-c)+n1*s, n3*n3*(1-c)+c   ]
    ])

    # For homogeneous arrays use einsum (fast, handles (3,) and (3,N,M))
    # For inhomogeneous plotting data fall back to element-wise broadcasting
    try:
        return np.einsum("ij,j...->i...", R, np.asarray(data, dtype=float))
    except (ValueError, TypeError):
        return [R[i, 0]*data[0] + R[i, 1]*data[1] + R[i, 2]*data[2] for i in range(3)]
