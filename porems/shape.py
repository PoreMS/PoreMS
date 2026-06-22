################################################################################
# Shape Pack                                                                   #
#                                                                              #
"""This file contains shape definitions to be cut out from the crystal block."""
################################################################################


import math
import numpy as np
import matplotlib.pyplot as plt

import porems.utils as utils
import porems.geometry as geometry


class Shape():
    """This class is a container for individual shape classes.

    Parameters
    ----------
    inp : dictionary
        Dictionary of necessary inputs
    """
    def __init__(self, inp):
        self._inp = inp

        # Calculate angle and normal vector for rotation
        self._angle = geometry.angle(inp["central"], geometry.main_axis("z"))
        self._normal = geometry.cross_product(inp["central"], geometry.main_axis("z"))

        # Calculate distance towards central axis start
        self._dist_start = geometry.vector(self._centroid, inp["centroid"])
        self._dist_zero = geometry.vector(geometry.rotate(inp["centroid"], self._normal, -self._angle, True), self._centroid)


    ##################
    # Helper Methods #
    ##################
    def convert(self, data, to_zero=True):
        """This helper method rotates the given data to match the global central
        axis and translates it so that the center point is aligned.

        Parameters
        ----------
        data : list
            Data input
        to_zero : bool
            True to convert data towards zero axis, False to convert data from
            zero axis to central axis.

        Returns
        -------
        data : list
            Converted input
        """
        # Rotate towards main axis to the zero axis
        data = geometry.rotate(data, self._normal, -self._angle if to_zero else self._angle, True)

        # Translate to zero or to start
        dist = self._dist_zero if to_zero else self._dist_start
        data = [data[i]+dist[i] for i in range(3)]

        return data

    def plot(self, inp=0, num=100, vec=[]):
        """Plot surface and rim.

        Parameters
        ----------
        inp : float, optional
            Position on the axis
        num : integer, optional
            Number of points
        vec : list, optional
            Vector on surface to test normal vector
        """
        fig = plt.gcf()
        ax = fig.add_subplot(111, projection="3d")

        # Surface
        ax.plot_surface(*self.surf(num=100), alpha=0.7)

        # Rim
        ax.plot3D(*[x[0] for x in self.rim(inp, num)])

        # Normal
        if vec:
            line = [self.convert([0, 0, 0], False),
                    vec,
                    self.convert(self.normal(vec), False)]
            ax.plot3D(*utils.column(line))


    ##########
    # Getter #
    ##########
    def get_inp(self):
        """Return full input dictionary.

        Returns
        -------
        inp : dictionary
            Dictionary of all inputs
        """
        return self._inp


class Cylinder(Shape):
    """This class defines a cylindrical shape. Needed inputs are

    * **central** - Central axis
    * **centroid** - Centroid of block
    * **length** - Cylinder length
    * **diameter** - Cylinder diameter

    Parameters
    ----------
    inp : dictionary
        Dictionary of necessary inputs
    """
    def __init__(self, inp):
        # Set centroid
        self._centroid = [0, 0, inp["length"]/2]

        # Call super class
        super(Cylinder, self).__init__(inp)


    ############
    # Function #
    ############
    def Phi(self, r, phi, z):
        """Surface function of a cylinder

        .. math::

            \\Phi(r,\\phi,z)=
            \\begin{bmatrix}r\\cos(\\phi)\\\\r\\sin(\\phi)\\\\z\\end{bmatrix}

        with radius :math:`r`, polar angle :math:`\\phi` and cylinder length
        :math:`z`.

        Parameters
        ----------
        r : float
            Radius
        phi : float
            Polar angle
        z : float
            Distance ins z-axis

        Returns
        -------
        pos : list
            Cartesian coordinates for given polar coordinates
        """
        x = np.outer(r, np.cos(phi))
        y = np.outer(r, np.sin(phi))
        z = np.outer(z, np.ones(len(z)))

        return self.convert([x, y, z], False)

    def d_Phi_phi(self, r, phi, z):
        """Derivative of the surface function considering the polar angle

        .. math::

            \\frac{\\partial\\Phi}{\\partial\\phi}(r,\\phi,z)=
            \\begin{bmatrix}-r\\sin(\\phi)\\\\r\\cos(\\phi)\\\\0\\end{bmatrix}

        with radius :math:`r`, polar angle :math:`\\phi` and cylinder length
        :math:`z`.

        Parameters
        ----------
        r : float
            Radius
        phi : float
            Polar angle
        z : float
            Distance ins z-axis

        Returns
        -------
        pos : list
            Cartesian coordinates for given polar coordinates
        """
        x = -r*np.sin(phi)
        y = r*np.cos(phi)
        z = 0

        return [x, y, z]

    def d_Phi_z(self, r, phi, z):
        """Derivative of the surface function considering the z-axis

        .. math::

            \\frac{\\partial\\Phi}{\\partial z}(r,\\phi,z)=
            \\begin{bmatrix}0\\\\0\\\\1\\end{bmatrix}

        with radius :math:`r`, polar angle :math:`\\phi` and cylinder length
        :math:`z`.

        Parameters
        ----------
        r : float
            Radius
        phi : float
            Polar angle
        z : float
            Distance ins z-axis

        Returns
        -------
        pos : list
            Cartesian coordinates for given polar coordinates
        """
        x = 0
        y = 0
        z = 1

        return [x, y, z]


    ############
    # Features #
    ############
    def normal(self, pos):
        """Calculate unit normal vector on surface for a given position

        .. math::

            \\frac{\\partial\\Phi}{\\partial\\phi}(r,\\phi,z)\\times
            \\frac{\\partial\\Phi}{\\partial z}(r,\\phi,z)=
            \\begin{bmatrix}r\\cos(\\phi)\\\\r\\sin(\\phi)\\\\0\\end{bmatrix}

        with radius :math:`r`, polar angle :math:`\\phi` and cylinder length
        :math:`z`.

        Parameters
        ----------
        pos : list
            Position

        Returns
        -------
        normal : list
            Normal vector
        """
        # Initialize
        x, y, z = self.convert(pos)

        # Cartesian to polar
        r = math.sqrt(x**2+y**2)
        phi = geometry.angle_polar([x, y, z])

        # Calculate derivatives
        d_Phi_phi = self.d_Phi_phi(r, phi, z)
        d_Phi_z = self.d_Phi_z(r, phi, z)

        # Calculate normal vector
        return geometry.cross_product(d_Phi_phi, d_Phi_z)

    def is_in(self, pos):
        """Check if given position is inside of shape.

        Parameters
        ----------
        pos : list
            Position

        Returns
        -------
        is_in : bool
            True if position is inside of shape
        """
        # Check if within shape
        if geometry.length(self.normal(pos)) < self._inp["diameter"]/2:
            pos_zero = self.convert(pos)
            return pos_zero[2]>0 and pos_zero[2]<self._inp["length"]
        else:
            return False


    #########
    # Shape #
    #########
    def rim(self, z, num=100):
        """Return x and y values for given z-position.

        Parameters
        ----------
        z : float
            Position on the axis
        num : integer, optional
            Number of points

        Returns
        -------
        positions : list
            x and y arrays of the surface rim on the z-position
        """
        phi = np.linspace(0, 2*np.pi, num)
        r = self._inp["diameter"]/2

        return self.Phi(r, phi, [z])

    def surf(self, num=100):
        """Return x, y and z values for the shape.

        Parameters
        ----------
        num : integer, optional
            Number of points

        Returns
        -------
        positions : list
            x, y and z arrays of the surface rim
        """
        phi = np.linspace(0, 2*np.pi, num)
        r = np.ones(num)*self._inp["diameter"]/2
        z = np.linspace(0, self._inp["length"], num)

        return self.Phi(r, phi, z)


    ##############
    # Properties #
    ##############
    def volume(self):
        """Calculate volume

        .. math::

            V=\\pi r^2l

        with radius :math:`r` and cylinder length :math:`l`.

        Returns
        -------
        volume : float
            Volume
        """
        return math.pi*(self._inp["diameter"]/2)**2*self._inp["length"]

    def surface(self):
        """Calculate inner surface

        .. math::

            S=2\\pi rl

        with radius :math:`r` and cylinder length :math:`l`.

        Returns
        -------
        surface : float
            Inner surface
        """
        return 2*math.pi*self._inp["diameter"]/2*self._inp["length"]


class Sphere(Shape):
    """This class defines a sphere shape. Needed inputs are

    * **central** - Central axis
    * **centroid** - Sphere Centroid
    * **diameter** - Sphere diameter

    Parameters
    ----------
    inp : dictionary
        Dictionary of necessary inputs
    """
    def __init__(self, inp):
        # Set centroid
        self._centroid = [0, 0, 0]

        # Call super class
        super(Sphere, self).__init__(inp)


    ############
    # Function #
    ############
    def Phi(self, r, theta, phi):
        """Surface function of a sphere

        .. math::

            \\Phi(r,\\phi,\\theta)=
            \\begin{bmatrix}r\\cos(\\phi)\\sin(\\theta)\\\\r\\sin(\\phi)\\sin(\\theta)\\\\r\\cos(\\theta)\\end{bmatrix}

        with radius :math:`r`, polar angle :math:`\\phi` and azimuthal angle
        :math:`\\theta`.

        Parameters
        ----------
        r : float
            Radius
        theta : float
            Azimuth angle
        phi : float
            Polar angle

        Returns
        -------
        pos : list
            Cartesian coordinates for given spherical coordinates
        """
        x = r*np.outer(np.cos(phi), np.sin(theta))
        y = r*np.outer(np.sin(phi), np.sin(theta))
        z = r*np.outer(np.ones(len(phi)), np.cos(theta))

        return self.convert([x, y, z], False)

    def d_Phi_phi(self, r, theta, phi):
        """Derivative of the surface function considering the polar angle
        :math:`\\phi`

        .. math::

            \\frac{\\partial\\Phi}{\\partial\\phi}(r,\\phi,\\theta)=
            \\begin{bmatrix}-r\\sin(\\phi)\\sin(\\theta)\\\\r\\cos(\\phi)\\sin(\\theta)\\\\0\\end{bmatrix}

        with radius :math:`r`, polar angle :math:`\\phi` and azimuthal angle
        :math:`\\theta`.

        Parameters
        ----------
        r : float
            Radius
        theta : float
            Azimuth angle
        phi : float
            Polar angle

        Returns
        -------
        pos : list
            Cartesian coordinates for given spherical coordinates
        """
        x = -r*np.sin(phi)*np.sin(theta)
        y = r*np.cos(phi)*np.sin(theta)
        z = 0

        return [x, y, z]

    def d_Phi_theta(self, r, theta, phi):
        """Derivative of the surface function considering the azimuthal angle
        :math:`\\theta`

        .. math::

            \\frac{\\partial\\Phi}{\\partial\\theta}(r,\\phi,\\theta)=
            \\begin{bmatrix}r\\cos(\\phi)\\cos(\\theta)\\\\r\\sin(\\phi)\\cos(\\theta)\\\\-r\\sin(\\theta)\\end{bmatrix}

        with radius :math:`r`, polar angle :math:`\\phi` and azimuthal angle
        :math:`\\theta`.

        Parameters
        ----------
        r : float
            Radius
        theta : float
            Azimuth angle
        phi : float
            Polar angle

        Returns
        -------
        pos : list
            Cartesian coordinates for given spherical coordinates
        """
        x = r*np.cos(phi)*np.cos(theta)
        y = r*np.sin(phi)*np.cos(theta)
        z = -r*np.sin(theta)

        return [x, y, z]


    ############
    # Features #
    ############
    def normal(self, pos):
        """Calculate unit normal vector on surface for a given position

        .. math::

            \\frac{\\partial\\Phi}{\\partial\\theta}(r,\\phi,\\theta)\\times
            \\frac{\\partial\\Phi}{\\partial\\phi}(r,\\phi,\\theta)=
            \\begin{bmatrix}
            -r^2\\cos(\\phi)\\sin(\\theta)^2\\\\
            r^2\\sin(\\phi)\\sin(\\theta)^2\\\\
            -r^2\\sin(\\theta)\\cos(\\theta)\\left[\\sin(\\phi)^2-\\cos(\\phi)^2\\right]\\\\
            \\end{bmatrix}

        with radius :math:`r`, polar angle :math:`\\phi` and cylinder length
        :math:`z`.

        Parameters
        ----------
        pos : list
            Position

        Returns
        -------
        normal : list
            Normal vector
        """
        # Initialize
        x, y, z = self.convert(pos)

        # Cartesian to polar
        r = math.sqrt(x**2+y**2+z**2)
        theta = geometry.angle_azi([x, y, z])
        phi = geometry.angle_polar([x, y, z])

        # Calculate derivatives
        d_Phi_theta = self.d_Phi_theta(r, theta, phi)
        d_Phi_phi = self.d_Phi_phi(r, theta, phi)

        # Calculate normal vector
        return geometry.cross_product(d_Phi_theta, d_Phi_phi)

    def is_in(self, pos):
        """Check if given position is inside of shape.

        Parameters
        ----------
        pos : list
            Position

        Returns
        -------
        is_in : bool
            True if position is inside of shape
        """
        # Check if within shape
        pos_zero = self.convert(pos)
        if geometry.length(geometry.vector(self._centroid, pos_zero)) < self._inp["diameter"]/2:
            return abs(pos_zero[2])<self._inp["diameter"]/2
        else:
            return False


    #########
    # Shape #
    #########
    def rim(self, phi, num=100):
        """Return x and y values for given polar angle.

        Parameters
        ----------
        phi : float
            Position on the axis
        num : integer, optional
            Number of points

        Returns
        -------
        positions : list
            x and y arrays of the surface rim on the z-position
        """
        r = self._inp["diameter"]/2
        theta = np.linspace(0, 2*np.pi, num)

        return self.Phi(r, theta, [phi])

    def surf(self, num=100):
        """Return x, y and z values for the shape.

        Parameters
        ----------
        num : integer, optional
            Number of points

        Returns
        -------
        positions : list
            x, y and z arrays of the surface rim
        """
        r = self._inp["diameter"]/2
        theta = np.linspace(0, np.pi, num)
        phi = np.linspace(0, 2*np.pi, num)

        return self.Phi(r, theta, phi)


    ##############
    # Properties #
    ##############
    def volume(self):
        """Calculate volume

        .. math::

            V=\\frac43\\pi r^3

        with radius :math:`r`.

        Returns
        -------
        volume : float
            Volume
        """
        return 4/3*math.pi*(self._inp["diameter"]/2)**3

    def surface(self):
        """Calculate inner surface

        .. math::

            S=4\\pi r^2

        with radius :math:`r`.

        Returns
        -------
        surface : float
            Inner surface
        """
        return 4*math.pi*(self._inp["diameter"]/2)**2


class Cuboid(Shape):
    """This class defines a cuboid shape. Needed inputs are

    * **central** - Central axis
    * **centroid** - Centroid of block
    * **length** - Cuboid length
    * **width** - Cuboid width
    * **height** - Cuboid height

    Parameters
    ----------
    inp : dictionary
        Dictionary of necessary inputs
    """
    def __init__(self, inp):
        # Set centroid
        self._centroid = [inp["width"]/2, inp["height"]/2, inp["length"]/2]

        # Call super class
        super(Cuboid, self).__init__(inp)


    ############
    # Function #
    ############
    def Phi(self, x, y, z):
        """Surface function of a cuboid.

        Parameters
        ----------
        x : float
            Width
        y : float
            Height
        z : float
            Length

        Returns
        -------
        pos : list
            Cartesian coordinates for given spherical coordinates
        """
        phi = np.arange(1,10,2)*np.pi/4
        Phi, Theta = np.meshgrid(phi, phi)

        x = x*np.cos(Phi)*np.sin(Theta)
        y = y*np.sin(Phi)*np.sin(Theta)
        z = z*np.cos(Theta)/np.sqrt(2)

        return self.convert([x, y, z], False)


    ############
    # Features #
    ############
    def normal(self, pos):
        """Calculate unit normal vector on surface for a given position.

        Parameters
        ----------
        pos : list
            Position

        Returns
        -------
        normal : list
            Unit normal vector
        """
        # Initialize
        x, y, z = self.convert(pos)

        # Calculate derivatives
        return [0, -1, 0] if y < self._centroid[1] else [0, 1, 0]

    def is_in(self, pos):
        """Check if given position is inside of shape.

        Parameters
        ----------
        pos : list
            Position

        Returns
        -------
        is_in : bool
            True if position is inside of shape
        """
        pos_zero = self.convert(pos)

        return pos_zero[1] > 0 and pos_zero[1] < self._inp["height"]


    #########
    # Shape #
    #########
    def rim(self, z, num=100):
        """Return x and y values for given length.

        Parameters
        ----------
        z : float
            Position on the axis
        num : integer, optional
            Number of points

        Returns
        -------
        positions : list
            x and y arrays of the surface rim on the z-position
        """
        x = self._inp["width"]
        y = self._inp["height"]

        return self.Phi(x, y, z)

    def surf(self, num=100):
        """Return x, y and z values for the shape.

        Parameters
        ----------
        num : integer, optional
            Number of points

        Returns
        -------
        positions : list
            x, y and z arrays of the surface rim
        """
        x = self._inp["width"]
        y = self._inp["height"]
        z = self._inp["length"]

        return self.Phi(x, y, z)


    ##############
    # Properties #
    ##############
    def volume(self):
        """Calculate volume

        .. math::

            V=w\\cdot h\\cdot l

        with width :math:`w`, height :math:`h` and length :math:`l`.

        Returns
        -------
        volume : float
            Volume
        """
        return self._inp["length"]*self._inp["width"]*self._inp["height"]

    def surface(self):
        """Calculate inner surface

        .. math::

            S=2\\cdot(w\\cdot h+w\\cdot l+h\\cdot l)

        with width :math:`w`, height :math:`h` and length :math:`l`.

        Returns
        -------
        surface : float
            Inner surface
        """
        return 2*(self._inp["length"]*self._inp["width"]+self._inp["length"]*self._inp["height"]+self._inp["width"]*self._inp["height"])


class Cone(Shape):
    """This class defines a conical shape. Needed inputs are

    * **central** - Central axis
    * **centroid** - Centroid of block
    * **length** - Cone length
    * **diameter_1** - Cone starting diameter
    * **diameter_2** - Cone ending diameter

    Parameters
    ----------
    inp : dictionary
        Dictionary of necessary inputs
    """
    def __init__(self, inp):
        # Set centroid
        self._centroid = [0, 0, inp["length"]/2]

        # Call super class
        super(Cone, self).__init__(inp)


    ############
    # Function #
    ############
    def Phi(self, r, phi, z):
        """Surface function of a cone

        .. math::

            \\Phi(r(z),\\phi,z)=
            \\begin{bmatrix}r(z)\\cos(\\phi)\\\\r(z)\\sin(\\phi)\\\\z\\end{bmatrix}

        with polar angle :math:`\\phi` and radius function along the :math:`z`-axis

        .. math::

            r(z)=r_1+\\frac{r_2-r_1}{l-1}(z-1)

        with radii :math:`r_1` and :math:`r_2` and cone length :math:`l`.

        Parameters
        ----------
        r : float
            Radius
        phi : float
            Polar angle
        z : float
            Distance ins z-axis

        Returns
        -------
        pos : list
            Cartesian coordinates for given polar coordinates
        """
        def r(z):
            z = np.array(z)
            r_1 = self._inp["diameter_1"]/2
            r_2 = self._inp["diameter_2"]/2
            l = self._inp["length"]
            return r_1+(r_2-r_1)/(l-1)*(z-1)

        x = np.outer(r(z), np.cos(phi))
        y = np.outer(r(z), np.sin(phi))
        z = np.outer(z, np.ones(len(z)))

        return self.convert([x, y, z], False)

    def d_Phi_phi(self, r, phi, z):
        """Derivative of the surface function considering the polar angle

        .. math::

            \\frac{\\partial\\Phi}{\\partial\\phi}(r(z),\\phi,z)=
            \\begin{bmatrix}-r(z)\\sin(\\phi)\\\\r(z)\\cos(\\phi)\\\\0\\end{bmatrix}

        with polar angle :math:`\\phi` and radius function along the :math:`z`-axis

        .. math::

            r(z)=r_1+\\frac{r_2-r_1}{l-1}(z-1)

        with radii :math:`r_1` and :math:`r_2` and cone length :math:`l`.

        Parameters
        ----------
        r : float
            Radius
        phi : float
            Polar angle
        z : float
            Distance ins z-axis

        Returns
        -------
        pos : list
            Cartesian coordinates for given polar coordinates
        """
        def r(z):
            r_1 = self._inp["diameter_1"]/2
            r_2 = self._inp["diameter_2"]/2
            l = self._inp["length"]
            return r_1+(r_2-r_1)/(l-1)*(z-1)

        x = -r(z)*np.sin(phi)
        y = r(z)*np.cos(phi)
        z = 0

        return [x, y, z]

    def d_Phi_z(self, r, phi, z):
        """Derivative of the surface function considering the z-axis

        .. math::

            \\frac{\\partial\\Phi}{\\partial z}(r(z),\\phi,z)=
            \\begin{bmatrix}r(z)\\cos(\\phi)\\\\r(z)\\sin(\\phi)\\\\1\\end{bmatrix}

        with polar angle :math:`\\phi` and radius function along the :math:`z`-axis

        .. math::

            r(z)=\\frac{r_2-r_1}{l-1}

        with radii :math:`r_1` and :math:`r_2` and cone length :math:`l`.

        Parameters
        ----------
        r : float
            Radius
        phi : float
            Polar angle
        z : float
            Distance ins z-axis

        Returns
        -------
        pos : list
            Cartesian coordinates for given polar coordinates
        """
        def r(z):
            r_1 = self._inp["diameter_1"]/2
            r_2 = self._inp["diameter_2"]/2
            l = self._inp["length"]
            return (r_2-r_1)/(l-1)

        x = r(z)*np.cos(phi)
        y = r(z)*np.sin(phi)
        z = 1

        return [x, y, z]


    ############
    # Features #
    ############
    def normal(self, pos):
        """Calculate unit normal vector on surface for a given position

        .. math::

            \\frac{\\partial\\Phi}{\\partial\\phi}(\\tilde r(z),\\phi,z)\\times
            \\frac{\\partial\\Phi}{\\partial z}(\\hat r(z),\\phi,z)=
            \\begin{bmatrix}\\tilde r(z)\\cos(\\phi)\\\\\\tilde r(z)\\sin(\\phi)\\\\\\tilde r(z)\\hat r(z)\\end{bmatrix}

        with polar angle :math:`\\phi` and radius functions along the :math:`z`-axis

        .. math::

            &\\tilde r(z)=r_1+\\frac{r_2-r_1}{l-1}(z-1)\\\\
            &\\hat r(z)=\\frac{r_2-r_1}{l-1},

        with radii :math:`r_1` and :math:`r_2` and cone length :math:`l`.

        Parameters
        ----------
        pos : list
            Position

        Returns
        -------
        normal : list
            Normal vector
        """
        # Initialize
        x, y, z = self.convert(pos)

        # Cartesian to polar
        r = math.sqrt(x**2+y**2)
        phi = geometry.angle_polar([x, y, z])

        # Calculate derivatives
        d_Phi_phi = self.d_Phi_phi(r, phi, z)
        d_Phi_z = self.d_Phi_z(r, phi, z)

        # Calculate normal vector
        return geometry.cross_product(d_Phi_phi, d_Phi_z)

    def is_in(self, pos):
        """Check if given position is inside of shape.

        Parameters
        ----------
        pos : list
            Position

        Returns
        -------
        is_in : bool
            True if position is inside of shape
        """
        def r(z):
            r_1 = self._inp["diameter_1"]/2
            r_2 = self._inp["diameter_2"]/2
            l = self._inp["length"]
            return r_1+(r_2-r_1)/(l-1)*(z-1)


        # Check if within shape
        pos_zero = self.convert(pos)
        length = geometry.length(geometry.cross_product(self._inp["central"], geometry.vector([0, 0, 0], pos_zero)))/geometry.length(self._inp["central"])

        if length < r(pos_zero[2]):
            return pos_zero[2]>0 and pos_zero[2]<self._inp["length"]
        else:
            return False


    #########
    # Shape #
    #########
    def rim(self, z, num=100):
        """Return x and y values for given z-position.

        Parameters
        ----------
        z : float
            Position on the axis
        num : integer, optional
            Number of points

        Returns
        -------
        positions : list
            x and y arrays of the surface rim on the z-position
        """
        phi = np.linspace(0, 2*np.pi, num)
        r = self._inp["diameter_1"]/2

        return self.Phi(r, phi, [z])

    def surf(self, num=100):
        """Return x, y and z values for the shape.

        Parameters
        ----------
        num : integer, optional
            Number of points

        Returns
        -------
        positions : list
            x, y and z arrays of the surface rim
        """
        phi = np.linspace(0, 2*np.pi, num)
        r = np.linspace(self._inp["diameter_1"]/2, self._inp["diameter_2"]/2, num)
        z = np.linspace(0, self._inp["length"], num)

        return self.Phi(r, phi, z)


    ##############
    # Properties #
    ##############
    def volume(self):
        """Calculate volume

        .. math::

            V=\\frac{1}{3}\\pi \\left[r_1^2+r_2^2+r_1r_2\\right]l

        with radii :math:`r_1` and :math:`r_2` and cone length :math:`l`.

        Returns
        -------
        volume : float
            Volume
        """
        r_1 = self._inp["diameter_1"]/2
        r_2 = self._inp["diameter_2"]/2
        l = self._inp["length"]
        return 1/3*math.pi*(r_1**2+r_2**2+r_1*r_2)*l

    def surface(self):
        """Calculate inner surface

        .. math::

            S=\\pi (r_1+r_2)l

        with radii :math:`r_1` and :math:`r_2` and cone length :math:`l`.

        Returns
        -------
        surface : float
            Inner surface
        """
        r_1 = self._inp["diameter_1"]/2
        r_2 = self._inp["diameter_2"]/2
        l = self._inp["length"]
        return math.pi*(r_1+r_2)*math.sqrt((r_1-r_2)**2+l**2)


class Hourglass(Shape):
    """Cylindrical shape whose radius varies as a cosine along the axis,
    narrowing to a minimum at the mid-point.

    Required inputs:

    * **central** - Central axis unit vector
    * **centroid** - Centroid of the block
    * **length** - Total length along the axis
    * **diameter_outer** - Diameter at both ends (widest)
    * **diameter_inner** - Diameter at the mid-point (narrowest)

    The radius profile is

    .. math::

        r(z) = r_i + (r_o - r_i)\\,\\frac{1 + \\cos\\!\\left(\\frac{2\\pi z}{l}\\right)}{2}

    with outer radius :math:`r_o`, inner radius :math:`r_i` and length
    :math:`l`.

    Parameters
    ----------
    inp : dictionary
        Dictionary of necessary inputs
    """
    def __init__(self, inp):
        self._centroid = [0, 0, inp["length"] / 2]
        super(Hourglass, self).__init__(inp)

    def _r(self, z):
        r_o = self._inp["diameter_outer"] / 2
        r_i = self._inp["diameter_inner"] / 2
        l = self._inp["length"]
        return r_i + (r_o - r_i) * (1 + np.cos(2 * np.pi * z / l)) / 2

    def _dr_dz(self, z):
        r_o = self._inp["diameter_outer"] / 2
        r_i = self._inp["diameter_inner"] / 2
        l = self._inp["length"]
        return -(r_o - r_i) * np.pi / l * np.sin(2 * np.pi * z / l)

    ############
    # Function #
    ############
    def Phi(self, r, phi, z):
        """Surface parametrisation.

        Parameters
        ----------
        r : ignored
            Kept for interface consistency; radius is determined by ``z``.
        phi : array-like
            Polar angles
        z : array-like
            Positions along the axis

        Returns
        -------
        pos : list
            x, y, z arrays of the surface
        """
        z_arr = np.asarray(z)
        r_arr = self._r(z_arr)
        phi_arr = np.asarray(phi)

        x = np.outer(r_arr, np.cos(phi_arr))
        y = np.outer(r_arr, np.sin(phi_arr))
        z_grid = np.outer(z_arr, np.ones(len(phi_arr)))

        return self.convert([x, y, z_grid], False)

    def d_Phi_phi(self, r, phi, z):
        """Azimuthal tangent vector.

        Parameters
        ----------
        r : float
            Radius at z
        phi : float
            Polar angle
        z : float
            Position along axis

        Returns
        -------
        tangent : list
        """
        return [-r * np.sin(phi), r * np.cos(phi), 0]

    def d_Phi_z(self, r, phi, z):
        """Axial tangent vector.

        Parameters
        ----------
        r : float
            Radius at z (unused; r'(z) is recomputed internally)
        phi : float
            Polar angle
        z : float
            Position along axis

        Returns
        -------
        tangent : list
        """
        dr = self._dr_dz(z)
        return [dr * np.cos(phi), dr * np.sin(phi), 1]

    ############
    # Features #
    ############
    def normal(self, pos):
        """Outward surface normal at ``pos``.

        Parameters
        ----------
        pos : list
            Global position on the surface

        Returns
        -------
        normal : list
            Normal vector (unnormalised)
        """
        x, y, z = self.convert(pos)
        r = math.sqrt(x ** 2 + y ** 2)
        phi = geometry.angle_polar([x, y, z])
        return geometry.cross_product(self.d_Phi_phi(r, phi, z), self.d_Phi_z(r, phi, z))

    def is_in(self, pos):
        """Return True if ``pos`` is inside the hourglass.

        Parameters
        ----------
        pos : list
            Global position

        Returns
        -------
        is_in : bool
        """
        pos_local = self.convert(pos)
        z = pos_local[2]
        if z <= 0 or z >= self._inp["length"]:
            return False
        r = math.sqrt(pos_local[0] ** 2 + pos_local[1] ** 2)
        return r < self._r(z)

    #########
    # Shape #
    #########
    def rim(self, z, num=100):
        """Circle of surface points at axial position ``z``.

        Parameters
        ----------
        z : float
            Position along the axis (local frame)
        num : int, optional
            Number of points

        Returns
        -------
        positions : list
            x, y, z arrays
        """
        phi = np.linspace(0, 2 * np.pi, num)
        return self.Phi(None, phi, [z])

    def surf(self, num=100):
        """Parametric surface mesh.

        Parameters
        ----------
        num : int, optional
            Number of points in each direction

        Returns
        -------
        positions : list
            x, y, z arrays
        """
        phi = np.linspace(0, 2 * np.pi, num)
        z = np.linspace(0, self._inp["length"], num)
        return self.Phi(None, phi, z)

    ##############
    # Properties #
    ##############
    def volume(self):
        """Volume enclosed by the hourglass.

        .. math::

            V = \\pi l \\left(r_i^2 + r_i\\,\\Delta r + \\frac{3\\,\\Delta r^2}{8}\\right)

        with :math:`\\Delta r = r_o - r_i`.

        Returns
        -------
        volume : float
        """
        r_o = self._inp["diameter_outer"] / 2
        r_i = self._inp["diameter_inner"] / 2
        l = self._inp["length"]
        dr = r_o - r_i
        return math.pi * l * (r_i ** 2 + r_i * dr + 3 * dr ** 2 / 8)

    def surface(self):
        """Lateral surface area (numerical integration).

        Returns
        -------
        surface : float
        """
        from scipy.integrate import quad
        l = self._inp["length"]

        def integrand(z):
            return 2 * math.pi * self._r(z) * math.sqrt(1 + self._dr_dz(z) ** 2)

        area, _ = quad(integrand, 0, l)
        return area
