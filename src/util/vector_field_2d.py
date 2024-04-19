import numpy as np
import cartopy.crs as ccrs


class VectorField2d:
    """
    2 dimensional vector field class. Constructor gives the mesh grid of x and y coordinate.
    """

    def __init__(self, _lx: float, _ly: float, _grid_width: float):
        """

        :param _lx: maximum x coordinate value, region is -_lx to _lx
        :param _ly: maximum y coordinate value, region is -_ly to _ly
        :param _grid_width: grid width
        """
        self.LX = _lx
        self.LY = _ly
        self.grid_width = _grid_width
        self.X, self.Y = np.meshgrid(np.arange(-self.LX, self.LX * 1.01, self.grid_width),
                                     np.arange(-self.LY, self.LY * 1.01, self.grid_width))

    def get_grid(self):
        """

        :return: array of mesh grid X and Y
        """
        return self.X, self.Y

    def draw_limits(self, _plt):
        """
        draw_limits function set the limits of x and y coordinate.
        :param _plt: plotter instance
        """
        _plt.set_aspect('equal')
        gl = _plt.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        _plt.set_extent([-self.LX * 1.1, self.LX * 1.1, -self.LY * 1.1, self.LY * 1.1], ccrs.PlateCarree())

    def draw(self, _u, _v, _plt, _title_string, _arrow_scale):
        _max = max(np.max(_u), -np.min(_u), np.max(_v), -np.min(_v))
        print(str(_max * 3600), end=',')
        _scale = _max / self.grid_width
        _arrow = _arrow_scale / _scale / 3600
        if _arrow_scale < 1e-3:
            arrow_string = "{:.1f}".format(_arrow_scale * 1e6) + "uas"
        elif _arrow_scale < 1:
            arrow_string = "{:.1f}".format(_arrow_scale * 1e3) + "mas"
        else:
            arrow_string = "{:.1f}".format(_arrow_scale) + "as"
        self.draw_limits(_plt)
        _plt.arrow(x=0, y=-self.LY, dx=_arrow, dy=0, head_width=0.25 * _arrow, length_includes_head=True)
        _plt.text(_arrow, -self.LY, arrow_string)
        _plt.quiver(self.X, self.Y, _u, _v, color='red', scale_units='xy', scale=_scale, transform=ccrs.PlateCarree())
        _plt.set_title(_title_string)
