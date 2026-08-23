import numpy as np
#import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import cartopy.feature as cfeature
import cartopy.crs as ccrs
#import datetime
from pyproj import Geod
from cartopy.mpl.contour import GeoContourSet
from matplotlib.animation import FuncAnimation

plt.rcParams["animation.embed_limit"] = 50


def draw_map(ax, lon, lat, rel_coords):
    """
    Draw a map with coastlines, borders, lakes, and grid on a Cartopy axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Cartopy axes object to draw on.
    lon : array_like
        Longitude coordinates of the domain.
    lat : array_like
        Latitude coordinates of the domain.
    rel_coords : list of lists
        Coordinates of the release points / boxes
        (lon1, lat1, lon2, lat2).

    Returns
    -------
    matplotlib.axes.Axes
        The modified axes object with map features.
    """

    ax.add_feature(cfeature.OCEAN.with_scale("110m"), facecolor=[1.0] * 3)
    ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor=[0.9] * 3)
    ax.add_feature(
        cfeature.BORDERS.with_scale("110m"),
        linewidth=0.2,
        edgecolor=[0.4] * 3,
    )
    lakes = cfeature.NaturalEarthFeature(
        "physical",
        "lakes",
        "110m",
        edgecolor="black",
        facecolor="white",
    )
    ax.add_feature(lakes, zorder=2)

    # plot release point/box
    for rel_coord in rel_coords:
        if rel_coord[2] == rel_coord[0]:
            plt.scatter(
                rel_coord[0],
                rel_coord[1],
                color="black",
                linewidth=1.5,
                marker="x",
                transform=ccrs.PlateCarree(),
                zorder=3,
            )
        else:
            plt.plot(
                [rel_coord[0], rel_coord[2], rel_coord[2], rel_coord[0], rel_coord[0]],
                [rel_coord[1], rel_coord[1], rel_coord[3], rel_coord[3], rel_coord[1]],
                color="black",
                linewidth=1.5,
                transform=ccrs.PlateCarree(),
                zorder=3,
            )

    # plot outgrid
    lons, lats = np.meshgrid(lon, lat)
    plt.plot(
        np.hstack([lons[:, 0], lons[-1, :], lons[::-1, -1], lons[0, ::-1]]),
        np.hstack([lats[:, 0], lats[-1, :], lats[::-1, -1], lats[0, ::-1]]),
        color="black",
        linewidth=0.5,
        transform=ccrs.PlateCarree(),
    )

    # plot gridlines with labels
    gl = ax.gridlines(
        draw_labels=True,
        color="lightgrey",
        linestyle="dashed",
        linewidth=0.5,
        xlocs=np.arange(-180, 180, 20),
        ylocs=np.arange(-80, 90, 10),
    )

    ax.coastlines(resolution="110m")

    return ax


def setup_fig(ds):
    """
    Set up a figure with orthographic projection and map.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing longitude, latitude, and release point information.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    ax : matplotlib.axes.Axes
        Cartopy axes with orthographic projection and map features.
    gs : matplotlib.gridspec.GridSpec
        GridSpec object for layout management.
    """
    gs = gridspec.GridSpec(1, 3, width_ratios=[40, 1, 2])

    central_lon = float(np.mean(ds["longitude"]))
    central_lat = float(np.mean(ds["latitude"]))

    proj = ccrs.Orthographic(
        central_longitude=central_lon,
        central_latitude=central_lat,
    )

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(gs[0, 0], projection=proj)

    rel_coords = set(
        zip(
            ds["RELLNG1"].values,
            ds["RELLAT1"].values,
            ds["RELLNG2"].values,
            ds["RELLAT2"].values,
        )
    )
    ax = draw_map(
        ax,
        ds["longitude"],
        ds["latitude"],
        rel_coords,
    )

    return fig, ax, gs


def plot_map(ds, cmap, levels, norm, lon1=None, lat1=None, lon2=None, lat2=None):
    """
    Plot a map with FLEXPART concentration data.

    Parameters
    ----------
    ds : xarray.Dataset
        FLEXPART dataset at one time step.
    cmap : str or matplotlib.colors.Colormap
        Colormap to use for plotting.
    levels : array_like
        Contour levels for the plot.
    norm : matplotlib.colors.Normalize
        Normalization for colormap scaling.
    lon1 : float, optional
        Starting longitude for cross-section line.
    lat1 : float, optional
        Starting latitude for cross-section line.
    lon2 : float, optional
        Ending longitude for cross-section line.
    lat2 : float, optional
        Ending latitude for cross-section line.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object containing the map plot.
    """
    # plot map
    fig, ax, gs = setup_fig(ds)

    # mask data
    data = ds["spec001_mr"].isel(time=0)
    data.values = np.ma.MaskedArray(data, mask=data == 0.0).filled(1e-10)

    # add title
    ax.set_title(ds["spec001_mr"].long_name, loc="left")
    ax.set_title(np.datetime_as_string(data["time"], unit="h"), loc="right")

    # plot data
    for r in data["pointspec"]:
        for a in data["nageclass"]:
            h1 = ax.contourf(
                data["longitude"],
                data["latitude"],
                data.sel(pointspec=r, nageclass=a),
                norm=norm,
                cmap=cmap,
                levels=levels,
                extend="max",
                transform=ccrs.PlateCarree(),
                zorder=2,
            )

    if all(l is not None for l in [lon1, lat1, lon2, lat2]):
        plt.plot(
            [lon1, lon2], [lat1, lat2], transform=ccrs.PlateCarree(), color="black"
        )
        plt.scatter(
            [lon1],
            [lat1],
            marker="^",
            transform=ccrs.PlateCarree(),
            color="black",
            zorder=3,
        )
        plt.scatter(
            [lon2],
            [lat2],
            marker="v",
            transform=ccrs.PlateCarree(),
            color="black",
            zorder=3,
        )

    cax = fig.add_subplot(gs[0, 2])
    cbar = plt.colorbar(h1, ax=ax, cax=cax, label=ds["spec001_mr"].units)

    nearest_int = np.rint(np.log10(levels))
    is_whole = np.isclose(np.log10(levels), nearest_int, atol=1e-12)
    cbar.set_ticks(levels[is_whole])

    return fig


def plot_map_anim(ds, cmap, levels, norm):
    """
    Create an animated map showing FLEXPART concentration evolution over time.

    Parameters
    ----------
    ds : xarray.Dataset
        FLEXPART dataset with time dimension.
    cmap : str or matplotlib.colors.Colormap
        Colormap to use for plotting.
    levels : array_like
        Contour levels for the plot.
    norm : matplotlib.colors.Normalize
        Normalization for colormap scaling.

    Returns
    -------
    matplotlib.animation.FuncAnimation
        Animation object showing concentration data over time.
    """
    # plot map
    fig, ax, gs = setup_fig(ds)

    # mask data
    data = ds["spec001_mr"].weighted(ds["height"]).mean(dim="height").sortby("time")
    data.values = np.ma.MaskedArray(data, mask=data == 0.0).filled(1e-10)

    # add title
    ax.set_title(ds["spec001_mr"].long_name, loc="left")

    # plot data
    def update(n, data, cmap, levels, norm):
        nall = data["time"].size - 1
        print(f"\rPlotting frame {n}/{nall}", end="", flush=True)

        h1_list = ax.findobj(lambda x: isinstance(x, GeoContourSet))
        if len(h1_list) > 0:
            for h1 in h1_list:
                h1.remove()
        for r in data["pointspec"]:
            for a in data["nageclass"]:
                h1 = ax.contourf(
                    data["longitude"],
                    data["latitude"],
                    data.isel(time=n,pointspec=r,nageclass=a),
                    norm=norm,
                    cmap=cmap,
                    levels=levels,
                    extend="max",
                    transform=ccrs.PlateCarree(),
                    zorder=2,
                )

        ax.set_title(
            np.datetime_as_string(data["time"].isel(time=n), unit="h"), loc="right"
        )

        return h1

    # initialize h1 for colorbar
    h1 = update(0, data, cmap, levels, norm)
    cax = fig.add_subplot(gs[0, 2])
    cbar = plt.colorbar(h1, ax=ax, cax=cax, label=ds["spec001_mr"].units)

    nearest_int = np.rint(np.log10(levels))
    is_whole = np.isclose(np.log10(levels), nearest_int, atol=1e-12)
    cbar.set_ticks(levels[is_whole])

    # plot animation
    ani = FuncAnimation(
        fig,
        update,
        fargs=(data, cmap, levels, norm),
        frames=data["time"].size,
        blit=False,
        interval=50,
    )
    plt.close(fig)

    return ani


def plot_cross_section(
    ds_int, distance, heights_exp, cmap, norm, etopo_int, etopo_int_avg
):
    """
    Plot a cross-section of FLEXPART concentration with topography.

    Parameters
    ----------
    ds_int : xarray.Dataset
        Interpolated FLEXPART data along the cross-section line.
    distance : numpy.ndarray
        Cumulative distance along the cross-section line in km.
    heights_exp : numpy.ndarray
        Height levels for the cross-section.
    cmap : str or matplotlib.colors.Colormap
        Colormap to use for plotting.
    norm : matplotlib.colors.Normalize
        Normalization for colormap scaling.
    etopo_int : xarray.Dataset
        Topography interpolated at all cross-section points.
    etopo_int_avg : xarray.Dataset
        Topography interpolated at midpoints between cross-section points.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object containing the cross-section plot.
    """
    # plot cross section

    distance2d = np.repeat(distance[None, :], heights_exp.size, axis=0)
    heightstopo = etopo_int["topo"].values[None, :] + heights_exp[:, None]

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    h1 = ax.pcolor(
        distance2d,
        heightstopo,
        ds_int["spec001_mr"].isel(time=0),
        cmap=cmap,
        shading="flat",
        norm=norm,
    )
    ax.plot(ds_int["distance"], etopo_int_avg["topo"], color="black")
    ax.fill_between(ds_int["distance"], etopo_int_avg["topo"], color=[0.9] * 3)
    ax.plot(distance2d.T, heightstopo.T, color="black", linewidth=0.5, alpha=0.1)
    ax.scatter(
        [0],
        [-0.07],
        marker="^",
        linestyle="None",
        transform=ax.transAxes,
        clip_on=False,
        color="black",
    )
    ax.scatter(
        [1],
        [-0.07],
        marker="v",
        linestyle="None",
        transform=ax.transAxes,
        clip_on=False,
        color="black",
    )

    plt.ylim((0.0, 7000))
    plt.xlabel("distance (km)")
    plt.ylabel("height (m)")
    plt.colorbar(h1, label=ds_int["spec001_mr"].units)

    return fig