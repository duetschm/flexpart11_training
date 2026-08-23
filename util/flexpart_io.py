import numpy as np
import xarray as xr
import glob
import f90nml
import datetime


def create_datelist(startdate, stopdate, step=datetime.timedelta(hours=1)):
    """
    Constructs a list of dates.

    Parameters
    ----------
    startdate : datetime.datetime
        First date
    stopdate : datetime.datetime
        Last date
    step : datetime.timedelta
        Difference between dates

    Returns
    -------
    datelist : list
        List of dates between startdate and stopdate
    """

    datelist = []
    date = startdate
    while date <= stopdate:
        datelist.append(date)
        date += step

    return datelist


def read_namelists(flxdir):
    """
    Read FLEXPART namelists.

    Parameters
    ----------
    flxdir : str
        Directory containing the FLEXPART namelists.

    Returns
    -------
    namelist
        f90nml namelist object with all FLEXPART namelists.
    """
    nml = f90nml.read(f"{flxdir}/RELEASES.namelist")
    for fn in glob.glob(f"{flxdir}/*.namelist"):
        if not "RELEASES" in fn:
            nml.update(f90nml.read(fn))

    return nml


def decode_binary(filenames, ntim, nrel, nage, nlev, nlat, nlon):
    """
    Decode FLEXPART binary sensitivity files.

    Parameters
    ----------
    filenames : sequence of str
        Paths to the FLEXPART binary output files.
    ntim : int
        Number of time steps to decode.
    nrel : int
        Number of release points.
    nage : int
        Number of age classes.
    nlev : int
        Number of vertical levels.
    nlat : int
        Number of latitude grid points.
    nlon : int
        Number of longitude grid points.

    Returns
    -------
    numpy.ndarray
        Decoded sensitivity array with shape
        (ntim, nage, nrel, nlev, nlat, nlon).
    """
    len_dummy_entries = 20
    sparse_dump_i = []
    sparse_dump_r = []
    for filename in filenames:
        offset = 3 * 4
        while True:
            try:
                sp_count_i = np.fromfile(
                    filename,
                    dtype="<i4",
                    count=1,
                    offset=offset + 4 * (len_dummy_entries + 1),
                )[0]
                sparse_dump_i.append(
                    np.fromfile(
                        filename,
                        dtype="<i4",
                        count=sp_count_i,
                        offset=offset + 4 * (len_dummy_entries + 4),
                    )
                )
                sp_count_r = np.fromfile(
                    filename,
                    dtype="<i4",
                    count=1,
                    offset=offset + 4 * (len_dummy_entries + 4 + sp_count_i + 2),
                )[0]
                sparse_dump_r.append(
                    np.fromfile(
                        filename,
                        dtype="<f4",
                        count=sp_count_r,
                        offset=offset + 4 * (len_dummy_entries + 4 + sp_count_i + 5),
                    )
                )
                offset += 4 * (len_dummy_entries + 4 + sp_count_i + 5 + sp_count_r + 1)
            except IndexError:  # end of file
                break

    sens_1d = np.zeros((ntim, nrel * nage, nlon * nlat * nlev))

    for t in range(ntim):
        for r in range(nrel * nage):
            rind = 0
            sind = nrel * nage * t + r
            for i in range(len(sparse_dump_i[sind])):
                sign = np.sign(sparse_dump_r[sind][rind])
                n = 0
                while (
                    rind + n < len(sparse_dump_r[sind]) - 1
                    and sparse_dump_r[sind][rind + n] * sign > 0
                ):
                    sens_1d[t, r, sparse_dump_i[sind][i] - (nlon * nlat) + n] = (
                        sign * sparse_dump_r[sind][rind + n]
                    )
                    n += 1
                rind += n

    sens = sens_1d.reshape((ntim, nage, nrel, nlev, nlat, nlon))

    return sens


def output_units(ldirect, ind_source, ind_receptor):
    """
    Return the units for FLEXPART output sensitivities.
    Taken the subroutine output_units in from netcdf_output_mod.90.

    Parameters
    ----------
    ldirect : int
        Simulation direction: 1 for forward simulations and -1
        for backward simulations.
    ind_source : int
        Source indicator. 1 for mass and 2 for mass mixing ratio.
    ind_receptor : int
        Receptor indicator. 1 for mass, 2 for mass mixing ratio,
        3 for wet deposition, and 4 for dry deposition.

    Returns
    -------
    str
        Units corresponding to the specified simulation, source, and
        receptor indicators.
    """
    if ldirect == 1:  # forward simulation
        if ind_source == 1:
            if ind_receptor == 1:
                units = "ng m-3"
            else:
                units = "ng kg-1"
        else:
            if ind_receptor == 1:
                units = "ng m-3"
            else:
                units = "ng kg-1"
    else:  # backward simulation
        if ind_source == 1:
            if ind_receptor == 1:
                units = "s"
            else:
                units = "s m3 kg-1"
        else:
            if ind_receptor == 1:
                units = "s kg m-3"
            else:
                units = "s"

    return units


def read_binary(flxdir, filepattern):
    """
    Load binary sensitivities into an xarray.Dataset.

    Parameters
    ----------
    flxdir : str
        Directory containing the FLEXPART namelists and binary output files.
    filepattern : str
        Prefix used to select the binary output files.

    Returns
    -------
    xarray.Dataset
        Dataset containing the decoded sensitivity field on the output grid.
    """

    nml = read_namelists(flxdir)

    if type(nml["release"]) == f90nml.namelist.Cogroup:
        nrel = len(nml["release"])
    else:
        nrel = 1

    ldirect = nml["command"]["ldirect"]
    bdate = datetime.datetime.strptime(
        f"{nml['command']['ibdate']:08d}{nml['command']['ibtime']:06d}", "%Y%m%d%H%M%S"
    )
    edate = datetime.datetime.strptime(
        f"{nml['command']['iedate']:08d}{nml['command']['ietime']:06d}", "%Y%m%d%H%M%S"
    )
    loutstep = nml["command"]["loutstep"]
    ntim = int((edate - bdate).total_seconds() / loutstep)
    times = np.array(
        create_datelist(bdate, edate, step=datetime.timedelta(seconds=loutstep))[
            ::ldirect
        ]
    )[1:]

    nage = nml["nage"]["nageclass"]

    nlev = len(nml["outgrid"]["outheights"])
    nlon = nml["outgrid"]["numxgrid"]
    nlat = nml["outgrid"]["numygrid"]

    lon0 = nml["outgrid"]["outlon0"]
    dx = nml["outgrid"]["dxout"]
    lat0 = nml["outgrid"]["outlat0"]
    dy = nml["outgrid"]["dyout"]

    lon = np.arange(lon0, lon0 + nlon * dx, dx)
    lat = np.arange(lat0, lat0 + nlat * dy, dy)

    filenames = sorted(glob.glob(f"{flxdir}/{filepattern}*[!nc]"))[::ldirect]
    sens = decode_binary(filenames, ntim, nrel, nage, nlev, nlat, nlon)
    units = output_units(
        nml["command"]["ldirect"],
        nml["command"]["ind_source"],
        nml["command"]["ind_receptor"],
    )

    ds = xr.Dataset(
        data_vars={
            "spec001_mr": (
                ["time", "pointspec", "nageclass", "height", "latitude", "longitude"],
                sens,
                {
                    "long_name": nml["species_params"]["pspecies"].rstrip(),
                    "units": units,
                },
            ),
            "RELLNG1": (
                ["pointspec"],
                np.array([nml["release"][i]["lon1"] for i in range(nrel)])
                if nrel > 1
                else nml["release"]["lon1"]
            ),
            "RELLNG2": (
                ["pointspec"],
                np.array([nml["release"][i]["lon2"] for i in range(nrel)])
                if nrel > 1
                else nml["release"]["lon2"]
            ),
            "RELLAT1": (
                ["pointspec"],
                np.array([nml["release"][i]["lat1"] for i in range(nrel)])
                if nrel > 1
                else nml["release"]["lat1"]
            ),
            "RELLAT2": (
                ["pointspec"],
                np.array([nml["release"][i]["lat2"] for i in range(nrel)])
                if nrel > 1
                else nml["release"]["lat2"]
            ),
        },
        coords={
            "longitude": ("longitude", lon),
            "latitude": ("latitude", lat),
            "height": ("height", nml["outgrid"]["outheights"]),
            "pointspec": ("pointspec", np.arange(nrel)),
            "nageclass": ("nageclass", np.arange(nage)),
            "time": ("time", times),
        },
    )

    return ds


def open_flexpart_grid(flxdir, fileformat="netcdf", filepattern="grid"):
    """
    Open and process a FLEXPART output file.

    Parameters
    ----------
    flxdir : str
        Directory containing the FLEXPART output files.
    fileformat : {"netcdf", "binary"}, default="netcdf"
        Format of the FLEXPART output to open.
    filepattern : str, default="grid"
        Prefix used to select the FLEXPART output files.

    Returns
    -------
    ds : xarray.Dataset
        Loaded FLEXPART dataset with reduced dimensions.
    heights_exp : numpy.ndarray
        Heights array with 0 prepended.
    """
    # load dataset
    if fileformat == "netcdf":
        ds = xr.open_dataset(
            glob.glob(f"{flxdir}/{filepattern}*.nc")[0], decode_timedelta=True
        )
    elif fileformat == "binary":
        ds = read_binary(flxdir, filepattern)

    # assign new height values
    heights_exp = np.concatenate([[0.0], ds["height"]])

    return ds, heights_exp


def open_topography(filename):
    """
    Open and process a topography dataset.

    Parameters
    ----------
    filename : str
        Path to the topography NetCDF file.

    Returns
    -------
    xarray.Dataset
        Topography dataset with masked negative values and updated coordinates.
    """
    etopo = xr.open_dataset(filename)
    etopo["topo"].values = np.ma.MaskedArray(
        etopo["topo"], mask=etopo["topo"] < 0.0
    ).filled(0.0)
    etopo = etopo.assign_coords(
        lon=("lon", etopo["topo_lon"].values), lat=("lat", etopo["topo_lat"].values)
    )

    return etopo


def interpolate_values(ds, etopo, lon1, lat1, lon2, lat2):
    """
    Interpolate FLEXPART and topography data along a line.

    Parameters
    ----------
    ds : xarray.Dataset
        FLEXPART dataset to interpolate.
    etopo : xarray.Dataset
        Topography dataset to interpolate.
    lon1 : float
        Starting longitude of the line.
    lat1 : float
        Starting latitude of the line.
    lon2 : float
        Ending longitude of the line.
    lat2 : float
        Ending latitude of the line.

    Returns
    -------
    ds_int : xarray.Dataset
        Interpolated FLEXPART data along the line.
    distance : numpy.ndarray
        Cumulative distance along the line in km.
    etopo_int_avg : xarray.Dataset
        Topography interpolated at midpoints between line segments.
    etopo_int : xarray.Dataset
        Topography interpolated at all line points.
    """

    # build points along the line
    npts = 201
    lons = np.linspace(lon1, lon2, npts)
    lats = np.linspace(lat1, lat2, npts)
    lons_avg = (lons[1:] + lons[:-1]) / 2.0
    lats_avg = (lats[1:] + lats[:-1]) / 2.0

    # cumulative distance (km)
    g = Geod(ellps="WGS84")
    d = [0.0]
    for i in range(1, npts):
        _, _, dist = g.inv(lons[i - 1], lats[i - 1], lons[i], lats[i])
        d.append(d[-1] + dist / 1000.0)
    distance = np.array(d)
    distance_avg = (distance[1:] + distance[:-1]) / 2.0

    # flexpart
    ds_int = ds.interp(longitude=("points", lons_avg), latitude=("points", lats_avg))
    ds_int = ds_int.assign_coords(distance=("points", distance_avg))
    ds_int["distance"].attrs.update(units="km", long_name="cumulative distance")

    # topography
    etopo_int_avg = etopo.interp(lon=("points", lons_avg), lat=("points", lats_avg))
    etopo_int = etopo.interp(lon=("points", lons), lat=("points", lats))

    return ds_int, distance, etopo_int_avg, etopo_int