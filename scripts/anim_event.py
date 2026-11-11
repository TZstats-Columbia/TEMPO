# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ce Bian
"""Animate one event of the North American heatwave-object catalogue over hourly ERA5 2 m temperature.

Reads the released dataset (Bian and Zheng 2026, doi:10.5281/zenodo.23966834): heatwave_events.csv and heatwave_event_labels.nc.  The shading is
hourly ERA5 2 m temperature, either from a local netCDF file or read anonymously from the public ARCO-ERA5 store on Google Cloud.

Map: 2 m temperature every --step hours; the chosen event's footprint on that UTC day (white outline, hatched inside); the footprints of
all other catalogued events on that day (translucent white); a marked point (e.g. a city).  Lower panel: hourly 2 m temperature at the
point, with the days on which the point lies inside the event's footprint shaded and a cursor moving with the map.

Examples
  # the July 2013 New York event, picked by the label at New York on 18 July, ERA5 from ARCO-ERA5
  python scripts/anim_event.py --data na025_heatwave_objects_v1 --point -74.01 40.71 --date 2013-07-18 --name NYC \\
      --t2m arco --start 2013-07-10 --end 2013-07-31 --box -97 -60 30 54 --out anim_nyc_2013.mp4
  # the same event by identifier, with a local ERA5 file (variable t2m in K, dims time/valid_time x latitude x longitude)
  python scripts/anim_event.py --data na025_heatwave_objects_v1 --event G34542 --point -74.01 40.71 --name NYC --t2m era5_t2m.nc --out anim.mp4

Outputs: the MP4 (--out), a poster frame (<out>_poster.png) and a strip of six frames (<out>_frames.png).
Requirements: numpy, pandas, xarray, netCDF4, matplotlib, cartopy, imageio-ffmpeg (or ffmpeg on the PATH); for --t2m arco also zarr, gcsfs.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.animation import FFMpegWriter, FuncAnimation

ARCO = "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
BG, FG, ACCENT = "black", "white", "#ff5a1f"


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="folder of the released dataset (heatwave_events.csv, heatwave_event_labels.nc)")
    ap.add_argument("--event", help="event identifier, e.g. G34542 (default: the event at --point on --date)")
    ap.add_argument("--point", type=float, nargs=2, metavar=("LON", "LAT"), required=True, help="marked location, e.g. -74.01 40.71")
    ap.add_argument("--date", help="UTC day used with --point to pick the event, e.g. 2013-07-18")
    ap.add_argument("--name", default="", help="label of the marked location, e.g. NYC")
    ap.add_argument("--t2m", default="arco", help="'arco' (read ERA5 from ARCO-ERA5) or a local netCDF of hourly 2 m temperature in K")
    ap.add_argument("--start", help="first UTC day shown (default: 3 days before onset)")
    ap.add_argument("--end", help="last UTC day shown (default: 9 days after the event's end)")
    ap.add_argument("--box", type=float, nargs=4, metavar=("LON0", "LON1", "LAT0", "LAT1"), help="map extent (default: the local file's, or the event's +-12 deg)")
    ap.add_argument("--utc-offset", type=int, default=-4, help="hours added to UTC for the clock (default -4)")
    ap.add_argument("--tz", default="EDT", help="time-zone label of the clock")
    ap.add_argument("--step", type=int, default=2, help="hours per frame")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--vmin", type=float, default=12); ap.add_argument("--vmax", type=float, default=40)
    ap.add_argument("--no-others", action="store_true", help="do not draw the other events of the day")
    ap.add_argument("--dpi", type=int, default=120)
    ap.add_argument("--out", default="anim_event.mp4")
    return ap.parse_args()


def pick_event(data, events, point, date):
    """Event identifier from --event, or the label found at --point on --date."""
    lab = xr.open_dataset(data / "heatwave_event_labels.nc").event_label
    v = int(lab.sel(time=date, latitude=point[1], longitude=point[0], method="nearest"))
    if v == 0:
        raise SystemExit(f"no catalogued event at {point} on {date}")
    return events.loc[events.label == v].iloc[0]


def read_t2m(src, box, start, end):
    """Hourly 2 m temperature (deg C) as (time, latitude, longitude), latitude north to south, longitude -180..180."""
    t0, t1 = pd.Timestamp(start), pd.Timestamp(end) + pd.Timedelta(hours=23)
    if src == "arco":
        ds = xr.open_zarr(ARCO, storage_options={"token": "anon"}, chunks=None)
        lo0, lo1 = box[0] % 360, box[1] % 360
        da = ds["2m_temperature"].sel(time=slice(t0, t1), latitude=slice(box[3], box[2]), longitude=slice(lo0, lo1)).load()
    else:
        ds = xr.open_dataset(src)
        tn = "valid_time" if "valid_time" in ds.dims else "time"
        da = ds["t2m"].rename({tn: "time"}).sel(time=slice(t0, t1))
    da = da.assign_coords(longitude=((da.longitude + 180) % 360) - 180).sortby("longitude").sortby("latitude", ascending=False)
    if box is not None:
        da = da.sel(latitude=slice(box[3], box[2]), longitude=slice(box[0], box[1]))
    return (da - 273.15).transpose("time", "latitude", "longitude")


def staircase(lat, lon, mask, up=8):
    """Upsampled coordinates and mask, so that contour(0.5) follows the edges of the 0.25-degree cells."""
    dlat, dlon = lat[1] - lat[0], lon[1] - lon[0]
    off = (np.arange(up) + .5) / up - .5
    la = (lat[:, None] + dlat * off[None]).ravel(); lo = (lon[:, None] + dlon * off[None]).ravel()
    return lo, la, np.kron(mask.astype(float), np.ones((up, up)))


def main():
    a = parse_args(); data = Path(a.data)
    events = pd.read_csv(data / "heatwave_events.csv", parse_dates=["onset_date", "end_date"])
    ev = events.loc[events.event_id == a.event].iloc[0] if a.event else pick_event(data, events, a.point, a.date)
    start = pd.Timestamp(a.start) if a.start else ev.onset_date - pd.Timedelta(days=3)
    end = pd.Timestamp(a.end) if a.end else ev.end_date + pd.Timedelta(days=9)
    box = a.box or (None if a.t2m != "arco" else [ev.peak_lon - 12, ev.peak_lon + 12, ev.peak_lat - 9, ev.peak_lat + 9])
    print(f"event {ev.event_id}: {ev.onset_date.date()} .. {ev.end_date.date()} ({ev.duration_days} days); window {start.date()} .. {end.date()}", flush=True)

    T = read_t2m(a.t2m, box, start, end)
    times = pd.DatetimeIndex(T.time.values); lat, lon = T.latitude.values, T.longitude.values
    series = T.sel(latitude=a.point[1], longitude=a.point[0], method="nearest").values
    lab = xr.open_dataset(data / "heatwave_event_labels.nc").event_label.sel(
        time=slice(start, end), latitude=slice(lat.max() + .2, lat.min() - .2), longitude=slice(lon.min() - .2, lon.max() + .2)).load()
    glat, glon = lab.latitude.values.astype(float), lab.longitude.values.astype(float)
    days = pd.DatetimeIndex(lab.time.values)
    fp = {d: lab.values[k] == ev.label for k, d in enumerate(days) if (lab.values[k] == ev.label).any()}
    other = {d: (lab.values[k] > 0) & (lab.values[k] != ev.label) for k, d in enumerate(days)}
    iy, ix = int(np.abs(glat - a.point[1]).argmin()), int(np.abs(glon - a.point[0]).argmin())
    inside = sorted(d for d, m in fp.items() if m[iy, ix])
    event_days = list(pd.date_range(ev.onset_date, ev.end_date, freq="D"))
    frames = np.arange(0, len(times), a.step); local = times + pd.Timedelta(hours=a.utc_offset)

    try:
        import imageio_ffmpeg; mpl.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        if shutil.which("ffmpeg") is None: raise SystemExit("needs ffmpeg: pip install imageio-ffmpeg")
    mpl.rcParams["hatch.color"] = FG; mpl.rcParams["hatch.linewidth"] = .6

    fig = plt.figure(figsize=(12, 8.2), facecolor=BG)
    mid_lon, mid_lat = float(lon.mean()), float(lat.mean())
    proj = ccrs.LambertConformal(central_longitude=mid_lon, central_latitude=mid_lat, standard_parallels=(mid_lat - 9, mid_lat + 3))
    ax = fig.add_axes([.02, .27, .86, .71], projection=proj); ax.set_facecolor(BG)
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], ccrs.PlateCarree())
    mesh = ax.pcolormesh(lon, lat, T.values[0], cmap="inferno", vmin=a.vmin, vmax=a.vmax, shading="nearest", transform=ccrs.PlateCarree(), rasterized=True)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), edgecolor="#d0d0d8", lw=.7)
    ax.add_feature(cfeature.STATES.with_scale("50m"), edgecolor="#8a8a92", lw=.35)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), edgecolor="#a0a0a8", lw=.5); ax.spines["geo"].set_edgecolor("#444")
    ax.plot(*a.point, marker="*", ms=16, mfc=FG, mec=BG, mew=1, transform=ccrs.PlateCarree(), zorder=6)
    if a.name:
        ax.text(a.point[0] + .6, a.point[1] - .9, a.name, color=FG, fontsize=13, weight="bold", transform=ccrs.PlateCarree(), zorder=6)
    clock = ax.text(.015, .975, "", transform=ax.transAxes, color=FG, fontsize=16, va="top", zorder=7, bbox=dict(fc=BG, ec="none", alpha=.6, pad=3))
    cax = fig.add_axes([.895, .33, .014, .6]); cb = fig.colorbar(mesh, cax=cax, extend="both")
    cb.set_label("2 m temperature (°C)", color=FG, fontsize=12); cb.ax.tick_params(colors=FG, labelsize=11); cb.outline.set_visible(False)

    bx = fig.add_axes([.07, .06, .81, .16], facecolor=BG)
    for d in inside:
        bx.axvspan(d + pd.Timedelta(hours=a.utc_offset), d + pd.Timedelta(hours=24 + a.utc_offset), color=ACCENT, alpha=.25, lw=0)
    bx.plot(local, series, color=FG, lw=1.3); cur = bx.axvline(local[0], color=ACCENT, lw=2); dot, = bx.plot([], [], "o", color=ACCENT, ms=7)
    bx.set_xlim(local[0], local[-1]); bx.set_ylabel(f"{a.name or 'point'} (°C)", color=FG, fontsize=12); bx.tick_params(colors=FG, labelsize=10)
    bx.xaxis.set_major_formatter(mpl.dates.DateFormatter("%b %d")); bx.xaxis.set_major_locator(mpl.dates.DayLocator(interval=3))
    for s in bx.spines.values(): s.set_color("#666")
    bx.grid(alpha=.2, color=FG)
    note = f"shaded: {a.name or 'point'} inside the footprint of {ev.event_id}" + ("" if a.no_others else "   |   map: hatched = this event, translucent = other events")
    bx.text(0, 1.03, note, transform=bx.transAxes, color="#ffb08a", fontsize=10, va="bottom")

    art = []
    def draw(k):
        i = frames[k]; day = times[i].normalize(); mesh.set_array(T.values[i].ravel())
        for x in art: x.remove()
        art.clear()
        if not a.no_others and day in other and other[day].any():
            st = staircase(glat, glon, other[day])
            art.append(ax.contourf(*st, levels=[.5, 1.5], colors=FG, alpha=.22, transform=ccrs.PlateCarree(), zorder=4))
            art.append(ax.contour(*st, levels=[.5], colors="#c8c8d0", linewidths=.5, transform=ccrs.PlateCarree(), zorder=4))
        if day in fp:
            st = staircase(glat, glon, fp[day])
            art.append(ax.contourf(*st, levels=[.5, 1.5], colors="none", hatches=["///"], transform=ccrs.PlateCarree(), zorder=5))
            art.append(ax.contour(*st, levels=[.5], colors=FG, linewidths=1.8, transform=ccrs.PlateCarree(), zorder=5))
        label = local[i].strftime(f"%Y-%m-%d  %H:00 {a.tz}")
        if day in event_days: label += f"   {ev.event_id} day {event_days.index(day) + 1}/{ev.duration_days}"
        clock.set_text(label); cur.set_xdata([local[i]]); dot.set_data([local[i]], [series[i]])

    out = Path(a.out); stem = out.with_suffix("")
    FuncAnimation(fig, draw, frames=len(frames)).save(out, writer=FFMpegWriter(fps=a.fps, bitrate=4000, extra_args=["-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p"]),
                                                      savefig_kwargs=dict(facecolor=BG), dpi=a.dpi)
    print("->", out)
    peak = pd.Timestamp(ev.peak_date) + pd.Timedelta(hours=20)                                    # poster: 20 UTC on the event's peak day
    draw(int(np.clip(np.searchsorted(times, peak) // a.step, 0, len(frames) - 1)))
    fig.savefig(f"{stem}_poster.png", facecolor=BG, dpi=a.dpi); print("->", f"{stem}_poster.png")
    strip = plt.figure(figsize=(18, 2.2), facecolor=BG)
    for k, j in enumerate(np.linspace(0, len(frames) - 1, 6).astype(int)):
        draw(j); fig.canvas.draw(); s = strip.add_axes([k / 6, 0, 1 / 6, 1]); s.imshow(np.asarray(fig.canvas.buffer_rgba())); s.axis("off")
    strip.savefig(f"{stem}_frames.png", facecolor=BG, dpi=a.dpi); print("->", f"{stem}_frames.png")


if __name__ == "__main__":
    main()
