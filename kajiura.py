import numpy as np
import netCDF4 as nc
from scipy import interpolate, fft
from math import sqrt, ceil
from concurrent.futures import ThreadPoolExecutor
from scipy.fft import fft, ifft


def G(r):
    # G(r) ≈ 1/π * Σ [(-1)^n * (2n + 1) / ((2n + 1)^2 + r^2)^(3/2)] for n = 0 to 10000
    result = 0.0
    for n in range(10001):
        frac = (-1) ** n * (2 * n + 1) / ((2 * n + 1) ** 2 + r ** 2) ** 1.5
        result += frac
    return result / np.pi


def precalculate_sigma(h_min, h_max, Δx, Δy, n_h=20.0):
    Δh = (h_max - h_min) / 10
    if Δh == 0.0:
        Δh = 1.0
    l_h = np.arange(max(0.0001, h_min), h_max + Δh, Δh)

    def σ(h):
        σ_inv = 0.0
        n_max = ceil(n_h * h / Δx)
        m_max = ceil(n_h * h / Δy)

        N = np.arange(1, n_max + 1)
        M = np.arange(1, m_max + 1)

        for n in N:
            for m in M:
                σ_inv += 4 * G(sqrt((n * Δx) ** 2 + (m * Δy) ** 2) / h)

        for n in N:
            σ_inv += 2 * G(n * Δx / h)

        for m in M:
            σ_inv += 2 * G(m * Δy / h)

        σ_inv += G(0)

        return h ** 2 / (σ_inv * Δx * Δy)

    l_σ = [σ(h) for h in l_h]

    itp = interpolate.interp1d(l_h, l_σ, kind='cubic', fill_value='extrapolate')
    return itp


def apply_kajiura(b, d, η, h_min, h_max, Δx, Δy, water_level=0.0, n_h=20.0, σ=None):
    filter_nx_half = int(np.ceil(n_h * h_max / Δx / 2))
    filter_ny_half = int(np.ceil(n_h * h_max / Δy / 2))

    nx, ny = η.shape

    filter_matrix = np.zeros((ny, nx), dtype=float)

    assert filter_nx_half < nx
    assert filter_ny_half < ny

    for x in range(-filter_nx_half, filter_nx_half + 1):
        for y in range(-filter_ny_half, filter_ny_half + 1):
            filter_matrix[(ny + y) % ny, (nx + x) % nx] = G(sqrt((x * Δx) ** 2 + (y * Δy) ** 2) / h_max)

    η_diff = np.zeros_like(η)

    for x in range(nx):
        for y in range(ny):
            h_yx = max(0.0, water_level - b[y, x])

            if h_yx != 0.0:
                σ_h = σ(h_yx)
                η_diff[y, x] = σ_h * Δx * Δy / h_yx ** 2 * d[y, x]

    η_complex = fft.ifftshift(η_diff)
    Η = fft.fft2(η_complex)
    Filter = fft.fft2(filter_matrix)
    Η *= Filter
    η_complex = fft.ifft2(Η)
    η[:] = np.real(η_complex)


def main():
    in_filename = "input.nc"
    out_filename = "output.nc"
    t_end = -1

    print(f"Processing {t_end + 1} timesteps.")

    with nc.Dataset(in_filename, 'r') as ncfile:
        d = ncfile['d'][:]

        lx = ncfile['x'][:]
        ly = ncfile['y'][:]

        Δx = lx[1] - lx[0]
        Δy = ly[1] - ly[0]

        nx = len(lx)
        ny = len(ly)

        l_times = ncfile['time'][:]
        n_times = len(l_times)
        last_timestamp = l_times[-1]

        if t_end == -1 or t_end > n_times - 1:
            t_end = n_times - 1

        b = ncfile['b'][:]

        current_disp = np.zeros((ny, nx))
        current_η_diff = np.zeros((ny, nx))
        current_d_diff = np.zeros((ny, nx))

        print("Creating output file")
        if nc.Dataset(out_filename, 'r'):
            nc.Dataset(out_filename, 'r').close()
        timeatts = {"units": "seconds"}
        xatts = {"units": "m"}
        yatts = {"units": "m"}

        with nc.Dataset(out_filename, 'w') as outfile:
            outfile.createDimension("x", nx)
            outfile.createDimension("y", ny)
            outfile.createDimension("time", None)

            x = outfile.createVariable("x", "f4", ("x",))
            x.units = "m"
            x[:] = lx

            y = outfile.createVariable("y", "f4", ("y",))
            y.units = "m"
            y[:] = ly

            time = outfile.createVariable("time", "f4", ("time",))
            time.units = "seconds"
            time[:] = l_times[:t_end + 1]

            b_var = outfile.createVariable("b", "f4", ("y", "x"))
            b_var[:, :] = b

            eta_diff_var = outfile.createVariable("eta_diff", "f4", ("time", "y", "x"))
            d_diff_var = outfile.createVariable("d_diff", "f4", ("time", "y", "x"))

            with ThreadPoolExecutor() as executor:
                for t in range(t_end + 1):
                    print(f"Working on timestep {t} of {t_end}")
                    current_η_diff[:] = 0.0
                    current_d_diff[:] = d[:, :, t] - current_disp
                    apply_kajiura(b + current_disp, current_d_diff, current_η_diff, -np.max(b), -np.min(b), Δx, Δy)
                    current_disp = d[:, :, t]

                    print("Writing output for timestep")
                    eta_diff_var[t, :, :] = current_η_diff
                    d_diff_var[t, :, :] = current_d_diff


if __name__ == "__main__":
    main()
