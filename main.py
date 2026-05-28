import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import odr
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

KPC_TO_M = 3.0856775814913673e19
KM2S2_PER_KPC_TO_MS2 = 1e6 / KPC_TO_M

U_DISK = 0.50
U_BUL = 0.70


def read_sparc_massmodels(path: str) -> pd.DataFrame:
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith((
                "Title:", "Authors:", "Table:", "Byte-by-byte",
                "Bytes Format", "Note", "----", "================================================================================"
            )):
                continue

            parts = s.split()
            if len(parts) < 10:
                continue

            try:
                rows.append({
                    "ID": parts[0],
                    "D": float(parts[1]),
                    "R": float(parts[2]),
                    "Vobs": float(parts[3]),
                    "e_Vobs": float(parts[4]),
                    "Vgas": float(parts[5]),
                    "Vdisk": float(parts[6]),
                    "Vbul": float(parts[7]),
                    "SBdisk": float(parts[8]),
                    "SBbul": float(parts[9]),
                })
            except ValueError:
                continue

    return pd.DataFrame(rows)


def add_accelerations(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    gobs = (out["Vobs"] ** 2) / out["R"]
    gbar = (out["Vgas"] ** 2 + U_DISK * out["Vdisk"] ** 2 + U_BUL * out["Vbul"] ** 2) / out["R"]

    out["gobs"] = gobs * KM2S2_PER_KPC_TO_MS2
    out["gbar"] = gbar * KM2S2_PER_KPC_TO_MS2

    out["log_gobs"] = np.log10(out["gobs"])
    out["log_gbar"] = np.log10(out["gbar"])
    out["sigma_log_gobs"] = 2.0 * (out["e_Vobs"] / out["Vobs"]) / np.log(10)

    return out


def rar_model_log(beta, log_gbar):
    gdag = 10 ** beta[0]
    gbar = 10 ** log_gbar
    gobs = gbar / (1.0 - np.exp(-np.sqrt(gbar / gdag)))
    return np.log10(gobs)


def fit_rar_odr(df: pd.DataFrame, sigma_x_dex: float = 0.12):
    x = df["log_gbar"].to_numpy()
    y = df["log_gobs"].to_numpy()
    sx = np.full_like(x, sigma_x_dex, dtype=float)
    sy = df["sigma_log_gobs"].to_numpy()

    data = odr.RealData(x, y, sx=sx, sy=sy)
    model = odr.Model(rar_model_log)
    odr_obj = odr.ODR(data, model, beta0=[np.log10(1.2e-10)])
    out = odr_obj.run()

    gdag = 10 ** out.beta[0]
    return gdag, out


def binned_stats(x, y, nbins=16):
    x = np.asarray(x)
    y = np.asarray(y)
    edges = np.linspace(np.nanmin(x), np.nanmax(x), nbins + 1)

    x_mean, y_mean, y_rms = [], [], []
    for i in range(nbins):
        m = (x >= edges[i]) & (x < edges[i + 1])
        if np.sum(m) == 0:
            continue
        x_mean.append(np.mean(x[m]))
        y_mean.append(np.mean(y[m]))
        y_rms.append(np.std(y[m], ddof=1) if np.sum(m) > 1 else 0.0)

    return np.array(x_mean), np.array(y_mean), np.array(y_rms), edges


def add_density_mosaic(ax, x, y, xlim, ylim, bins=55, cmap="Blues"):
    """
    方形蓝色马赛克：二维直方图 + pcolormesh
    """
    H, xedges, yedges = np.histogram2d(x, y, bins=bins, range=[xlim, ylim])

    # 让颜色更像“密度块”，可按需要调整 norm
    H = np.ma.masked_where(H == 0, H)

    pcm = ax.pcolormesh(
        xedges, yedges, H.T,
        cmap=cmap,
        shading="auto",
        alpha=0.75
    )
    return pcm


def plot_rar_figure3_like(path: str, exclude_ids=None, nbins=16, mosaic_bins=55):
    df = read_sparc_massmodels(path)

    if exclude_ids is not None:
        exclude_ids = set(exclude_ids)
        df = df[~df["ID"].isin(exclude_ids)].copy()

    df = df[(df["Vobs"] > 0) & (df["e_Vobs"] / df["Vobs"] <= 0.10)].copy()
    df = add_accelerations(df)
    df = df[(df["gobs"] > 0) & (df["gbar"] > 0)].copy()

    gdag_fit, odr_out = fit_rar_odr(df, sigma_x_dex=0.12)

    # 这里的横纵坐标都直接用 log 值，不再设置 log scale
    x = df["log_gbar"].to_numpy()
    y = df["log_gobs"].to_numpy()
    resid = y - rar_model_log([np.log10(gdag_fit)], x)

    xlim = (x.min() - 0.15, x.max() + 0.15)
    ylim_top = (y.min() - 0.25, y.max() + 0.25)
    ylim_bot = (-0.65, 0.65)

    fig = plt.figure(figsize=(8.6, 9.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.25], hspace=0.06)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)

    # ========= 子图1：蓝色马赛克 + binned means =========
    add_density_mosaic(ax1, x, y, xlim, ylim_top, bins=mosaic_bins, cmap="Blues")

    x_mean, y_mean, y_rms, _ = binned_stats(x, y, nbins=nbins)

    ax1.errorbar(
        x_mean, y_mean, yerr=y_rms,
        fmt="s", ms=6.5, mfc="red", mec="red",
        ecolor="red", elinewidth=1.0, capsize=2, zorder=5
    )

    # 1:1 line
    xx = np.linspace(*xlim, 400)
    ax1.plot(xx, xx, ls=":", lw=1.0, c="0.3")

    # 拟合线
    yy_fit = rar_model_log([np.log10(gdag_fit)], 10 ** xx)
    ax1.plot(xx, yy_fit, lw=2.0, c="black")

    # ridge width
    ax1.plot(x_mean, y_mean + y_rms, ls="--", lw=1.0, c="black")
    ax1.plot(x_mean, y_mean - y_rms, ls="--", lw=1.0, c="black")

    ax1.set_xlim(*xlim)
    ax1.set_ylim(*ylim_top)
    ax1.set_ylabel(r"$\log_{10}(g_{\rm obs})$")

    # 右下角 inset：Residuals histogram
    axins = inset_axes(ax1, width="33%", height="30%", loc="lower right", borderpad=2.5)
    bins = np.linspace(resid.min(), resid.max(), 30)
    axins.hist(resid, bins=bins, color="0.65", edgecolor="none")
    mu = resid.mean()
    sig = resid.std(ddof=1)
    xxr = np.linspace(bins.min(), bins.max(), 300)
    pdf = np.exp(-0.5 * ((xxr - mu) / sig) ** 2)
    pdf = pdf / pdf.max() * axins.get_ylim()[1] * 0.9
    axins.plot(xxr, pdf, lw=1.1, c="black")
    axins.set_xlabel("Residuals", fontsize=8)
    axins.set_ylabel("Measurements", fontsize=8)
    axins.tick_params(axis="both", labelsize=8)



    # ========= 子图2：log(gbar) vs residuals =========
    # 先把 residual 按 log_gbar 分箱
    x_mean2, y_mean2, y_rms2, _ = binned_stats(x, resid, nbins=nbins)

    # 蓝色马赛克：二维密度图
    add_density_mosaic(ax2, x, resid, xlim, ylim_bot, bins=mosaic_bins, cmap="Blues")

    # 分箱误差棒：均值 ± rms
    ax2.errorbar(
        x_mean2, y_mean2, yerr=y_rms2,
        fmt="s", ms=6.5,
        mfc="red", mec="red",
        ecolor="red", elinewidth=1.0, capsize=2,
        zorder=5
    )

    # 画 dashed lines：mean ± rms
    ax2.plot(x_mean2, y_mean2 + y_rms2, ls="--", lw=1.0, c="black")
    ax2.plot(x_mean2, y_mean2 - y_rms2, ls="--", lw=1.0, c="black")

    # 中心线
    ax2.axhline(0.0, ls=":", lw=1.0, c="0.3")

    ax2.set_xlim(*xlim)
    ax2.set_ylim(*ylim_bot)
    ax2.set_xlabel(r"$\log_{10}(g_{\rm bar})$")
    ax2.set_ylabel(r"Residuals [dex]")


    fig.suptitle(rf"RAR Fig. 3 style, $g_\dagger \approx {gdag_fit:.2e}\,\mathrm{{m\,s^{{-2}}}}$", y=0.985)
    plt.show()

    return df, gdag_fit, odr_out


# 用法
file_path = "MassModels_Lelli2016c.mrt.txt"
df_used, gdag_fit, odr_out = plot_rar_figure3_like(file_path, exclude_ids=None)
print("Fitted gdag =", gdag_fit)
