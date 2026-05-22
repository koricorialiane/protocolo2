import math
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import streamlit as st
from scipy.special import j0

C = 299_792_458.0
KBM = -174.0
GSM_TIMESLOT_US = 577.0
GSM_CHANNEL_KHZ = 200.0

st.set_page_config(
    page_title="Protocolo Titán | GSM/EDGE Tactical Network",
    page_icon="📡",
    layout="wide",
    initial_sidebar_state="expanded",
)

CSS = """
<style>
[data-testid="stAppViewContainer"] {background: radial-gradient(circle at 15% 15%, #10284a 0, #07111f 36%, #030711 100%);}
.block-container {padding-top: 1.7rem; padding-bottom: 3rem; max-width: 1420px;}
.metric-card {border: 1px solid rgba(120,170,255,.22); background: linear-gradient(135deg, rgba(16,28,46,.94), rgba(7,17,31,.94)); border-radius: 18px; padding: 18px 20px; box-shadow: 0 10px 30px rgba(0,0,0,.25);}
.hero {border: 1px solid rgba(120,170,255,.26); border-radius: 24px; padding: 28px 32px; background: linear-gradient(135deg, rgba(23,50,88,.88), rgba(7,17,31,.78)); box-shadow: 0 20px 50px rgba(0,0,0,.32); margin-bottom: 18px;}
.hero h1 {font-size: 2.35rem; margin: 0 0 8px 0; letter-spacing: -0.02em;}
.hero p {font-size: 1.04rem; color: #C7D7F2; max-width: 980px;}
.badge {display:inline-block; border:1px solid rgba(78,161,255,.45); color:#B9DAFF; background:rgba(78,161,255,.10); padding:6px 10px; border-radius:999px; margin-right:8px; font-size:.82rem;}
.section-note {border-left: 4px solid #4EA1FF; padding: 12px 16px; background: rgba(78,161,255,.08); border-radius: 10px; color:#D7E8FF;}
.small {font-size: .9rem; color:#AFC3DE;}
hr {border-color: rgba(255,255,255,.12);}
</style>
"""
st.markdown(CSS, unsafe_allow_html=True)

@dataclass
class Params:
    fc_mhz: float
    tx_power_dbm: float
    tx_gain_dbi: float
    rx_gain_dbi: float
    cable_loss_db: float
    rx_sensitivity_dbm: float
    noise_figure_db: float
    cell_radius_a_km: float
    cell_radius_b_km: float
    hb_m: float
    hm_m: float
    n_cluster: int
    total_arfcn: int
    gos: float
    users_per_cell: int
    traffic_per_user: float
    reuse_tiers: int
    pathloss_exp: float
    shadow_sigma_db: float
    rician_k_db: float


def db_to_lin(x):
    return 10 ** (x / 10)


def lin_to_db(x):
    x = np.maximum(np.asarray(x), 1e-30)
    return 10 * np.log10(x)


def kmh_to_ms(v):
    return np.asarray(v) / 3.6


def doppler_hz(v_kmh, fc_mhz):
    return kmh_to_ms(v_kmh) * fc_mhz * 1e6 / C


def coherence_time_s(fd):
    return 0.423 / np.maximum(fd, 1e-12)


def fspl_db(d_km, fc_mhz):
    d_km = np.maximum(np.asarray(d_km), 1e-6)
    return 32.44 + 20 * np.log10(d_km) + 20 * np.log10(fc_mhz)


def hata_urban_db(d_km, fc_mhz, hb_m, hm_m):
    d_km = np.maximum(np.asarray(d_km), 1e-3)
    a_hm = (1.1 * np.log10(fc_mhz) - 0.7) * hm_m - (1.56 * np.log10(fc_mhz) - 0.8)
    return 69.55 + 26.16 * np.log10(fc_mhz) - 13.82 * np.log10(hb_m) - a_hm + (44.9 - 6.55 * np.log10(hb_m)) * np.log10(d_km)


def link_rx_dbm(p: Params, loss_db):
    return p.tx_power_dbm + p.tx_gain_dbi + p.rx_gain_dbi - p.cable_loss_db - loss_db


def thermal_noise_dbm(rbw_hz, nf_db):
    return KBM + 10 * np.log10(np.asarray(rbw_hz)) + nf_db


def erlang_b(A, m):
    B = 1.0
    for i in range(1, int(m) + 1):
        B = (A * B) / (i + A * B)
    return B


def inv_erlang_b(m, target_b):
    lo, hi = 0.0, max(1.0, m * 2.5)
    while erlang_b(hi, m) < target_b:
        hi *= 2
    for _ in range(80):
        mid = (lo + hi) / 2
        if erlang_b(mid, m) > target_b:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def c_i_hex_db(N, gamma=3.5, tiers=1):
    # Aproximación clásica: 6 interferentes primer anillo. Se amplía de forma conservadora con más anillos.
    D_over_R = math.sqrt(3 * N)
    total_i = 0.0
    for k in range(1, tiers + 1):
        interferers = 6 * k
        distance_factor = D_over_R * k
        total_i += interferers * (distance_factor ** (-gamma))
    return 10 * math.log10(1 / total_i)


def jakes_fading_db(v_kmh, fc_mhz, duration_ms=45, fs_hz=4000, kind="Rayleigh", k_db=6.0, seed=1):
    rng = np.random.default_rng(seed)
    fd = float(doppler_hz(v_kmh, fc_mhz))
    n = int(duration_ms / 1000 * fs_hz)
    t = np.arange(n) / fs_hz
    M = 48
    alpha = rng.uniform(0, 2 * np.pi, M)
    phi_i = rng.uniform(0, 2 * np.pi, M)
    phi_q = rng.uniform(0, 2 * np.pi, M)
    wi = 2 * np.pi * fd * np.cos(alpha)
    i = np.sum(np.cos(np.outer(t, wi) + phi_i), axis=1) / np.sqrt(M)
    q = np.sum(np.sin(np.outer(t, wi) + phi_q), axis=1) / np.sqrt(M)
    h = (i + 1j * q) / np.sqrt(2)
    if kind == "Rician":
        K = db_to_lin(k_db)
        los = np.exp(1j * 2 * np.pi * fd * t)
        h = np.sqrt(K / (K + 1)) * los + np.sqrt(1 / (K + 1)) * h
    amp = np.abs(h)
    amp = amp / np.sqrt(np.mean(amp**2))
    fade_db = 20 * np.log10(np.maximum(amp, 1e-6))
    return t * 1000, fade_db, fd


def arfcn_plan(total, N, first=1):
    arfcns = list(range(first, first + total))
    rows = []
    per_cell = math.floor(total / N)
    extras = total % N
    idx = 0
    for c in range(N):
        n = per_cell + (1 if c < extras else 0)
        cell = chr(ord('A') + c)
        assigned = arfcns[idx:idx+n]
        idx += n
        for j, a in enumerate(assigned):
            rows.append({
                "Celda": cell,
                "ARFCN": a,
                "Tipo": "BCCH" if j == 0 else "TCH / Datos EDGE",
                "Downlink MHz": 935.0 + 0.2 * a,
                "Uplink MHz": 890.0 + 0.2 * a,
                "Observación": "Potencia fija, sin hopping" if j == 0 else "Puede usar hopping controlado",
            })
    return pd.DataFrame(rows)


def cluster_points(N):
    # Coordenadas simples para visualización profesional del clúster.
    angles = np.linspace(0, 2*np.pi, N, endpoint=False)
    r = 1.0 if N <= 4 else 1.25
    pts = []
    for i, a in enumerate(angles):
        pts.append((chr(ord('A') + i), r*np.cos(a), r*np.sin(a)))
    return pd.DataFrame(pts, columns=["Celda", "x", "y"])


def sidebar() -> Params:
    st.sidebar.title("⚙️ Parámetros de ingeniería")
    with st.sidebar.expander("Radio GSM/EDGE", expanded=True):
        fc_mhz = st.slider("Frecuencia central [MHz]", 850.0, 960.0, 900.0, 1.0)
        tx_power_dbm = st.slider("Potencia BTS [dBm]", 20.0, 50.0, 43.0, 1.0)
        tx_gain_dbi = st.slider("Ganancia antena BTS [dBi]", 0.0, 18.0, 9.0, 0.5)
        rx_gain_dbi = st.slider("Ganancia terminal [dBi]", -3.0, 6.0, 0.0, 0.5)
        cable_loss_db = st.slider("Pérdidas feeder/conectores [dB]", 0.0, 6.0, 2.0, 0.5)
        rx_sensitivity_dbm = st.slider("Sensibilidad receptor [dBm]", -115.0, -85.0, -104.0, 1.0)
        noise_figure_db = st.slider("Figura de ruido analizador [dB]", 0.0, 15.0, 6.0, 0.5)
    with st.sidebar.expander("Escenarios", expanded=True):
        cell_radius_a_km = st.slider("Radio celda convoy A [km]", 0.5, 10.0, 3.0, 0.1)
        cell_radius_b_km = st.slider("Radio celda campamento B [km]", 0.3, 5.0, 1.5, 0.1)
        hb_m = st.slider("Altura antena BTS [m]", 10.0, 80.0, 35.0, 1.0)
        hm_m = st.slider("Altura terminal [m]", 1.0, 5.0, 1.5, 0.1)
        shadow_sigma_db = st.slider("Desvanecimiento lento σ [dB]", 0.0, 12.0, 6.0, 0.5)
        rician_k_db = st.slider("Factor K Rician [dB]", 0.0, 15.0, 6.0, 0.5)
    with st.sidebar.expander("Planificación y capacidad", expanded=True):
        n_cluster = st.selectbox("Tamaño de clúster N", [3, 4, 7, 9, 12], index=1)
        total_arfcn = st.slider("Portadoras disponibles", 4, 64, 24, 1)
        pathloss_exp = st.slider("Exponente de propagación γ para C/I", 2.5, 5.0, 3.7, 0.1)
        reuse_tiers = st.slider("Anillos co-canal considerados", 1, 3, 1, 1)
        gos = st.select_slider("GoS bloqueo Erlang B", options=[0.005, 0.01, 0.02, 0.05], value=0.02)
        users_per_cell = st.slider("Usuarios por celda", 10, 500, 120, 10)
        traffic_per_user = st.slider("Tráfico por usuario [Erlang]", 0.005, 0.12, 0.025, 0.005)
    return Params(fc_mhz, tx_power_dbm, tx_gain_dbi, rx_gain_dbi, cable_loss_db, rx_sensitivity_dbm, noise_figure_db,
                  cell_radius_a_km, cell_radius_b_km, hb_m, hm_m, n_cluster, total_arfcn, gos, users_per_cell,
                  traffic_per_user, reuse_tiers, pathloss_exp, shadow_sigma_db, rician_k_db)


def plot_template(fig, height=430):
    fig.update_layout(
        template="plotly_dark",
        height=height,
        margin=dict(l=30, r=30, t=55, b=35),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="Inter, Segoe UI, Arial", size=13),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    fig.update_xaxes(showgrid=True, gridcolor="rgba(255,255,255,.08)", zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor="rgba(255,255,255,.08)", zeroline=False)
    return fig


def mobility_tab(p: Params):
    st.subheader("Escenario A — movilidad extrema, Doppler y coherencia")
    speeds = np.linspace(1, 320, 500)
    fd = doppler_hz(speeds, p.fc_mhz)
    tc_us = coherence_time_s(fd) * 1e6

    base_speeds = np.array([50, 250])
    df = pd.DataFrame({
        "Velocidad [km/h]": base_speeds,
        "Velocidad [m/s]": kmh_to_ms(base_speeds),
        "Doppler máximo [Hz]": doppler_hz(base_speeds, p.fc_mhz),
        "Tiempo de coherencia [ms]": coherence_time_s(doppler_hz(base_speeds, p.fc_mhz))*1e3,
        "Tc / timeslot GSM": coherence_time_s(doppler_hz(base_speeds, p.fc_mhz))/(GSM_TIMESLOT_US*1e-6),
    })
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Doppler a 250 km/h", f"{df.loc[1,'Doppler máximo [Hz]']:.1f} Hz")
    c2.metric("Tc a 250 km/h", f"{df.loc[1,'Tiempo de coherencia [ms]']:.2f} ms")
    c3.metric("Timeslot GSM", f"{GSM_TIMESLOT_US:.0f} µs")
    c4.metric("Relación Tc/TS", f"{df.loc[1,'Tc / timeslot GSM']:.1f}×")

    col1, col2 = st.columns(2)
    fig1 = go.Figure()
    fig1.add_trace(go.Scatter(x=speeds, y=fd, name="fD máximo", line=dict(width=3)))
    fig1.add_trace(go.Scatter(x=base_speeds, y=doppler_hz(base_speeds, p.fc_mhz), mode="markers+text", text=["50", "250"], textposition="top center", name="Casos guía"))
    fig1.update_layout(title="Desviación Doppler máxima frente a velocidad")
    fig1.update_xaxes(title="Velocidad [km/h]")
    fig1.update_yaxes(title="fD [Hz]")
    col1.plotly_chart(plot_template(fig1), use_container_width=True)

    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(x=speeds, y=tc_us, name="Tc", line=dict(width=3)))
    fig2.add_hline(y=GSM_TIMESLOT_US, line_dash="dash", annotation_text="Timeslot GSM 577 µs")
    fig2.update_layout(title="Tiempo de coherencia comparado con el timeslot GSM")
    fig2.update_xaxes(title="Velocidad [km/h]")
    fig2.update_yaxes(title="Tc [µs]", type="log")
    col2.plotly_chart(plot_template(fig2), use_container_width=True)

    st.dataframe(df.style.format({"Velocidad [m/s]":"{:.2f}", "Doppler máximo [Hz]":"{:.2f}", "Tiempo de coherencia [ms]":"{:.3f}", "Tc / timeslot GSM":"{:.2f}"}), use_container_width=True)
    st.markdown("""
<div class='section-note'>
<b>Lectura de ingeniería:</b> aunque el Doppler aumenta linealmente con la velocidad, el tiempo de coherencia sigue siendo mayor que una ráfaga GSM. Eso permite defender que la ráfaga puede considerarse aproximadamente cuasiestática, pero a 250 km/h el margen temporal baja y conviene justificar codificación, entrelazado y diseño conservador.
</div>
""", unsafe_allow_html=True)


def fading_tab(p: Params):
    st.subheader("Canal multitrayecto — Rayleigh, Rician y correlación temporal")
    v = st.slider("Velocidad de simulación [km/h]", 10, 320, 250, 5)
    duration = st.slider("Duración de traza [ms]", 20, 200, 80, 10)
    col1, col2 = st.columns(2)
    for kind, col, seed in [("Rayleigh", col1, 4), ("Rician", col2, 8)]:
        t, fdb, fd_val = jakes_fading_db(v, p.fc_mhz, duration_ms=duration, kind=kind, k_db=p.rician_k_db, seed=seed)
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=t, y=fdb, name=kind, line=dict(width=2)))
        fig.add_hline(y=-10, line_dash="dot", annotation_text="Fade -10 dB")
        fig.update_layout(title=f"Traza de fading {kind} | fD={fd_val:.1f} Hz")
        fig.update_xaxes(title="Tiempo [ms]")
        fig.update_yaxes(title="Amplitud normalizada [dB]")
        col.plotly_chart(plot_template(fig, 390), use_container_width=True)

    tau_ms = np.linspace(0, 20, 500)
    fd_val = float(doppler_hz(v, p.fc_mhz))
    corr = j0(2*np.pi*fd_val*tau_ms/1000)
    figc = go.Figure()
    figc.add_trace(go.Scatter(x=tau_ms, y=corr, name="J0(2πfDτ)", line=dict(width=3)))
    figc.add_hline(y=0.5, line_dash="dash", annotation_text="correlación 0.5")
    figc.update_layout(title="Autocorrelación temporal aproximada del canal móvil")
    figc.update_xaxes(title="Retardo τ [ms]")
    figc.update_yaxes(title="Correlación")
    st.plotly_chart(plot_template(figc, 390), use_container_width=True)

    st.markdown("""
<div class='section-note'>
<b>Por qué ahora no son gráficas rectas:</b> el fading se genera con una suma de osciladores tipo Jakes. Rayleigh representa entorno sin línea de vista dominante; Rician añade una componente directa, por eso las caídas profundas son menos frecuentes.
</div>
""", unsafe_allow_html=True)


def link_budget_tab(p: Params):
    st.subheader("Presupuesto de enlace — FSPL, Okumura-Hata y margen operativo")
    d = np.linspace(0.05, max(8, p.cell_radius_a_km*1.6), 500)
    loss_fspl = fspl_db(d, p.fc_mhz)
    loss_hata = hata_urban_db(d, p.fc_mhz, p.hb_m, p.hm_m)
    rx_fspl = link_rx_dbm(p, loss_fspl)
    rx_hata = link_rx_dbm(p, loss_hata)
    margin_hata = rx_hata - p.rx_sensitivity_dbm - p.shadow_sigma_db

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=d, y=rx_fspl, name="Potencia recibida FSPL", line=dict(width=2)))
    fig.add_trace(go.Scatter(x=d, y=rx_hata, name="Potencia recibida Hata urbano", line=dict(width=3)))
    fig.add_hline(y=p.rx_sensitivity_dbm, line_dash="dash", annotation_text="Sensibilidad receptor")
    fig.update_layout(title="Potencia recibida frente a distancia")
    fig.update_xaxes(title="Distancia [km]")
    fig.update_yaxes(title="Prx [dBm]")
    st.plotly_chart(plot_template(fig), use_container_width=True)

    col1, col2 = st.columns(2)
    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(x=d, y=margin_hata, name="Margen con σ shadowing", fill="tozeroy", line=dict(width=3)))
    fig2.add_hline(y=0, line_dash="dash", annotation_text="Margen cero")
    fig2.update_layout(title="Margen de enlace Hata incluyendo reserva por shadowing")
    fig2.update_xaxes(title="Distancia [km]")
    fig2.update_yaxes(title="Margen [dB]")
    col1.plotly_chart(plot_template(fig2, 390), use_container_width=True)

    d_points = np.array([0.5, 1.5, 3.0, p.cell_radius_a_km])
    tbl = pd.DataFrame({
        "Distancia [km]": d_points,
        "FSPL [dB]": fspl_db(d_points, p.fc_mhz),
        "Hata urbano [dB]": hata_urban_db(d_points, p.fc_mhz, p.hb_m, p.hm_m),
        "Prx Hata [dBm]": link_rx_dbm(p, hata_urban_db(d_points, p.fc_mhz, p.hb_m, p.hm_m)),
        "Margen neto [dB]": link_rx_dbm(p, hata_urban_db(d_points, p.fc_mhz, p.hb_m, p.hm_m)) - p.rx_sensitivity_dbm - p.shadow_sigma_db,
    })
    col2.dataframe(tbl.style.format("{:.2f}"), use_container_width=True)


def planning_tab(p: Params):
    st.subheader("Escenario B — planificación espectral, reutilización e interferencia")
    D = p.cell_radius_b_km * math.sqrt(3*p.n_cluster)
    ci = c_i_hex_db(p.n_cluster, p.pathloss_exp, p.reuse_tiers)
    per_cell = p.total_arfcn / p.n_cluster
    traffic = p.users_per_cell * p.traffic_per_user
    trx_per_cell = max(1, math.floor(per_cell))
    # 1 BCCH TS + signaling. Aproximación conservadora: 7 TCH por TRX, restamos 1 canal de control total.
    voice_channels = max(1, trx_per_cell * 8 - 1)
    carried_capacity = inv_erlang_b(voice_channels, p.gos)

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("D/R", f"{math.sqrt(3*p.n_cluster):.2f}")
    c2.metric("D reutilización", f"{D:.2f} km")
    c3.metric("C/I estimado", f"{ci:.1f} dB")
    c4.metric("Carga/capacidad", f"{traffic:.1f}/{carried_capacity:.1f} Erl")

    plan = arfcn_plan(p.total_arfcn, p.n_cluster)
    st.dataframe(plan.style.format({"Downlink MHz":"{:.1f}", "Uplink MHz":"{:.1f}"}), use_container_width=True)

    col1, col2 = st.columns(2)
    pts = cluster_points(p.n_cluster)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=pts.x, y=pts.y, mode="markers+text", text=pts.Celda, textposition="middle center", marker=dict(size=70, opacity=.9), name="Celdas"))
    for _, r in pts.iterrows():
        fig.add_shape(type="circle", x0=r.x-.42, y0=r.y-.42, x1=r.x+.42, y1=r.y+.42, line=dict(width=1, color="rgba(78,161,255,.45)"))
    fig.update_layout(title="Visualización conceptual del clúster de frecuencias", xaxis_visible=False, yaxis_visible=False)
    col1.plotly_chart(plot_template(fig, 420), use_container_width=True)

    Ns = np.array([3,4,7,9,12,16,21])
    ci_vals = [c_i_hex_db(int(n), p.pathloss_exp, p.reuse_tiers) for n in Ns]
    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(x=Ns, y=ci_vals, mode="lines+markers", name="C/I", line=dict(width=3)))
    fig2.add_hline(y=9, line_dash="dot", annotation_text="Referencia GSM C/I≈9 dB")
    fig2.update_layout(title="Compromiso capacidad-interferencia al variar N")
    fig2.update_xaxes(title="Tamaño de clúster N")
    fig2.update_yaxes(title="C/I [dB]")
    col2.plotly_chart(plot_template(fig2, 420), use_container_width=True)

    if traffic <= carried_capacity:
        st.success("La capacidad estimada es suficiente para la carga ofrecida con el GoS seleccionado.")
    else:
        st.error("La carga ofrecida supera la capacidad Erlang B. Harían falta más TRX/portadoras, menor GoS o sectorización.")


def certification_tab(p: Params):
    st.subheader("Certificación y laboratorio — RBW, ruido integrado y barrido")
    rbw = np.array([100, 300, 1_000, 3_000, 10_000, 30_000, 100_000, 200_000, 300_000])
    noise = thermal_noise_dbm(rbw, p.noise_figure_db)
    sweep_rel = rbw.max()/rbw
    df = pd.DataFrame({"RBW [Hz]": rbw, "RBW [kHz]": rbw/1000, "Ruido integrado [dBm]": noise, "Tiempo relativo de barrido": sweep_rel})

    col1, col2 = st.columns(2)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=rbw, y=noise, mode="lines+markers", name="Ruido", line=dict(width=3)))
    fig.update_layout(title="Suelo de ruido del analizador frente a RBW")
    fig.update_xaxes(title="RBW [Hz]", type="log")
    fig.update_yaxes(title="Ruido [dBm]")
    col1.plotly_chart(plot_template(fig, 410), use_container_width=True)

    fig2 = go.Figure()
    fig2.add_trace(go.Bar(x=df["RBW [kHz]"].astype(str), y=df["Tiempo relativo de barrido"], name="Tiempo relativo"))
    fig2.update_layout(title="Coste experimental: menor RBW implica barrido más lento")
    fig2.update_xaxes(title="RBW [kHz]")
    fig2.update_yaxes(title="Tiempo relativo")
    col2.plotly_chart(plot_template(fig2, 410), use_container_width=True)

    st.dataframe(df.style.format({"RBW [Hz]":"{:.0f}", "RBW [kHz]":"{:.1f}", "Ruido integrado [dBm]":"{:.2f}", "Tiempo relativo de barrido":"{:.1f}×"}), use_container_width=True)

    st.markdown("""
<div class='section-note'>
<b>Criterio de certificación:</b> reducir RBW no reduce la señal real emitida; reduce el ruido integrado que ve el instrumento. Esto ayuda a detectar señales débiles, pero aumenta el tiempo de barrido y exige más cuidado con span, detector, VBW y promediado.
</div>
""", unsafe_allow_html=True)


def synthesis_tab(p: Params):
    st.subheader("Conclusión integrada y defensa técnica")
    fd250 = float(doppler_hz(250, p.fc_mhz))
    tc250_ms = float(coherence_time_s(fd250)*1e3)
    D = p.cell_radius_b_km * math.sqrt(3*p.n_cluster)
    ci = c_i_hex_db(p.n_cluster, p.pathloss_exp, p.reuse_tiers)
    noise100 = float(thermal_noise_dbm(100_000, p.noise_figure_db))
    noise1 = float(thermal_noise_dbm(1_000, p.noise_figure_db))

    st.markdown(f"""
<div class='hero'>
<h1>Síntesis de ingeniería</h1>
<p>Con una frecuencia central de <b>{p.fc_mhz:.0f} MHz</b>, a <b>250 km/h</b> aparece un Doppler máximo de <b>{fd250:.1f} Hz</b> y un tiempo de coherencia aproximado de <b>{tc250_ms:.2f} ms</b>. Como el timeslot GSM dura 0.577 ms, la ráfaga se mantiene razonablemente estable, aunque se recomienda diseño conservador.</p>
<p>En el campamento, con <b>N={p.n_cluster}</b>, la distancia de reutilización es <b>{D:.2f} km</b> y el C/I aproximado es <b>{ci:.1f} dB</b>. En laboratorio, bajar RBW de 100 kHz a 1 kHz desplaza el suelo de ruido de <b>{noise100:.1f} dBm</b> a <b>{noise1:.1f} dBm</b>, mejorando la visibilidad de emisiones débiles a costa de un barrido más lento.</p>
</div>
""", unsafe_allow_html=True)

    checklist = pd.DataFrame([
        ["Movilidad", "Doppler y Tc calculados; comparación explícita con timeslot", "Cumplido"],
        ["Canal", "Rayleigh/Rician, correlación temporal y fading rápido", "Cumplido"],
        ["Cobertura", "FSPL + Okumura-Hata + margen de shadowing", "Cumplido"],
        ["Planificación", "ARFCN, BCCH/TCH, clúster y distancia de reutilización", "Cumplido"],
        ["Interferencia", "C/I co-canal y compromiso N-capacidad", "Cumplido"],
        ["Capacidad", "Erlang B con usuarios, tráfico y GoS", "Cumplido"],
        ["Certificación", "RBW, ruido integrado y criterio de medida", "Cumplido"],
    ], columns=["Bloque", "Evidencia técnica", "Estado"])
    st.dataframe(checklist, use_container_width=True)


def main():
    p = sidebar()
    st.markdown("""
<div class='hero'>
<span class='badge'>GSM-900</span><span class='badge'>EDGE</span><span class='badge'>Movilidad</span><span class='badge'>Certificación RF</span>
<h1>Protocolo Titán — Simulador profesional de red táctica GSM/EDGE</h1>
<p>Dashboard técnico para diseñar, justificar y defender una red móvil clásica en escenarios críticos: convoy de alta velocidad y campamento base de emergencias. Incluye cálculos de movilidad, fading, presupuesto de enlace, reutilización, interferencia, capacidad y medidas de laboratorio.</p>
</div>
""", unsafe_allow_html=True)

    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "01 Movilidad", "02 Fading", "03 Enlace", "04 Planificación", "05 Certificación", "06 Síntesis"
    ])
    with tab1:
        mobility_tab(p)
    with tab2:
        fading_tab(p)
    with tab3:
        link_budget_tab(p)
    with tab4:
        planning_tab(p)
    with tab5:
        certification_tab(p)
    with tab6:
        synthesis_tab(p)

    st.divider()
    st.caption("Modelo académico-profesional: las fórmulas base son defendibles para el reto; las simulaciones aportan interpretación visual y criterio de ingeniería, no sustituyen una campaña real de medida.")

if __name__ == "__main__":
    main()
