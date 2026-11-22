# =====================================================================================
# Proyecto final de Analítica de Datos 2026-2S · Parte de IA (ítem 4: 4.2 a 4.5)
#
# Red funcional que reconstruye la curva anual de lluvia de una estación a partir de la
# anomalía de la curva satelital de su píxel en ese año:
#     anomalía satelital     X̃(t) = X(t) − X̄_s(t)   (X̄_s: curva satelital media de la estación en entrenamiento)
#     codificador funcional  z_k = ∫ w_k(t) X̃(t) dt = a_kᵀ W c̃
#     cabeza                 lineal | MLP | KAN   (make_head de la clase)
#     decodificador FPCA     corrección(t) = μ_r(t) + Σ_k ξ_k φ_k(t)   (FPCADecoder de la clase)
#     conexión residual      ŷ(t) = Ȳ_s(t) + corrección(t)
#                            (Ȳ_s: climatología de la estación, su curva media de entrenamiento;
#                             la red aprende cuánto se aparta cada año de esa curva media)
# Entradas: anomalía satelital y ONI (la altitud no entra: es fija por estación)
#
# Estructuras tomadas de los cuadernos de clase:
#   S14 (LSTM, GRU, Transformer, redes funcionales y KAN): set_seed, make_loader, count_params,
#       predict, train_model, KANLayer, fit_fpca, FPCADecoder, make_head, ForecasterBase,
#       save_checkpoint, KANRegressor (prueba sintética)
#   S13 (ARIMA, SARIMA, TCN): bucle de semillas y resumen media ± desviación
#   PDF de inferencia funcional: integrales por cuadratura con la regla del trapecio
# Entrada: la carpeta Salidas_imputacion/ del script de R, junto a este archivo.
# Convención del curso: comentarios en español; variables, clases y funciones en inglés.
# =====================================================================================

# ---- Librerías estándar ----
import os, json, time, copy, math, random, hashlib, itertools
from pathlib import Path

# ---- Cómputo numérico y datos ----
import numpy as np
import pandas as pd
from scipy.interpolate import BSpline

# ---- PyTorch ----
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader

# ---- Visualización ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3, "font.size": 10})

# ---- Configuración global ----
SEED = 42
SEEDS = [0, 1, 2]                                       # 3 semillas por configuración (rúbrica 4.3)
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "Salidas_imputacion"              # salidas del script de R (con relleno por Kriging)
OUT_DIR = BASE_DIR / "Salidas_IA"; OUT_DIR.mkdir(exist_ok=True)
FIG_DIR = BASE_DIR / "Figuras_IA"; FIG_DIR.mkdir(exist_ok=True)
CKPT_DIR = OUT_DIR / "checkpoints"; CKPT_DIR.mkdir(exist_ok=True)


def set_seed(seed: int = SEED) -> None:
    # Fija todas las fuentes de aleatoriedad para reproducibilidad
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


set_seed()
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if DEVICE.type == "cpu":
    torch.set_num_threads(max(1, os.cpu_count() or 1))
LOADER_KW = dict(num_workers=0, pin_memory=False)   # requerido en Windows / 4 GB VRAM

print(f"PyTorch {torch.__version__} | dispositivo: {DEVICE}")


# =====================================================================================
# 1. Utilidades de la clase (S14, celda "Utilidades compartidas")
# =====================================================================================
def make_loader(X, Y, batch_size=128, shuffle=False):
    ds = TensorDataset(torch.as_tensor(X), torch.as_tensor(Y))
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle, **LOADER_KW)


def count_params(model: nn.Module) -> int:
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def regression_metrics(y_true, y_pred, y_insample=None, season: int = 1) -> dict:
    # MAE, RMSE, sMAPE (%), MASE y R² sobre arreglos con la misma forma
    y_true = np.asarray(y_true, dtype=np.float64); y_pred = np.asarray(y_pred, dtype=np.float64)
    err = y_pred - y_true
    mae = np.mean(np.abs(err)); rmse = np.sqrt(np.mean(err ** 2))
    smape = 100 * np.mean(2 * np.abs(err) / (np.abs(y_true) + np.abs(y_pred) + 1e-8))
    r2 = 1 - np.sum(err ** 2) / np.sum((y_true - y_true.mean()) ** 2)
    out = {"MAE": mae, "RMSE": rmse, "sMAPE": smape, "R2": r2}
    if y_insample is not None:
        # MASE: error relativo al pronóstico ingenuo estacional en entrenamiento (<1 = mejor que ingenuo)
        scale = np.mean(np.abs(y_insample[season:] - y_insample[:-season]))
        out["MASE"] = mae / scale
    return out


@torch.no_grad()
def predict(model, loader):
    model.eval()
    preds, trues = [], []
    for xb, yb in loader:
        preds.append(model(xb.to(DEVICE)).cpu()); trues.append(yb)
    return torch.cat(preds).numpy(), torch.cat(trues).numpy()


def train_model(model, train_loader, val_loader, *, epochs=30, lr=1e-3, weight_decay=1e-4,
                patience=6, clip=1.0, scheduler="onecycle", loss_name="mse", verbose=True, params=None,
                loss_fn=None, penalty=None, lam=0.0):
    # Bucle de entrenamiento de la clase:
    #   AdamW (regularización L2 desacoplada) · scheduler · gradient clipping · early stopping
    # Agregado para el proyecto: loss_fn (pérdida L² funcional) y penalty·lam (rugosidad de w_k)
    model.to(DEVICE)
    params = list(params) if params is not None else [p for p in model.parameters() if p.requires_grad]
    opt = torch.optim.AdamW(params, lr=lr, weight_decay=weight_decay)
    if scheduler == "onecycle":
        sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=lr, epochs=epochs,
                                                    steps_per_epoch=len(train_loader), pct_start=0.3)
    elif scheduler == "plateau":
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=2)
    else:
        sched = None
    if loss_fn is None:
        loss_fn = nn.MSELoss() if loss_name == "mse" else nn.HuberLoss(delta=1.0)
    hist = {"train": [], "val": []}
    best_val, best_state, bad = float("inf"), None, 0
    for ep in range(1, epochs + 1):
        model.train(); run = 0.0
        for xb, yb in train_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            loss = loss_fn(model(xb), yb)
            if penalty is not None:
                loss = loss + lam * penalty()                 # penalización de rugosidad
            opt.zero_grad(set_to_none=True)
            loss.backward()
            nn.utils.clip_grad_norm_(params, clip)          # evita gradientes explosivos
            opt.step()
            if scheduler == "onecycle":
                sched.step()                                  # OneCycle avanza por batch
            run += loss.item() * len(xb)
        train_loss = run / len(train_loader.dataset)
        p, t = predict(model, val_loader)
        val_loss = loss_fn(torch.tensor(p), torch.tensor(t)).item()
        if scheduler == "plateau":
            sched.step(val_loss)                              # ¡no olvidar pasar la métrica!
        hist["train"].append(train_loss); hist["val"].append(val_loss)
        if val_loss < best_val - 1e-6:
            best_val, best_state, bad = val_loss, copy.deepcopy(model.state_dict()), 0
        else:
            bad += 1
        if verbose and (ep == 1 or ep % 5 == 0 or bad == 0):
            print(f"  época {ep:3d} | train {train_loss:.4f} | val {val_loss:.4f}{' *' if bad == 0 else ''}")
        if bad >= patience:
            if verbose: print(f"  early stopping en época {ep} (mejor val={best_val:.4f})")
            break
    model.load_state_dict(best_state)
    return hist


# =====================================================================================
# 2. Capa KAN de la clase (S14, Sección 3.5): FastKAN, bases radiales + SiLU
# =====================================================================================
class KANLayer(nn.Module):
    # Capa KAN estilo FastKAN: cada arista (i → o) tiene su función φ_oi(x) aprendible:
    #   φ_oi(x) = w_b[o,i]·SiLU(x) + Σ_g c[o,i,g]·exp(-((x-μ_g)/h)²)
    def __init__(self, in_features, out_features, grid_size=8, grid_range=(-2.0, 2.0),
                 use_layernorm=True, spline_init_std=0.1):
        super().__init__()
        self.in_features, self.out_features, self.grid_size = in_features, out_features, grid_size
        grid = torch.linspace(grid_range[0], grid_range[1], grid_size)
        self.register_buffer("grid", grid)                                  # centros μ_g (fijos)
        self.h = (grid_range[1] - grid_range[0]) / (grid_size - 1)          # separación entre centros
        self.norm = nn.LayerNorm(in_features) if use_layernorm else nn.Identity()
        self.base_weight = nn.Parameter(torch.empty(out_features, in_features))
        self.spline_weight = nn.Parameter(torch.randn(out_features, in_features, grid_size) * spline_init_std)
        nn.init.kaiming_uniform_(self.base_weight, a=math.sqrt(5))

    def basis(self, x):
        # x: (B, in) → bases radiales (B, in, G)
        return torch.exp(-((x[..., None] - self.grid) / self.h) ** 2)

    def forward(self, x):
        x = self.norm(x)
        base = F.silu(x) @ self.base_weight.T                                     # (B, out)
        spline = torch.einsum("big,oig->bo", self.basis(x), self.spline_weight)   # (B, out)
        return base + spline

    @torch.no_grad()
    def edge_function(self, out_idx, in_idx, x_grid):
        # Evalúa φ_{out,in} sobre una rejilla 1D (sin LayerNorm) para graficar
        b = self.basis(x_grid[:, None])[:, 0, :]                                  # (N, G)
        return self.base_weight[out_idx, in_idx] * F.silu(x_grid) + b @ self.spline_weight[out_idx, in_idx]

    def spline_l1(self):
        # Regularizador L1 sobre los coeficientes: favorece funciones simples (Liu et al., 2024)
        return self.spline_weight.abs().mean()


class KANRegressor(nn.Module):
    # KAN de una capa sin LayerNorm (las entradas ya están en la rejilla): ŷ = Σ_i φ_i(x_i)
    def __init__(self, n_in, grid_size=12, grid_range=(-2.0, 2.0)):
        super().__init__(); self.kan = KANLayer(n_in, 1, grid_size, grid_range, use_layernorm=False)
    def forward(self, x): return self.kan(x).squeeze(-1)


# =====================================================================================
# 3. FPCA y decodificador de la clase (S14, Sección 3.7)
# =====================================================================================
def bspline_basis(n_points, n_basis, degree=3):
    # Matriz B (n_points, n_basis) de B-splines cúbicos con nodos uniformes sobre [0, 1]
    x = np.linspace(0, 1, n_points)
    inner = np.linspace(0, 1, n_basis - degree + 1)
    knots = np.r_[[0.0] * degree, inner, [1.0] * degree]
    return BSpline.design_matrix(x, knots, degree).toarray()


def smoother_matrix(n_points, n_basis=12, lam=1e-2):
    # P-spline: c = (BᵀB + λ DᵀD)⁻¹ Bᵀ y  →  y_suave = S y  con  S = B (BᵀB + λ DᵀD)⁻¹ Bᵀ
    B = bspline_basis(n_points, n_basis)
    D = np.diff(np.eye(n_basis), n=2, axis=0)            # penalización de rugosidad (2ª diferencia)
    return B @ np.linalg.solve(B.T @ B + lam * D.T @ D, B.T)


def fit_fpca(curves, n_components, n_basis=12, lam=1e-2):
    # FPCA discretizada: suavizado P-spline → media μ(t) → SVD → autofunciones φ_k(t) y varianza explicada
    # En el proyecto se llama con n_basis=None: las curvas ya vienen suavizadas con la base del ítem 1 (R)
    Y = np.asarray(curves, dtype=np.float64)
    S = smoother_matrix(Y.shape[1], n_basis, lam) if n_basis else np.eye(Y.shape[1])
    Ys = Y @ S.T
    mean = Ys.mean(0)
    _, sing, Vt = np.linalg.svd(Ys - mean, full_matrices=False)
    Vt = Vt * np.sign(Vt.sum(1, keepdims=True) + 1e-12)   # fija el signo (la SVD es ambigua en ±)
    explained = sing ** 2 / np.sum(sing ** 2)
    return {"mean": mean.astype(np.float32), "basis": Vt[:n_components].astype(np.float32),
            "explained": explained, "all_basis": Vt.astype(np.float32)}


class FPCADecoder(nn.Module):
    # Decodificador funcional: puntajes ξ (B, K) → curva ŷ(t) = μ(t) + Σ_k ξ_k φ_k(t)  (B, H)
    # μ y φ son buffers: viajan en el state_dict pero NO se entrenan
    def __init__(self, horizon, n_components):
        super().__init__()
        self.register_buffer("mean", torch.zeros(horizon))
        self.register_buffer("basis", torch.zeros(n_components, horizon))

    def set_fpca(self, fpca):
        self.mean.copy_(torch.as_tensor(fpca["mean"])); self.basis.copy_(torch.as_tensor(fpca["basis"]))
        return self

    def forward(self, scores):
        return self.mean + scores @ self.basis

    def project(self, curves):
        # curva → puntajes (producto interno discreto con cada φ_k)
        return (curves - self.mean) @ self.basis.T


# =====================================================================================
# 4. Cabezas y modelo base de la clase (S14, Sección 4.1)
# =====================================================================================
def make_head(kind, in_features, out_features, hidden=32, dropout=0.1, kan_grid=8):
    # Cabeza de regresión intercambiable. 'mlp' y 'kan' comparten la proyección inicial:
    # la ÚNICA diferencia es el bloque no lineal (GELU fija vs funciones φ aprendidas)
    if kind == "linear":
        return nn.Linear(in_features, out_features)
    if kind == "mlp":
        return nn.Sequential(nn.Linear(in_features, hidden), nn.GELU(), nn.Dropout(dropout), nn.Linear(hidden, out_features))
    if kind == "kan":
        return nn.Sequential(nn.Linear(in_features, hidden), nn.Dropout(dropout), KANLayer(hidden, out_features, kan_grid))
    raise ValueError(kind)


class ForecasterBase(nn.Module):
    # Lógica común: encoder → cabeza → [decodificador FPCA]
    # (la RevIN de la clase no se usa: las curvas no son ventanas de una serie)
    def build_output(self, feat_dim, horizon, head, decoder, n_components, dropout, **head_kw):
        self.feat_dim, self.decoder_kind = feat_dim, decoder
        self.decoder = FPCADecoder(horizon, n_components) if decoder == "fpca" else None
        self.n_out = n_components if decoder == "fpca" else horizon       # la cabeza predice K puntajes o H valores
        self.head = make_head(head, feat_dim, self.n_out, dropout=dropout, **head_kw)

    def replace_head(self, kind, **kw):
        # Reemplazo de la cabeza de regresión conservando el encoder (y el decodificador)
        self.head = make_head(kind, self.feat_dim, self.n_out, **kw).to(next(self.parameters()).device)
        self.config["head"] = kind
        return self

    def encode(self, x):
        raise NotImplementedError

    def forward(self, x, return_scores=False):
        out = self.head(self.encode(x))                       # (B, K) puntajes o (B, H) valores
        y = self.decoder(out) if self.decoder is not None else out
        return (y, out) if return_scores else y


# =====================================================================================
# 5. Lo que pide el enunciado y la clase no trae: base del ítem 1, Gram por cuadratura,
#    codificador funcional y pérdida L² funcional
# =====================================================================================
def bspline_basis_knots(t, inner_knots, domain=(0.0, 365.0), degree=3):
    # Igual que bspline_basis de la clase (BSpline.design_matrix), pero con los nodos del ítem 1
    knots = np.r_[[domain[0]] * (degree + 1), inner_knots, [domain[1]] * (degree + 1)]
    return BSpline.design_matrix(np.clip(t, domain[0], domain[1] - 1e-9), knots, degree).toarray()


def trapezoid_weights(t):
    # Pesos de cuadratura w_m de la regla del trapecio en la malla t_1 < ... < t_M (PDF de inferencia)
    w = np.zeros(len(t))
    dt = np.diff(t)
    w[:-1] += dt / 2; w[1:] += dt / 2
    return w


def gram_matrix(B, w):
    # W = ∫ φ(t) φ(t)ᵀ dt ≈ Σ_m w_m φ(t_m) φ(t_m)ᵀ = Bᵀ D B,  D = diag(w)
    return B.T @ (w[:, None] * B)


class L2Functional(nn.Module):
    # Pérdida L² funcional: ‖ŷ − y‖² ≈ Σ_m w_m (ŷ(t_m) − y(t_m))²  (promedio sobre curvas)
    # Como ŷ y y viven en la base de estaciones, es igual a (ĉ − c)ᵀ W (ĉ − c)
    def __init__(self, weights):
        super().__init__()
        self.register_buffer("w", torch.as_tensor(weights, dtype=torch.float32))

    def forward(self, pred, true):
        return (((pred - true) ** 2) * self.w.to(pred.device)).sum(1).mean()


class FunctionalNet(ForecasterBase):
    # Red funcional: la entrada x = [coeficientes c̃ de la anomalía satelital | ONI | coeficientes de Ȳ_s]
    #   encode: z_k = a_kᵀ W c̃ = ∫ w_k(t) X̃(t) dt,  w_k(t) = Σ_j a_kj φ_j(t)
    #   salida: ŷ(t) = Ȳ_s(t) + μ_r(t) + Σ_k ξ_k φ_k(t)   (conexión residual: la red corrige la
    #           climatología de la estación; el decodificador FPCA reconstruye esa corrección)
    #   Los coeficientes de Ȳ_s viajan en la entrada solo para sumar la curva base: no pasan
    #   por el codificador ni por la cabeza
    def __init__(self, n_basis_x, gram, n_basis_y, basis_y, n_neurons=6, n_scalars=1, horizon=365,
                 head="linear", n_components=6, dropout=0.1, hidden=32):
        super().__init__()
        self.config = dict(n_basis_x=n_basis_x, n_basis_y=n_basis_y, n_neurons=n_neurons, n_scalars=n_scalars,
                           horizon=horizon, head=head, n_components=n_components, dropout=dropout, hidden=hidden)
        self.n_basis_x, self.n_scalars, self.n_basis_y = n_basis_x, n_scalars, n_basis_y
        self.register_buffer("gram", torch.as_tensor(gram, dtype=torch.float32))
        # Para la conexión residual: base de estaciones en la malla diaria
        self.register_buffer("basis_y", torch.as_tensor(basis_y, dtype=torch.float32))   # (K_Y, H)
        self.weight_coef = nn.Parameter(torch.randn(n_neurons, n_basis_x) * 0.1)   # a_k
        # Penalización de rugosidad P-spline de la clase (segundas diferencias de los coeficientes)
        self.register_buffer("D2", torch.as_tensor(np.diff(np.eye(n_basis_x), n=2, axis=0), dtype=torch.float32))
        kw = {} if head == "linear" else dict(hidden=hidden)
        self.build_output(n_neurons + n_scalars, horizon, head, "fpca", n_components, dropout, **kw)

    def encode(self, x):
        c, s = x[:, :self.n_basis_x], x[:, self.n_basis_x:self.n_basis_x + self.n_scalars]
        z = c @ (self.weight_coef @ self.gram).T            # z = A W c
        return torch.cat([z, s], dim=1)

    def roughness(self):
        # Σ_k ‖D₂ a_k‖²: rugosidad de las funciones de pesos w_k(t)
        return ((self.weight_coef @ self.D2.T) ** 2).sum()

    def base_curve(self, x):
        # Ȳ_s(t) en la malla diaria a partir de sus coeficientes (últimas columnas de la entrada)
        return x[:, -self.n_basis_y:] @ self.basis_y

    def forward(self, x, return_scores=False):
        # Conexión residual: ŷ(t) = Ȳ_s(t) + corrección(t)
        y, scores = super().forward(x, return_scores=True)
        y = y + self.base_curve(x)
        return (y, scores) if return_scores else y


# =====================================================================================
# 6. Checkpoints de la clase (S14, Sección 4.2)
# =====================================================================================
def fixed_batch(model):
    # Lote fijo de huella con el tamaño de entrada del modelo (K_X coeficientes + ONI + K_Y
    # coeficientes de la climatología), que depende de las bases que eligió R
    n_in = model.config["n_basis_x"] + model.config["n_scalars"] + model.config["n_basis_y"]
    return torch.randn(8, n_in, generator=torch.Generator().manual_seed(999))


@torch.no_grad()
def model_fingerprint(model, batch=None):
    model.eval()
    batch = fixed_batch(model) if batch is None else batch
    preds = model(batch.to(DEVICE)).cpu().numpy().round(5)              # redondeo: estable entre máquinas
    return hashlib.sha256(preds.tobytes()).hexdigest()


def save_checkpoint(model, path, meta=None, batch=None):
    torch.save({"class": type(model).__name__, "config": model.config, "state_dict": model.state_dict(),
                "fingerprint": model_fingerprint(model, batch), "meta": meta or {}}, path)


# =====================================================================================
# 7. Métricas del enunciado (4.4) en mm: error L² relativo, KGE, sesgo del acumulado anual
#    y error en la fecha de los dos picos (días). Las oficiales con bootstrap y Wilcoxon se
#    calculan en R.
# =====================================================================================
def curve_metrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=np.float64); y_pred = np.asarray(y_pred, dtype=np.float64)
    l2 = np.linalg.norm(y_pred - y_true, axis=1) / np.linalg.norm(y_true, axis=1)
    r = np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(y_pred, y_true)])
    alpha = y_pred.std(1) / y_true.std(1); beta = y_pred.mean(1) / y_true.mean(1)
    kge = 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
    bias = 100 * (y_pred.sum(1) - y_true.sum(1)) / y_true.sum(1)
    half = 181                                             # pico 1 en ene-jun, pico 2 en jul-dic
    p1 = np.abs(y_pred[:, :half].argmax(1) - y_true[:, :half].argmax(1))
    p2 = np.abs(y_pred[:, half:].argmax(1) - y_true[:, half:].argmax(1))
    return pd.DataFrame({"L2_rel": l2, "KGE": kge, "sesgo_pct": bias, "err_pico1": p1, "err_pico2": p2})


def interannual_metrics(y_true, y_pred, station_ids):
    # ¿El modelo capta cómo cambia cada año, o solo la curva media de cada estación?
    #   alfa_intra: variabilidad de la curva dentro del año, predicha / observada
    #   alfa_entre: variabilidad entre años (anomalías respecto a la curva media de cada estación), predicha / observada
    #   corr_anom_curva: correlación de las anomalías día a día
    #   corr_anom_total: correlación de las anomalías del total anual (¿acierta qué año fue más lluvioso o seco?)
    #   pico_pred_obs: altura del pico anual predicho / observado
    y_true = np.asarray(y_true, dtype=np.float64); y_pred = np.asarray(y_pred, dtype=np.float64)
    a_t, a_p = np.zeros_like(y_true), np.zeros_like(y_pred)
    keep = np.zeros(len(y_true), dtype=bool)
    for s_ in np.unique(station_ids):
        i = station_ids == s_
        if i.sum() >= 2:                                    # se necesitan al menos 2 años por estación
            a_t[i] = y_true[i] - y_true[i].mean(0); a_p[i] = y_pred[i] - y_pred[i].mean(0); keep |= i
    tot_t, tot_p = a_t[keep].sum(1), a_p[keep].sum(1)
    return {"alfa_intra": np.mean(y_pred.std(1) / y_true.std(1)),
            "alfa_entre": a_p[keep].std() / a_t[keep].std(),
            "corr_anom_curva": np.corrcoef(a_p[keep].ravel(), a_t[keep].ravel())[0, 1] if a_p[keep].std() > 0 else 0.0,
            "corr_anom_total": np.corrcoef(tot_p, tot_t)[0, 1] if tot_p.std() > 0 else 0.0,
            "pico_pred_obs": np.mean(y_pred.max(1) / y_true.max(1))}


# =====================================================================================
# 8. Programa principal
# =====================================================================================
if __name__ == "__main__":
    t_start = time.time()

    # ---- 8.1 Datos del ítem 1 (R) ----
    spec = json.loads((DATA_DIR / "base_bspline.json").read_text(encoding="utf-8"))
    dom = tuple(spec["dominio"])
    knots_y, knots_x = spec["estaciones"]["nodos_interiores"], spec["satelite"]["nodos_interiores"]
    K_Y, K_X = spec["estaciones"]["K"], spec["satelite"]["K"]
    days = np.arange(365) + 0.5                                    # malla diaria (la misma de R)
    By = bspline_basis_knots(days, knots_y, dom); Bx = bspline_basis_knots(days, knots_x, dom)

    # Verificación: la base de Python reproduce la curva de control de R
    control = pd.read_csv(DATA_DIR / "curva_control_satelite.csv")
    coef_x = pd.read_csv(DATA_DIR / "coeficientes_satelite.csv")
    cx_cols = [f"c{j + 1}" for j in range(K_X)]
    diff_control = np.abs(bspline_basis_knots(control.t.values, knots_x, dom) @ coef_x.iloc[0][cx_cols].values.astype(float)
                          - control.valor.values).max()
    print(f"Base B-spline contra R (curva de control): diferencia máxima {diff_control:.1e}")

    pairs = pd.read_csv(DATA_DIR / "particion.csv")
    coef_y = pd.read_csv(DATA_DIR / "coeficientes_estaciones.csv")
    oni = pd.read_csv(DATA_DIR / "oni_anual.csv")
    cy_cols = [f"c{j + 1}" for j in range(K_Y)]
    df = (pairs[["STATION_ID", "YEAR", "pixel_id", "Altitud", "Conjunto"]]
          .merge(coef_y[["STATION_ID", "YEAR"] + cy_cols], on=["STATION_ID", "YEAR"])
          .merge(coef_x[["pixel_id", "YEAR"] + cx_cols].rename(columns={c: "x" + c for c in cx_cols}), on=["pixel_id", "YEAR"])
          .merge(oni[["YEAR", "ONI", "Fase"]], on="YEAR"))
    print(f"Pares estación-año: {len(df)} | {df.Conjunto.value_counts().to_dict()}")

    # Curvas en mm/día sobre la malla diaria; se recortan a 0 los valores negativos del suavizado
    Y_mm = np.maximum(df[cy_cols].values @ By.T, 0)                # estación (observada)
    X_mm = np.maximum(df[["x" + c for c in cx_cols]].values @ Bx.T, 0)   # satélite sin corregir
    print(f"Días recortados a 0 en las curvas de R: estaciones {100 * (df[cy_cols].values @ By.T < 0).mean():.2f} %")

    is_tr = (df.Conjunto == "entrenamiento").values
    is_va = (df.Conjunto == "validacion").values
    is_te = (df.Conjunto == "prueba").values

    # ---- 8.2 Cuadratura (regla del trapecio) y matriz de Gram ----
    t_fine = np.linspace(dom[0], dom[1], 3651)                     # malla fina de 0.1 días
    W_x = gram_matrix(bspline_basis_knots(t_fine, knots_x, dom), trapezoid_weights(t_fine))
    w_days = trapezoid_weights(days)                               # pesos de la pérdida L² en la malla diaria
    W_y = gram_matrix(By, w_days)
    print(f"Gram satélite: 1ᵀW1 = {W_x.sum():.3f} (longitud del dominio = {dom[1] - dom[0]:.0f})")

    # ---- 8.3 Escala: √P o log(1+P) (el enunciado pide modelar en una de las dos) ----
    TRANSFORMS = {"sqrt": (np.sqrt, np.square), "log1p": (np.log1p, np.expm1)}

    def to_scale(name):
        # La transformación se aplica a la curva suavizada en mm y se vuelve a proyectar en la
        # misma base (mínimos cuadrados, como en el ítem 1)
        f = TRANSFORMS[name][0]
        cy = np.linalg.solve(By.T @ By, By.T @ f(Y_mm).T).T
        cx = np.linalg.solve(Bx.T @ Bx, Bx.T @ f(X_mm).T).T
        return cy, cx

    # ---- 8.4 Conjuntos: entrada x = [anomalía satelital | ONI estandarizado | coef. de Ȳ_s] ----
    # La altitud NO entra: es fija por estación y la red la usaba como identificador de estación.
    # La parte fija de cada estación queda en su climatología Ȳ_s (curva media de sus años de
    # entrenamiento) y la red solo ve lo que cambia cada año: la anomalía del satélite y el ONI
    def station_means(coef):
        # Media de los coeficientes por estación, SOLO con años de entrenamiento; una estación sin
        # años de entrenamiento usa la media global de entrenamiento
        st = df.STATION_ID.values
        tab = pd.DataFrame(coef[is_tr]).groupby(st[is_tr]).mean()
        glob = coef[is_tr].mean(0)
        return np.vstack([tab.loc[s_].values if s_ in tab.index else glob for s_ in st])

    def build_sets(name):
        cy, cx = to_scale(name)
        cy_clim, cx_clim = station_means(cy), station_means(cx)
        sc = df[["ONI"]].values.astype(float)
        sc_mu, sc_sd = sc[is_tr].mean(0), sc[is_tr].std(0)          # SOLO entrenamiento
        X = np.hstack([cx - cx_clim, (sc - sc_mu) / sc_sd, cy_clim]).astype(np.float32)
        Y = (cy @ By.T).astype(np.float32)                         # curva objetivo en la malla diaria
        Xs = (cy_clim @ By.T).astype(np.float32)                   # curva base Ȳ_s en la misma escala
        return X, Y, cy, cx, Xs

    def inverse(name, y):
        return TRANSFORMS[name][1](np.maximum(y, 0))               # vuelta a mm, no negativa

    def make_net(head, fp, n_neurons=6, hidden=32):
        # Red funcional con conexión residual; el decodificador FPCA se ajusta a la corrección
        m = FunctionalNet(K_X, W_x, K_Y, By.T, n_neurons=n_neurons, head=head, n_components=fp["basis"].shape[0], hidden=hidden)
        m.decoder.set_fpca(fp)                                     # μ_r y φ_k fijos, solo entrenamiento
        return m

    # Elección de la escala con la cabeza lineal y validación (nunca con prueba)
    scale_rows = []
    for name in TRANSFORMS:
        X, Y, cy, cx, Xs = build_sets(name)
        fp = fit_fpca(Y[is_tr] - Xs[is_tr], 6, n_basis=None)       # FPCA de la corrección (estación − climatología)
        set_seed(SEED)
        m = make_net("linear", fp)
        train_model(m, make_loader(X[is_tr], Y[is_tr], 64, shuffle=True), make_loader(X[is_va], Y[is_va], 512),
                    epochs=150, lr=3e-3, patience=15, verbose=False, loss_fn=L2Functional(w_days),
                    penalty=m.roughness, lam=0.1)
        p, _ = predict(m, make_loader(X[is_va], Y[is_va], 512))
        y_rt = inverse(name, Y[is_tr])                              # ida y vuelta de la curva verdadera
        scale_rows.append({"escala": name, "L2_rel_val_lineal": curve_metrics(Y_mm[is_va], inverse(name, p)).L2_rel.mean(),
                           "ida_vuelta_L2_rel_mediano": np.median(np.linalg.norm(y_rt - Y_mm[is_tr], axis=1) / np.linalg.norm(Y_mm[is_tr], axis=1))})
    scale_df = pd.DataFrame(scale_rows)
    SCALE = scale_df.loc[scale_df.L2_rel_val_lineal.idxmin(), "escala"]
    print(scale_df.round(4).to_string(index=False)); print(f"Escala elegida: {SCALE}")
    scale_df.to_csv(OUT_DIR / "eleccion_escala.csv", index=False)

    X, Y, cy, cx, Xs = build_sets(SCALE)
    R_tr = Y[is_tr] - Xs[is_tr]                                    # corrección a aprender: anomalía de la estación

    # ---- 8.5 FPCA de la corrección, solo con entrenamiento, y elección de K (piso, S14 8.1) ----
    fpca_rows = []
    for K in [4, 5, 6]:
        f = fit_fpca(R_tr, K, n_basis=None)
        dec = FPCADecoder(365, K).set_fpca(f)
        rec = Xs[is_va] + dec(dec.project(torch.tensor(Y[is_va] - Xs[is_va]))).numpy()   # oráculo
        fpca_rows.append({"K": K, "var_explicada_acum": f["explained"][:K].sum(),
                          "piso_L2_rel_val (oráculo)": curve_metrics(Y_mm[is_va], inverse(SCALE, rec)).L2_rel.mean()})
    fpca_tab = pd.DataFrame(fpca_rows).set_index("K")
    print(fpca_tab.round(4))
    fpca_tab.to_csv(OUT_DIR / "fpca_eleccion_K.csv")
    K_FPCA = 6                                                     # máximo permitido: menor piso
    fpca = fit_fpca(R_tr, K_FPCA, n_basis=None)

    # Verificación: la pérdida en la malla con pesos del trapecio = (ĉ − c)ᵀ W (ĉ − c)
    # ŷ − y vive en la base de estaciones (la climatología y la corrección también); se usa la
    # base con la unión de los nodos de estaciones y satélite, que la contiene, y W es su Gram
    knots_u = sorted(set(knots_y) | set(knots_x))
    Bu = bspline_basis_knots(days, knots_u, dom)
    W_u = gram_matrix(Bu, w_days)
    dec = FPCADecoder(365, K_FPCA).set_fpca(fpca)
    y_hat = (Xs[is_va][:50] + dec(dec.project(torch.tensor(Y[is_va][:50] - Xs[is_va][:50]))).numpy()).astype(np.float64)
    e = np.linalg.solve(Bu.T @ Bu, Bu.T @ (y_hat - Y[is_va][:50].astype(np.float64)).T).T
    loss_grid = L2Functional(w_days)(torch.tensor(y_hat), torch.tensor(Y[is_va][:50].astype(np.float64))).item()
    loss_W = np.mean(np.einsum("nj,jl,nl->n", e, W_u, e))
    print(f"Pérdida en la malla {loss_grid:.6f} | (ĉ − c)ᵀW(ĉ − c) {loss_W:.6f}")

    # ---- 8.5b Número de neuronas funcionales (6 a 8): cabeza lineal, 3 semillas, validación ----
    # Regla de la clase (S14, 7.1): si dos configuraciones difieren menos que su desviación entre
    # semillas, la diferencia no es concluyente y se prefiere la más simple
    neuron_rows = []
    for n_neu in [6, 7, 8]:
        for seed in SEEDS:
            set_seed(seed)
            m = make_net("linear", fpca, n_neurons=n_neu)
            h = train_model(m, make_loader(X[is_tr], Y[is_tr], 64, shuffle=True), make_loader(X[is_va], Y[is_va], 512),
                            epochs=150, lr=3e-3, patience=15, verbose=False, loss_fn=L2Functional(w_days),
                            penalty=m.roughness, lam=10.0)
            neuron_rows.append({"neuronas": n_neu, "semilla": seed, "perdida_val": min(h["val"])})
    neuron_df = pd.DataFrame(neuron_rows).groupby("neuronas").perdida_val.agg(["mean", "std"])
    best_n = neuron_df["mean"].idxmin()
    N_NEURONS = int(min(n for n in neuron_df.index
                        if neuron_df.loc[n, "mean"] - neuron_df.loc[best_n, "mean"]
                        < max(neuron_df.loc[n, "std"], neuron_df.loc[best_n, "std"])))
    neuron_df.to_csv(OUT_DIR / "eleccion_neuronas.csv")
    print(neuron_df.round(3)); print(f"Neuronas funcionales elegidas: {N_NEURONS}")

    # ---- 8.6 Prueba de la KAN con una función sintética (S14, Sección 3.8, sin ruido) ----
    set_seed()
    n = 3000
    Xtoy = torch.rand(n, 3) * 4 - 2                          # x ~ U(-2, 2)
    ytoy = torch.sin(math.pi * Xtoy[:, 0]) + Xtoy[:, 1] ** 2
    tr_toy = make_loader(Xtoy[:2400].numpy(), ytoy[:2400].numpy(), 128, shuffle=True)
    va_toy = make_loader(Xtoy[2400:].numpy(), ytoy[2400:].numpy(), 512)
    # Rejilla de centros un poco más ancha que el rango de x (−2.5 a 2.5): con la rejilla en
    # (−2, 2) el error se concentra en los bordes (|x| > 1.8) y el RMSE no baja de 0.02
    kan_toy = KANRegressor(3, grid_size=16, grid_range=(-2.5, 2.5))
    _ = train_model(kan_toy, tr_toy, va_toy, epochs=400, lr=1e-2, patience=40, verbose=False)
    p, t = predict(kan_toy, va_toy)
    rmse_toy = float(np.sqrt(np.mean((p - t) ** 2)))
    print(f"KAN sintética sin(πx₁) + x₂² (+ 0·x₃): RMSE de validación = {rmse_toy:.4f} (umbral 1e-2)")

    # ---- 8.7 Ablación de la cabeza: mismo codificador, mismo decodificador, mismo presupuesto ----
    tr_loader = lambda: make_loader(X[is_tr], Y[is_tr], 64, shuffle=True)
    va_loader = make_loader(X[is_va], Y[is_va], 512)
    te_loader = make_loader(X[is_te], Y[is_te], 512)

    # Parámetros: la KAN usa la configuración de la clase (hidden=32, kan_grid=8); el MLP se
    # ajusta para tener el mismo número de parámetros (±10 %)
    n_kan = count_params(make_net("kan", fpca, n_neurons=N_NEURONS).head)
    N_IN = N_NEURONS + 1                                           # z_1..z_K y el ONI
    H_MLP = round((n_kan - K_FPCA) / (N_IN + 1 + K_FPCA))
    HIDDEN = {"linear": None, "mlp": H_MLP, "kan": 32}

    def new_model(kind):
        return make_net(kind, fpca, n_neurons=N_NEURONS, hidden=HIDDEN[kind] or 32)

    param_df = pd.DataFrame({kind: {"params_cabeza": count_params(new_model(kind).head),
                                    "params_total": count_params(new_model(kind)),
                                    "hidden": HIDDEN[kind]} for kind in HIDDEN}).T
    param_df["dif_vs_KAN_pct"] = 100 * (param_df.params_cabeza - param_df.loc["kan", "params_cabeza"]) / param_df.loc["kan", "params_cabeza"]
    print(param_df)
    param_df.to_csv(OUT_DIR / "parametros.csv")

    # Búsqueda de hiperparámetros: mismas pruebas, épocas y paciencia para las tres cabezas
    grid_hp = list(itertools.product([1e-3, 3e-3], [0.1, 10.0]))      # (lr, λ de rugosidad)
    EPOCHS, PATIENCE = 150, 15
    rows, models = [], {}
    t0 = time.time()
    for kind in HIDDEN:
        for j, (lr, lam) in enumerate(grid_hp):
            for seed in SEEDS:
                set_seed(seed)
                m = new_model(kind)
                h = train_model(m, tr_loader(), va_loader, epochs=EPOCHS, lr=lr, patience=PATIENCE, verbose=False,
                                loss_fn=L2Functional(w_days), penalty=m.roughness, lam=lam)
                p, _ = predict(m, va_loader)
                rows.append({"cabeza": kind, "prueba": j, "lr": lr, "lam": lam, "semilla": seed,
                             "perdida_val": min(h["val"]), "epocas": len(h["val"]),
                             "L2_rel_val_mm": curve_metrics(Y_mm[is_va], inverse(SCALE, p)).L2_rel.mean()})
                models[(kind, j, seed)] = m
        print(f"  {kind}: listo ({time.time() - t0:.0f} s)")
    search_df = pd.DataFrame(rows)
    search_df.to_csv(OUT_DIR / "busqueda_hiperparametros.csv", index=False)
    summary = search_df.groupby(["cabeza", "prueba", "lr", "lam"]).perdida_val.agg(["mean", "std"]).reset_index()
    best = summary.loc[summary.groupby("cabeza")["mean"].idxmin()].set_index("cabeza")
    print(best.round(4))
    final = {kind: [models[(kind, int(best.loc[kind, "prueba"]), s)] for s in SEEDS] for kind in HIDDEN}

    # ---- 8.8 Equivalencia de la cabeza lineal con la regresión funcional ----
    # ŷ(t) = Ȳ_s(t) + μ_r(t) + Σ_k ξ_k φ_k(t), ξ = B z + b, z_k = ∫ w_k(s) X̃(s) ds  ⇒
    # ŷ(t) − Ȳ_s(t) = α(t) + ∫ β(s,t) X̃(s) ds + γ(t)·ONI,  β(s,t) = Σ_k Σ_m w_k(s) B_mk φ_m(t)
    # (la cabeza lineal es una regresión funcional de la anomalía de la estación sobre la anomalía satelital)
    Bx_fine = bspline_basis_knots(t_fine, knots_x, dom); w_fine = trapezoid_weights(t_fine)
    eq_diffs = []
    for m in final["linear"]:
        A = m.weight_coef.detach().double().numpy(); Bh = m.head.weight.detach().double().numpy()
        b0 = m.head.bias.detach().double().numpy()
        mu_t = fpca["mean"].astype(np.float64); phi = fpca["basis"].astype(np.float64)
        w_s = Bx_fine @ A.T                                         # w_k(s) en la malla fina
        beta = w_s @ Bh[:, :N_NEURONS].T @ phi                      # β(s, t): (n_s, 365)
        X_s = X[:, :K_X].astype(np.float64) @ Bx_fine.T             # X̃(s) en la malla fina
        alpha = mu_t + b0 @ phi
        gamma = Bh[:, N_NEURONS:].T @ phi                           # γ(t) del ONI
        Y_base = X[:, -K_Y:].astype(np.float64) @ By.T              # Ȳ_s(t): término de la conexión residual
        y_reg = Y_base + alpha + (X_s * w_fine) @ beta + X[:, K_X:K_X + 1].astype(np.float64) @ gamma
        m64 = copy.deepcopy(m).double()
        with torch.no_grad():
            y_net = m64(torch.tensor(X, dtype=torch.float64)).numpy()
        eq_diffs.append(np.abs(y_net - y_reg).max())
    print(f"Equivalencia cabeza lineal = regresión funcional: diferencia máxima {max(eq_diffs):.1e} (umbral 1e-4)")

    # ---- 8.9 Evaluación en mm: media ± desviación entre semillas (patrón de S13) ----
    # Modelo de referencia: la climatología de cada estación (su curva media de entrenamiento).
    # Si una cabeza no le gana, solo está aprendiendo la media y no la variación de cada año.
    clim_station = pd.DataFrame(Y_mm[is_tr]).groupby(df.loc[is_tr, "STATION_ID"].values).mean()
    clim_global = Y_mm[is_tr].mean(0)
    def climatology(mask):
        return np.vstack([clim_station.loc[s_].values if s_ in clim_station.index else clim_global
                          for s_ in df.loc[mask, "STATION_ID"]])
    no_train = sorted(set(df.loc[is_va | is_te, "STATION_ID"]) - set(clim_station.index))
    print(f"Estaciones de validación/prueba sin años de entrenamiento (usan la media global): {len(no_train)}")

    res_rows, per_curve, inter_rows = [], [], []
    for split, mask, loader in [("validacion", is_va, va_loader), ("prueba", is_te, te_loader)]:
        info = df.loc[mask, ["STATION_ID", "YEAR", "Fase"]].reset_index(drop=True)
        st_ids = df.loc[mask, "STATION_ID"].values
        for ref_name, ref_curves in [("satelite_sin_corregir", X_mm[mask]), ("climatologia_estacion", climatology(mask))]:
            base = curve_metrics(Y_mm[mask], ref_curves)
            res_rows.append({"conjunto": split, "modelo": ref_name, **{k + "_media": v for k, v in base.mean().items()}})
            per_curve.append(pd.concat([info, base], axis=1).assign(conjunto=split, modelo=ref_name))
            inter_rows.append({"conjunto": split, "modelo": ref_name, **interannual_metrics(Y_mm[mask], ref_curves, st_ids)})
        for kind, ms in final.items():
            mets, clipped, inter = [], [], []
            for m in ms:
                p, _ = predict(m, loader)
                clipped.append(100 * (p < 0).mean())
                mets.append(curve_metrics(Y_mm[mask], inverse(SCALE, p)))
                inter.append(interannual_metrics(Y_mm[mask], inverse(SCALE, p), st_ids))
            inter_rows.append({"conjunto": split, "modelo": kind, **pd.DataFrame(inter).mean().to_dict()})
            means = pd.DataFrame([x.mean() for x in mets])
            row = {"conjunto": split, "modelo": kind, "recorte_pct": np.mean(clipped)}
            for k in means:
                row[k + "_media"], row[k + "_DE"] = means[k].mean(), means[k].std(ddof=1)
            res_rows.append(row)
            per_curve.append(pd.concat([info, sum(mets) / len(mets)], axis=1).assign(conjunto=split, modelo=kind))
    res_df = pd.DataFrame(res_rows); per_curve = pd.concat(per_curve, ignore_index=True)
    inter_df = pd.DataFrame(inter_rows); inter_df.to_csv(OUT_DIR / "variabilidad_interanual.csv", index=False)
    print(inter_df.round(3).to_string(index=False))
    res_df.to_csv(OUT_DIR / "resultados_ablacion.csv", index=False)
    per_curve.to_csv(OUT_DIR / "metricas_por_curva.csv", index=False)
    print(res_df[["conjunto", "modelo"] + [c for c in res_df if c.endswith("_media")] + ["recorte_pct"]].round(3).to_string(index=False))

    # ---- 8.9b Intervalos bootstrap por estación (4.4) ----
    # Cada métrica se promedia por estación y se remuestrean ESTACIONES con reemplazo: las curvas de
    # una misma estación no son independientes. (La prueba de Wilcoxon pareada por estación se hace en R.)
    N_BOOT = 2000
    rng = np.random.default_rng(SEED)
    met_cols = ["L2_rel", "KGE", "sesgo_pct", "err_pico1", "err_pico2"]
    by_station = per_curve.groupby(["conjunto", "modelo", "STATION_ID"])[met_cols].mean().reset_index()
    ci_rows = []
    for (split, model_name), g in by_station.groupby(["conjunto", "modelo"]):
        for met in met_cols:
            v = g[met].values
            boots = v[rng.integers(0, len(v), (N_BOOT, len(v)))].mean(1)
            lo, hi = np.percentile(boots, [2.5, 97.5])
            ci_rows.append({"conjunto": split, "modelo": model_name, "metrica": met, "media": v.mean(),
                            "ic_inf": lo, "ic_sup": hi, "estaciones": len(v)})
    ci_df = pd.DataFrame(ci_rows); ci_df.to_csv(OUT_DIR / "ic_bootstrap_por_estacion.csv", index=False)

    # Diferencias pareadas entre modelos: una cabeza supera a otra solo si el IC del 95 % excluye el cero
    comparisons = [("kan", "linear"), ("mlp", "linear"), ("mlp", "kan"),
                   ("linear", "satelite_sin_corregir"), ("kan", "satelite_sin_corregir"), ("mlp", "satelite_sin_corregir"),
                   ("linear", "climatologia_estacion"), ("kan", "climatologia_estacion"), ("mlp", "climatologia_estacion")]
    diff_rows = []
    for split in ["validacion", "prueba"]:
        for met in ["L2_rel", "KGE"]:
            piv = by_station[by_station.conjunto == split].pivot(index="STATION_ID", columns="modelo", values=met)
            for a, b in comparisons:
                d = (piv[a] - piv[b]).values
                boots = d[rng.integers(0, len(d), (N_BOOT, len(d)))].mean(1)
                lo, hi = np.percentile(boots, [2.5, 97.5])
                diff_rows.append({"conjunto": split, "metrica": met, "comparacion": f"{a} − {b}", "dif_media": d.mean(),
                                  "ic_inf": lo, "ic_sup": hi, "excluye_cero": bool(lo > 0 or hi < 0), "estaciones": len(d)})
    diff_df = pd.DataFrame(diff_rows); diff_df.to_csv(OUT_DIR / "diferencias_bootstrap_por_estacion.csv", index=False)
    print(diff_df[diff_df.metrica == "L2_rel"].round(4).to_string(index=False))

    # Resultados por fase ENSO (ONI de la NOAA, clasificación del script de R)
    enso_df = per_curve.groupby(["conjunto", "Fase", "modelo"])[met_cols].mean()
    enso_df["curvas"] = per_curve.groupby(["conjunto", "Fase", "modelo"]).size()
    enso_df.round(4).to_csv(OUT_DIR / "resultados_por_fase_ENSO.csv")
    print(enso_df["L2_rel"].unstack("modelo").round(3))

    # Sesgo del satélite por conjunto: explica la sobrestimación en los años de prueba
    bias_rows = []
    for split, mask in [("entrenamiento", is_tr), ("validacion", is_va), ("prueba", is_te)]:
        est, sat = Y_mm[mask].sum(1).mean(), X_mm[mask].sum(1).mean()
        bias_rows.append({"conjunto": split, "curvas": int(mask.sum()), "ONI_medio": df.loc[mask, "ONI"].mean(),
                          "lluvia_anual_estacion_mm": est, "lluvia_anual_satelite_mm": sat,
                          "sesgo_satelite_pct": 100 * (sat - est) / est})
    bias_df = pd.DataFrame(bias_rows); bias_df.to_csv(OUT_DIR / "sesgo_satelite_por_conjunto.csv", index=False)
    print(bias_df.round(1).to_string(index=False))

    # =====================================================================================
    # 8.9c Interpretación (4.5) con las gráficas de la clase (S14: 3.7, 7.3 y 8.2)
    # =====================================================================================
    MONTHS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    MONTH_LBL = ["E", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
    MONTH_ES = ["ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic"]
    def day_label(d):
        # día del año (0-365) → "dd-mes"
        m_ = max(i for i in range(12) if MONTHS[i] <= d)
        return f"{int(d - MONTHS[m_]) + 1:02d}-{MONTH_ES[m_]}"
    seasons = pd.read_csv(DATA_DIR / "temporadas_lluvia.csv")
    seasons = seasons[seasons.Fuente == "Estaciones"][["t_inicio", "t_fin"]].values   # temporadas del ítem 2
    def month_axis(ax):
        ax.set_xticks(MONTHS); ax.set_xticklabels(MONTH_LBL)
        for a_, b_ in seasons:
            ax.axvspan(a_, b_, color="tab:blue", alpha=0.07)

    m_lin, m_kan, m_mlp = final["linear"][0], final["kan"][0], final["mlp"][0]
    for m in (m_lin, m_kan, m_mlp):
        m.eval()

    # (a) FPCA de la corrección estación − climatología: varianza, autofunciones y curvas medias (S14, 3.7)
    fig, axes = plt.subplots(1, 3, figsize=(17, 3.6))
    axes[0].bar(range(1, 9), 100 * fpca["explained"][:8]); axes[0].plot(range(1, 9), 100 * np.cumsum(fpca["explained"][:8]), "ro-")
    axes[0].set_title("Varianza explicada de la corrección (%)"); axes[0].set_xlabel("componente k")
    axes[1].plot(days, fpca["mean"], "k--", lw=1.5, label="μ_r(t): corrección media")
    for k in range(3):
        axes[1].plot(days, fpca["basis"][k], label=f"φ{k + 1} ({100 * fpca['explained'][k]:.1f}%)")
    axes[1].axhline(0, c="k", lw=0.5); axes[1].set_title("Corrección media y autofunciones (escala √P)"); axes[1].legend(fontsize=7); month_axis(axes[1])
    axes[2].plot(days, Y_mm[is_tr].mean(0), "k", lw=2, label="estaciones")
    axes[2].plot(days, X_mm[is_tr].mean(0), color="gray", ls="--", label="satélite")
    axes[2].plot(days, inverse(SCALE, Xs[is_tr] + fpca["mean"]).mean(0), "r:", label="climatología + corrección media")
    axes[2].set_title("Curvas medias de entrenamiento (mm/día)"); axes[2].legend(fontsize=7); month_axis(axes[2])
    plt.tight_layout(); fig.savefig(FIG_DIR / "Fig_IA_01_FPCA.png"); plt.close(fig)

    # (b) Funciones de pesos w_k(t): qué ventana del año usa cada neurona funcional
    wk_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(16, 3.8))
    for ax, (name, m) in zip(axes, [("lineal", m_lin), ("KAN", m_kan)]):
        A = m.weight_coef.detach().cpu().numpy()
        w = Bx @ A.T                                            # w_k(t) en la malla diaria
        for k in range(A.shape[0]):
            ax.plot(days, w[:, k], label=f"w{k + 1}")
            on = np.abs(w[:, k]) >= 0.5 * np.abs(w[:, k]).max()  # ventana: |w_k| ≥ 50 % de su máximo
            i_pk = np.abs(w[:, k]).argmax()
            wk_rows.append({"cabeza": name, "neurona": k + 1, "pico": day_label(days[i_pk]),
                            "signo_pico": int(np.sign(w[i_pk, k])), "dias_en_ventana": int(on.sum()),
                            "primer_dia_ventana": day_label(days[on].min()), "ultimo_dia_ventana": day_label(days[on].max())})
        ax.axhline(0, c="k", lw=0.5); month_axis(ax)
        ax.set_title(f"Funciones de pesos w_k(t), cabeza {name} (sombreado: temporadas de lluvia)"); ax.legend(fontsize=7, ncol=3)
    plt.tight_layout(); fig.savefig(FIG_DIR / "Fig_IA_02_pesos_wk.png"); plt.close(fig)
    wk_df = pd.DataFrame(wk_rows); wk_df.to_csv(OUT_DIR / "ventanas_wk.csv", index=False)

    # (c) Aristas de la cabeza KAN hacia cada puntaje, en el rango de los datos (S14, 8.2)
    kan_l = m_kan.head[-1]                                      # KANLayer(32 → K)
    with torch.no_grad():
        u_tr = kan_l.norm(m_kan.head[0](m_kan.encode(torch.tensor(X[is_tr]).to(DEVICE)))).cpu().numpy()
    lo_u, hi_u = np.percentile(u_tr, [1, 99])                   # rango de las entradas de la KAN
    g = torch.linspace(float(lo_u), float(hi_u), 200).to(DEVICE)
    amp = np.zeros((K_FPCA, kan_l.in_features)); edge_curves = {}
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.2))
    for k_out in range(K_FPCA):
        c_ = torch.stack([kan_l.edge_function(k_out, i, g) for i in range(kan_l.in_features)]).cpu().numpy()
        amp[k_out] = c_.max(1) - c_.min(1)
        for i in range(kan_l.in_features):
            edge_curves[(k_out, i)] = c_[i]
        if k_out < 3:
            for i in np.argsort(amp[k_out])[::-1][:4]:
                axes[k_out].plot(g.cpu(), c_[i] - c_[i].mean(), label=f"u{i}")
            axes[k_out].set_title(f"Aristas hacia ξ{k_out + 1} (top 4)"); axes[k_out].set_xlabel("entrada latente normalizada"); axes[k_out].legend(fontsize=7)
    plt.tight_layout(); fig.savefig(FIG_DIR / "Fig_IA_03_aristas_KAN.png"); plt.close(fig)
    pd.DataFrame(amp, index=[f"xi{k + 1}" for k in range(K_FPCA)], columns=[f"u{i}" for i in range(kan_l.in_features)]).round(4).to_csv(OUT_DIR / "amplitud_aristas_KAN.csv")

    # (d) Ajuste simbólico de las 3 aristas de mayor amplitud: a + b·g(x) por mínimos cuadrados, con R²
    xg = g.cpu().numpy().astype(np.float64)
    CANDIDATES = {"a + b·x": [xg], "a + b·x²": [xg ** 2], "a + b·x + c·x²": [xg, xg ** 2],
                  "a + b·x + c·x² + d·x³": [xg, xg ** 2, xg ** 3], "a + b·sin(x)": [np.sin(xg)],
                  "a + b·tanh(x)": [np.tanh(xg)], "a + b·exp(x)": [np.exp(xg)], "a + b·|x|": [np.abs(xg)]}
    sym_rows = []
    top_edges = sorted(edge_curves, key=lambda e_: -amp[e_[0], e_[1]])[:3]
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.2))
    for ax, (k_out, i) in zip(axes, top_edges):
        yv = edge_curves[(k_out, i)].astype(np.float64)
        fits = []
        for formula, cols in CANDIDATES.items():
            Fm = np.column_stack([np.ones_like(xg)] + cols)
            coef = np.linalg.solve(Fm.T @ Fm, Fm.T @ yv)
            fits.append((formula, coef, regression_metrics(yv, Fm @ coef)["R2"], len(cols), Fm @ coef))
        # la expresión más simple con R² ≥ 0.99; si ninguna llega, la de mayor R²
        good = [f_ for f_ in fits if f_[2] >= 0.99]
        best_fit = min(good, key=lambda f_: f_[3]) if good else max(fits, key=lambda f_: f_[2])
        for formula, coef, r2, n_terms, _ in fits:
            sym_rows.append({"arista": f"u{i} → ξ{k_out + 1}", "amplitud": amp[k_out, i], "expresion": formula,
                             "coeficientes": np.round(coef, 4).tolist(), "R2": r2, "elegida": formula == best_fit[0]})
        ax.plot(xg, yv, "k", lw=2, label="φ aprendida")
        ax.plot(xg, best_fit[4], "r--", label=f"{best_fit[0]} (R²={best_fit[2]:.3f})")
        ax.set_title(f"Arista u{i} → ξ{k_out + 1}"); ax.set_xlabel("entrada latente normalizada"); ax.legend(fontsize=7)
    plt.tight_layout(); fig.savefig(FIG_DIR / "Fig_IA_04_ajuste_simbolico.png"); plt.close(fig)
    sym_df = pd.DataFrame(sym_rows); sym_df.to_csv(OUT_DIR / "ajuste_simbolico.csv", index=False)
    print(sym_df[sym_df.elegida][["arista", "expresion", "coeficientes", "R2"]].to_string(index=False))

    # (e) Contraste con la cabeza lineal: efecto parcial del ONI sobre la lluvia anual
    #     (misma entrada en las tres cabezas; se varía el ONI y se promedia sobre las curvas de
    #      entrenamiento y sobre las 3 semillas de cada cabeza)
    oni_tr = df.loc[is_tr, "ONI"].values.astype(float)
    oni_mu, oni_sd = oni_tr.mean(), oni_tr.std()
    grid_v = np.linspace(*np.percentile(oni_tr, [1, 99]), 25)
    pd_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(13, 3.8))
    for name, kind in [("lineal", "linear"), ("KAN", "kan"), ("MLP", "mlp")]:
        tot = []
        for v in grid_v:
            Xv = X[is_tr].copy(); Xv[:, K_X] = (v - oni_mu) / oni_sd
            per_seed = [inverse(SCALE, predict(m, make_loader(Xv, Y[is_tr], 512))[0]).sum(1).mean() for m in final[kind]]
            tot.append(np.mean(per_seed))
            pd_rows.append({"ONI": v, "cabeza": name, "lluvia_anual_mm": tot[-1], "DE_semillas": np.std(per_seed, ddof=1)})
        axes[0].plot(grid_v, tot, marker="o", ms=3, label=name)
    axes[0].set_xlabel("ONI"); axes[0].set_ylabel("lluvia anual media (mm)")
    axes[0].set_title("Efecto parcial del ONI (promedio de 3 semillas)"); axes[0].legend(fontsize=8)
    # Coeficientes de la cabeza lineal (semilla 0): cuánto mueve cada entrada a cada puntaje
    Wlin = m_lin.head.weight.detach().cpu().numpy()
    im = axes[1].imshow(Wlin, cmap="RdBu_r", vmin=-np.abs(Wlin).max(), vmax=np.abs(Wlin).max())
    axes[1].set_xticks(range(Wlin.shape[1])); axes[1].set_xticklabels([f"z{k + 1}" for k in range(N_NEURONS)] + ["ONI"])
    axes[1].set_yticks(range(K_FPCA)); axes[1].set_yticklabels([f"ξ{k + 1}" for k in range(K_FPCA)])
    axes[1].set_title("Coeficientes de la cabeza lineal (semilla 0)"); plt.colorbar(im, ax=axes[1])
    plt.tight_layout(); fig.savefig(FIG_DIR / "Fig_IA_05_contraste_lineal.png"); plt.close(fig)
    pd.DataFrame(pd_rows).to_csv(OUT_DIR / "efecto_parcial_ONI.csv", index=False)
    pd.DataFrame(Wlin, index=[f"xi{k + 1}" for k in range(K_FPCA)],
                 columns=[f"z{k + 1}" for k in range(N_NEURONS)] + ["ONI"]).round(4).to_csv(OUT_DIR / "coeficientes_cabeza_lineal.csv")

    # (f) Puntajes reales contra predichos en prueba (S14, 8.2)
    true_scores = (Y[is_te] - Xs[is_te] - fpca["mean"]) @ fpca["basis"].T   # puntajes de la corrección real
    with torch.no_grad():
        _, pred_scores = m_kan(torch.tensor(X[is_te]).to(DEVICE), return_scores=True)
    pred_scores = pred_scores.cpu().numpy()
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.6))
    for k, ax in enumerate(axes):
        r = np.corrcoef(true_scores[:, k], pred_scores[:, k])[0, 1]
        ax.scatter(true_scores[:, k], pred_scores[:, k], s=4, alpha=0.4)
        lim = np.percentile(true_scores[:, k], [1, 99]); ax.plot(lim, lim, "k--")
        ax.set_title(f"ξ{k + 1}: real vs predicho, KAN (r = {r:.2f})"); ax.set_xlabel("puntaje real"); ax.set_ylabel("puntaje predicho")
    plt.tight_layout(); fig.savefig(FIG_DIR / "Fig_IA_06_puntajes.png"); plt.close(fig)

    # (g) Peores casos en prueba según la KAN (promedio de las 3 semillas), con sus curvas
    pc_te = per_curve[(per_curve.conjunto == "prueba")]
    kan_te = pc_te[pc_te.modelo == "kan"].reset_index(drop=True)
    worst = kan_te.L2_rel.sort_values(ascending=False).index[:4]
    # Curvas predichas: promedio de las 3 semillas de cada cabeza
    preds_te = {name: np.mean([inverse(SCALE, predict(m, te_loader)[0]) for m in final[kind]], axis=0)
                for name, kind in [("lineal", "linear"), ("KAN", "kan"), ("MLP", "mlp")]}
    clim_te = climatology(is_te)
    Yte_mm, Xte_mm = Y_mm[is_te], X_mm[is_te]
    worst_rows = []
    fig, axes = plt.subplots(1, 4, figsize=(19, 3.6))
    for ax, i in zip(axes, worst):
        info_i = kan_te.loc[i]
        ax.plot(days, Yte_mm[i], "k", lw=2, label="estación")
        ax.plot(days, Xte_mm[i], color="gray", ls="--", label="satélite")
        ax.plot(days, clim_te[i], color="purple", ls=":", label="climatología")
        for name, c in [("lineal", "tab:blue"), ("KAN", "tab:orange"), ("MLP", "tab:green")]:
            ax.plot(days, preds_te[name][i], color=c, lw=1, label=name)
        month_axis(ax); ax.set_title(f"{info_i.STATION_ID} · {info_i.YEAR} · {info_i.Fase} · L²={info_i.L2_rel:.2f}", fontsize=9)
        ax.set_ylabel("mm/día"); ax.legend(fontsize=6)
        worst_rows.append({"STATION_ID": info_i.STATION_ID, "YEAR": info_i.YEAR, "Fase": info_i.Fase,
                           "L2_rel_KAN": info_i.L2_rel, "sesgo_pct_KAN": info_i.sesgo_pct,
                           "lluvia_anual_estacion_mm": Yte_mm[i].sum(), "lluvia_anual_satelite_mm": Xte_mm[i].sum()})
    plt.tight_layout(); fig.savefig(FIG_DIR / "Fig_IA_07_peores_casos.png"); plt.close(fig)
    pd.DataFrame(worst_rows).to_csv(OUT_DIR / "peores_casos.csv", index=False)

    # ---- 8.10 Predicciones para R (formato station_id, year, day, value_mm) y checkpoints ----
    pred_dir = OUT_DIR / "predicciones"; pred_dir.mkdir(exist_ok=True)
    mask = is_va | is_te
    def write_pred(curves, name):
        info = df.loc[mask]
        pd.DataFrame({"station_id": np.repeat(info.STATION_ID.values, 365), "year": np.repeat(info.YEAR.values, 365),
                      "day": np.tile(np.arange(1, 366), mask.sum()), "value_mm": np.round(curves.ravel(), 4)}).to_csv(pred_dir / name, index=False)
    write_pred(Y_mm[mask], "obs_estaciones_mm.csv")
    write_pred(X_mm[mask], "pred_satelite_sin_corregir.csv")
    all_loader = make_loader(X[mask], Y[mask], 512)
    for kind, ms in final.items():
        for s, m in zip(SEEDS, ms):
            p, _ = predict(m, all_loader)
            write_pred(inverse(SCALE, p), f"pred_{kind}_s{s}.csv")
            save_checkpoint(m, CKPT_DIR / f"{kind}_s{s}.pth", meta={"escala": SCALE, "K": K_FPCA, "seed": s})

    # ---- 8.11 Resumen de verificaciones ----
    checks = pd.DataFrame([
        {"verificacion": "Base B-spline de Python contra R", "valor": diff_control, "umbral": "≈ 0"},
        {"verificacion": "Gram satélite 1ᵀW1 (longitud del dominio)", "valor": W_x.sum(), "umbral": 365},
        {"verificacion": "Pérdida en la malla − (ĉ−c)ᵀW(ĉ−c)", "valor": abs(loss_grid - loss_W), "umbral": "≈ 0"},
        {"verificacion": "KAN sintética: RMSE de validación", "valor": rmse_toy, "umbral": "≤ 1e-2"},
        {"verificacion": "Equivalencia lineal: diferencia máxima", "valor": max(eq_diffs), "umbral": "≤ 1e-4"},
        {"verificacion": "Parámetros MLP vs KAN (%)", "valor": param_df.loc["mlp", "dif_vs_KAN_pct"], "umbral": "≤ 10"},
        {"verificacion": "Neuronas funcionales elegidas", "valor": N_NEURONS, "umbral": "6 a 8"},
        {"verificacion": "Remuestreos bootstrap por estación", "valor": N_BOOT, "umbral": "≥ 1000"},
    ])
    checks.to_csv(OUT_DIR / "verificaciones.csv", index=False)
    print(checks.to_string(index=False))
    print(f"Tiempo total: {time.time() - t_start:.0f} s")
