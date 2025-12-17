#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helmholtz 2D em (0,1)x(0,1) com Dirichlet imposto pela solução exata (3 ondas planas)
Método: Diferenças finitas centradas de 2ª ordem (stencil 5 pontos) em malha uniforme.

Este script roda PARA TODOS OS GRUPOS (1..6) por padrão, gerando:
- solução numérica (complexa)
- erro relativo L2 discreto
- gráficos: superfície 3D + cortes em x=0.5 e y=0.5 (numérico vs exato)
- métricas em CSV por grupo e um CSV consolidado

Malhas (N) ajustadas para capturar bem k=100:
- Por padrão: Ns = [64, 128, 192, 256]
  (para k=100, N=256 dá ~16 pontos por comprimento de onda)

Uso:
  # Script conveniente (recomendado)
  ./run.sh  # usa config.yaml
  ./run.sh --config meu_config.yaml
  ./run.sh --groups 1 3 6 --Ns 64 128 192 256 --ks 1 20 40 100
  ./run.sh --solver auto --plot_what real

  # Ou ativar ambiente virtual manualmente
  source .venv/bin/activate && python main.py

  # Ou usar python3 diretamente (se PyYAML estiver instalado globalmente)
  python3 main.py

Saídas:
  outputs/metrics_all_groups.csv
  outputs/metrics_group1.csv ... metrics_group6.csv
  figures/group_1/*.pdf ... figures/group_6/*.pdf
"""

import argparse
import csv
import os
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from functools import lru_cache

import numpy as np
import matplotlib.pyplot as plt
import yaml

from scipy.sparse import csc_matrix, diags, eye, kron
from scipy.sparse.linalg import splu, spsolve, gmres, spilu, LinearOperator




@dataclass
class CaseResult:
    group: int
    N: int
    k: float
    h: float
    err_rel_l2: float
    t_build_s: float
    t_solve_s: float
    solver: str
    gmres_info: int


def ensure_dirs(config: Dict) -> None:
    os.makedirs("figures", exist_ok=True)
    os.makedirs("outputs", exist_ok=True)
    group_thetas = config.get('group_thetas', {})
    for g in group_thetas.keys():
        os.makedirs(f"figures/group_{g}", exist_ok=True)


def u_exact(x: np.ndarray, y: np.ndarray, k: float, thetas: Tuple[float, float, float]) -> np.ndarray:
    """u_ex(x,y) = sum exp(i*k*(cos(theta_i)*x + sin(theta_i)*y))"""
    val = np.zeros(np.broadcast(x, y).shape, dtype=np.complex128)
    for th in thetas:
        val += np.exp(1j * k * (np.cos(th) * x + np.sin(th) * y))
    return val


# Cache global para matrizes Laplacianas (reutilização por N)
_laplacian_cache: Dict[int, Tuple[csc_matrix, float]] = {}


def build_laplacian_interior(N: int, use_cache: bool = True) -> Tuple[csc_matrix, float]:
    """
    Convenção:
      - N subintervalos => N+1 pontos
      - h = 1/N
      - internos: 1..N-1 (total (N-1)^2 incógnitas)
    
    Com cache: reutiliza matrizes para mesmo N (otimização).
    """
    if use_cache and N in _laplacian_cache:
        return _laplacian_cache[N]
    
    if N < 2:
        raise ValueError("N deve ser >= 2.")
    n1 = N - 1
    h = 1.0 / N

    main = -2.0 * np.ones(n1)
    off = 1.0 * np.ones(n1 - 1)
    T = diags([off, main, off], offsets=[-1, 0, 1], format="csc")
    I = eye(n1, format="csc")
    L = (kron(T, I, format="csc") + kron(I, T, format="csc")) / (h * h)
    
    result = (L, h)
    if use_cache:
        _laplacian_cache[N] = result
    return result


def build_rhs_from_dirichlet(N: int, h: float, g: np.ndarray) -> np.ndarray:
    """
    RHS para Δu + k^2 u = 0 com Dirichlet via contorno conhecido:
      vizinho de contorno entra como -(1/h^2)*g_vizinho no RHS.
    Ordenação do vetor: order='F' (compatível com kron(T,I)+kron(I,T)).
    """
    n1 = N - 1
    inv_h2 = 1.0 / (h * h)
    b = np.zeros((n1, n1), dtype=np.complex128)

    b[0, :]        -= inv_h2 * g[0, 1:N]   # x=0
    b[n1 - 1, :]   -= inv_h2 * g[N, 1:N]   # x=1
    b[:, 0]        -= inv_h2 * g[1:N, 0]   # y=0
    b[:, n1 - 1]   -= inv_h2 * g[1:N, N]   # y=1

    return b.reshape((n1 * n1,), order="F")


def solve_linear_system(A: csc_matrix, b: np.ndarray, solver: str) -> Tuple[np.ndarray, int, str]:
    """
    Retorna (x, gmres_info, solver_used).
    solver:
      - "splu": LU
      - "spsolve": solve direto
      - "gmres_ilu": GMRES com ILU
      - "auto": splu para sistemas menores; gmres_ilu para maiores
    """
    n = A.shape[0]

    if solver == "auto":
        # Heurística conservadora: LU até ~40k incógnitas; acima, GMRES+ILU.
        solver = "splu" if n <= 40000 else "gmres_ilu"

    if solver == "splu":
        lu = splu(A)
        x = lu.solve(b)
        return x, 0, "splu"

    if solver == "spsolve":
        x = spsolve(A, b)
        return x, 0, "spsolve"

    if solver == "gmres_ilu":
        # ILU como pré-condicionador
        A = A.astype(np.complex128)

        try:
            # Parâmetros ILU otimizados para melhor desempenho
            ilu = spilu(A, drop_tol=1e-3, fill_factor=20)
            Mx = lambda v: ilu.solve(v)
            M = LinearOperator(A.shape, Mx, dtype=np.complex128)

            # Tolerâncias mais razoáveis para números complexos
            x, info = gmres(A, b, M=M, rtol=1e-8, atol=1e-12, restart=100, maxiter=1000)
            return x, int(info), "gmres_ilu"
        except RuntimeError as e:
                if "singular" in str(e).lower():
                    # Se ILU falha por singularidade, tentar GMRES sem pré-condicionador
                    # ou fazer fallback para SPLU se o sistema não for muito grande
                    if n <= 100000:  # Limite razoável para SPLU
                        try:
                            lu = splu(A)
                            x = lu.solve(b)
                            return x, 0, "splu_fallback"
                        except (RuntimeError, ValueError):
                            pass
                    
                    # Último recurso: GMRES sem pré-condicionador
                    x, info = gmres(A, b, rtol=1e-8, atol=1e-12, restart=100, maxiter=1000)
                    return x, int(info), "gmres_noprec"
                else:
                    raise

    raise ValueError("solver inválido. Use: splu | spsolve | gmres_ilu | auto")


# Cache global para coordenadas (reutilização por N)
_coords_cache: Dict[int, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}


def get_coordinates(N: int, use_cache: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Retorna (x, y, X, Y) para dado N, com cache."""
    if use_cache and N in _coords_cache:
        return _coords_cache[N]
    
    x = np.linspace(0.0, 1.0, N + 1)
    y = np.linspace(0.0, 1.0, N + 1)
    X, Y = np.meshgrid(x, y, indexing="ij")
    
    result = (x, y, X, Y)
    if use_cache:
        _coords_cache[N] = result
    return result


def solve_case(group: int, N: int, k: float, solver: str, config: Dict, use_cache: bool = True) -> Tuple[np.ndarray, np.ndarray, CaseResult]:
    """
    Resolve caso (grupo, N, k), retorna:
      - u_num_full (N+1 x N+1)
      - u_ex_full  (N+1 x N+1)
      - métricas
    
    Com cache: reutiliza matrizes Laplacianas e coordenadas para mesmo N.
    """
    group_thetas = config.get('group_thetas', {})
    thetas = group_thetas[group]

    t0 = time.perf_counter()
    L, h = build_laplacian_interior(N, use_cache=use_cache)

    x, y, X, Y = get_coordinates(N, use_cache=use_cache)

    u_ex_full = u_exact(X, Y, k, thetas)

    # contorno Dirichlet
    g = np.zeros_like(u_ex_full)
    g[0, :] = u_ex_full[0, :]
    g[N, :] = u_ex_full[N, :]
    g[:, 0] = u_ex_full[:, 0]
    g[:, N] = u_ex_full[:, N]

    b = build_rhs_from_dirichlet(N, h, g)

    n = (N - 1) * (N - 1)
    #A = (L + (k * k) * eye(n, format="csc")).tocsc()

    A = (L.astype(np.complex128) + (k * k) * eye(n, format="csc", dtype=np.complex128)).tocsc()

    t_build = time.perf_counter() - t0

    t1 = time.perf_counter()
    u_int, info, solver_used = solve_linear_system(A, b, solver=solver)
    t_solve = time.perf_counter() - t1

    if solver_used.startswith("gmres") and info != 0:
        raise RuntimeError(f"GMRES não convergiu (info={info}) para group={group}, N={N}, k={k}.")

    # solução completa
    u_num_full = np.zeros((N + 1, N + 1), dtype=np.complex128)
    u_num_full[0, :] = g[0, :]
    u_num_full[N, :] = g[N, :]
    u_num_full[:, 0] = g[:, 0]
    u_num_full[:, N] = g[:, N]

    u_int_mat = np.reshape(u_int, (N - 1, N - 1), order="F")
    u_num_full[1:N, 1:N] = u_int_mat

    # erro relativo L2 discreto com peso h^2
    diff = u_num_full - u_ex_full
    num = np.sqrt(np.sum(np.abs(diff) ** 2) * (h * h))
    den = np.sqrt(np.sum(np.abs(u_ex_full) ** 2) * (h * h))
    err_rel = float(num / den) if den != 0 else np.nan

    res = CaseResult(
        group=group,
        N=N,
        k=float(k),
        h=float(h),
        err_rel_l2=err_rel,
        t_build_s=float(t_build),
        t_solve_s=float(t_solve),
        solver=solver_used,
        gmres_info=int(info),
    )
    return u_num_full, u_ex_full, res


def save_metrics(rows: List[CaseResult], path: str) -> None:
    header = ["group", "N", "k", "h", "err_rel_l2", "t_build_s", "t_solve_s", "solver", "gmres_info"]
    write_header = not os.path.exists(path)
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if write_header:
            w.writerow(header)
        for r in rows:
            w.writerow([r.group, r.N, r.k, r.h, r.err_rel_l2, r.t_build_s, r.t_solve_s, r.solver, r.gmres_info])


def plot_surface_3d(U: np.ndarray, group: int, N: int, k: float, what: str) -> str:
    x = np.linspace(0.0, 1.0, N + 1)
    y = np.linspace(0.0, 1.0, N + 1)
    X, Y = np.meshgrid(x, y, indexing="ij")

    if what == "real":
        Z = np.real(U)
        zlabel = "Re(u)"
        title = f"Grupo {group} | Superfície Re(u) | k={k}, N={N}"
    elif what == "abs":
        Z = np.abs(U)
        zlabel = "|u|"
        title = f"Grupo {group} | Superfície |u| | k={k}, N={N}"
    else:
        raise ValueError("plot_what deve ser 'real' ou 'abs'.")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z, linewidth=0, antialiased=True)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel(zlabel)
    ax.set_title(title)

    out = f"figures/group_{group}/surf_k{int(k)}_N{N}_{what}.pdf"
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def calculate_convergence_rate(results: List[CaseResult], group: int, k: float) -> Optional[float]:
    """
    Calcula taxa de convergência para um grupo e k específicos.
    Taxa p: erro(h) ≈ C * h^p
    """
    # Filtrar resultados para grupo e k específicos, ordenar por N
    filtered = [r for r in results if r.group == group and abs(r.k - k) < 1e-6]
    filtered.sort(key=lambda x: x.N)
    
    if len(filtered) < 2:
        return None
    
    # Calcular taxa usando pares consecutivos
    # Taxa p: erro(h) ≈ C * h^p, então p = log(erro1/erro2) / log(h1/h2)
    rates = []
    for i in range(len(filtered) - 1):
        r1, r2 = filtered[i], filtered[i + 1]
        if r1.err_rel_l2 > 0 and r2.err_rel_l2 > 0 and r1.h > r2.h:
            p = np.log(r1.err_rel_l2 / r2.err_rel_l2) / np.log(r1.h / r2.h)
            rates.append(p)
    
    return np.mean(rates) if rates else None


def plot_convergence(results: List[CaseResult], group: int, k: float, output_dir: str = "figures") -> Optional[str]:
    """
    Gera gráfico de convergência (erro vs h) em escala log-log.
    """
    # Filtrar e ordenar
    filtered = [r for r in results if r.group == group and abs(r.k - k) < 1e-6]
    filtered.sort(key=lambda x: x.N)
    
    if len(filtered) < 2:
        return None
    
    h_vals = [r.h for r in filtered]
    err_vals = [r.err_rel_l2 for r in filtered]
    
    # Calcular taxa de convergência
    rate = calculate_convergence_rate(results, group, k)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog(h_vals, err_vals, 'o-', linewidth=2, markersize=8, label=f'Dados (taxa ≈ {rate:.2f})' if rate else 'Dados')
    
    # Linha de referência para ordem 2
    h_ref = np.array(h_vals)
    err_ref = h_ref**2 * (err_vals[0] / h_vals[0]**2)
    ax.loglog(h_ref, err_ref, '--', alpha=0.5, label='Ordem 2 (referência)')
    
    ax.set_xlabel('$h$ (tamanho da malha)', fontsize=12)
    ax.set_ylabel('Erro relativo L2', fontsize=12)
    ax.set_title(f'Convergência: Grupo {group}, $k={int(k)}$', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    os.makedirs(f"{output_dir}/group_{group}", exist_ok=True)
    out = f"{output_dir}/group_{group}/convergence_k{int(k)}.pdf"
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_cuts(U_num: np.ndarray, U_ex: np.ndarray, group: int, N: int, k: float, what: str) -> Tuple[str, str]:
    x = np.linspace(0.0, 1.0, N + 1)
    y = np.linspace(0.0, 1.0, N + 1)

    i0 = int(np.argmin(np.abs(x - 0.5)))
    j0 = int(np.argmin(np.abs(y - 0.5)))

    if what == "real":
        f = np.real
        ylabel = "Re(u)"
    elif what == "abs":
        f = np.abs
        ylabel = "|u|"
    else:
        raise ValueError("plot_what deve ser 'real' ou 'abs'.")

    # corte x ~ 0.5
    fig1 = plt.figure()
    plt.plot(y, f(U_num[i0, :]), label="Numérico")
    plt.plot(y, f(U_ex[i0, :]), label="Exato", linestyle="--")
    plt.xlabel("y")
    plt.ylabel(ylabel)
    plt.title(f"Grupo {group} | Corte em x={x[i0]:.4f} (~0.5) | k={k}, N={N}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    out_x = f"figures/group_{group}/cut_x05_k{int(k)}_N{N}_{what}.pdf"
    fig1.tight_layout()
    fig1.savefig(out_x)
    plt.close(fig1)

    # corte y ~ 0.5
    fig2 = plt.figure()
    plt.plot(x, f(U_num[:, j0]), label="Numérico")
    plt.plot(x, f(U_ex[:, j0]), label="Exato", linestyle="--")
    plt.xlabel("x")
    plt.ylabel(ylabel)
    plt.title(f"Grupo {group} | Corte em y={y[j0]:.4f} (~0.5) | k={k}, N={N}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    out_y = f"figures/group_{group}/cut_y05_k{int(k)}_N{N}_{what}.pdf"
    fig2.tight_layout()
    fig2.savefig(out_y)
    plt.close(fig2)

    return out_x, out_y


def load_config(config_path: str) -> Dict:
    """Carrega configuração do arquivo YAML."""
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        if config and 'group_thetas' in config:
            # Converter listas para tuplas para compatibilidade com código existente
            config['group_thetas'] = {int(k): tuple(v) for k, v in config['group_thetas'].items()}

        return config if config else {}
    except FileNotFoundError:
        print(f"Aviso: Arquivo de configuração '{config_path}' não encontrado. Usando valores padrão.")
        return {}
    except yaml.YAMLError as e:
        raise ValueError(f"Erro ao ler arquivo de configuração '{config_path}': {e}")


def parse_args() -> argparse.Namespace:
    # Primeiro, parse apenas o config para saber qual arquivo usar
    temp_parser = argparse.ArgumentParser(add_help=False)
    temp_parser.add_argument("--config", default="config.yaml")
    temp_args, remaining_argv = temp_parser.parse_known_args()

    # Carrega configuração do arquivo
    config = load_config(temp_args.config)

    # Cria parser completo com defaults da configuração
    p = argparse.ArgumentParser(description="Helmholtz 2D FD (rodar todos os grupos)")

    p.add_argument(
        "--config",
        type=str,
        default=temp_args.config,
        help="Arquivo de configuração YAML (default: config.yaml).",
    )
    p.add_argument(
        "--groups",
        type=int,
        nargs="+",
        default=config.get("groups", [1, 2, 3, 4, 5, 6]),
        help="Grupos a rodar (default do config ou fallback: %(default)s).",
    )
    p.add_argument(
        "--Ns",
        type=int,
        nargs="+",
        default=config.get("Ns", [64, 128, 192, 256]),
        help="Malhas: N subintervalos (pontos = N+1). Default do config ou fallback: %(default)s.",
    )
    p.add_argument(
        "--ks",
        type=float,
        nargs="+",
        default=config.get("ks", [1, 20, 40, 100]),
        help="Valores de k (default do config ou fallback: %(default)s).",
    )
    p.add_argument(
        "--solver",
        type=str,
        default=config.get("solver", "auto"),
        choices=["auto", "splu", "spsolve", "gmres_ilu"],
        help="Solver linear (default do config ou fallback: %(default)s).",
    )
    p.add_argument(
        "--plot_what",
        type=str,
        default=config.get("plot_what", "real"),
        choices=["real", "abs"],
        help="O que plotar: real (Re(u)) ou abs (|u|). Default do config ou fallback: %(default)s.",
    )
    p.add_argument(
        "--no_plots",
        action="store_true",
        default=config.get("no_plots", False),
        help="Não gerar figuras (somente CSV). Default do config ou fallback: %(default)s.",
    )
    p.add_argument(
        "--parallel",
        type=int,
        default=config.get("parallel", 1),
        help="Número de processos paralelos (1=sequencial). Default do config ou fallback: %(default)s.",
    )
    p.add_argument(
        "--no_cache",
        action="store_true",
        default=config.get("no_cache", False),
        help="Desabilitar cache de matrizes e coordenadas.",
    )

    # Parse com todos os argumentos restantes
    return p.parse_args(remaining_argv)


def _process_case_wrapper(args_tuple):
    """Wrapper para multiprocessing (deve estar no nível do módulo)."""
    g, N, k, solver, config_dict, use_cache, no_plots = args_tuple
    u_num, u_ex, res = solve_case(group=g, N=N, k=float(k), solver=solver, config=config_dict, use_cache=use_cache)
    
    if not no_plots:
        return (g, N, k, res, u_num, u_ex)
    return (g, N, k, res, None, None)


def main() -> None:
    args = parse_args()

    # Carregar configuração completa (incluindo group_thetas)
    config = load_config(args.config)
    ensure_dirs(config)

    groups = args.groups
    Ns = args.Ns
    ks = args.ks
    group_thetas = config.get('group_thetas', {})
    use_cache = not args.no_cache
    parallel_workers = max(1, args.parallel)

    # validações simples
    for g in groups:
        if g not in group_thetas:
            available_groups = list(group_thetas.keys())
            raise ValueError(f"Grupo inválido: {g}. Grupos disponíveis: {available_groups}.")
    for N in Ns:
        if N < 2:
            raise ValueError("Todos os N devem ser >= 2.")

    all_rows: List[CaseResult] = []
    use_parallel = parallel_workers > 1

    print("=== Helmholtz 2D (FD 2ª ordem) | Todos os grupos ===")
    print(f"Configuração: {args.config}")
    print(f"Grupos: {groups}")
    print(f"Ns: {Ns}")
    print(f"ks: {ks}")
    print(f"Solver: {args.solver}")
    print(f"Plot: {args.plot_what}")
    print(f"Cache: {'habilitado' if use_cache else 'desabilitado'}")
    if use_parallel:
        print(f"Paralelização: {parallel_workers} processos")
    print()

    # Pré-aquecer cache (construir matrizes e coordenadas uma vez)
    if use_cache:
        print("Pré-aquecendo cache de matrizes...")
        for N in Ns:
            build_laplacian_interior(N, use_cache=True)
            get_coordinates(N, use_cache=True)
        print("Cache aquecido.\n")

    if use_parallel:
        from multiprocessing import Pool
        
        # Preparar todos os casos como tuplas
        cases = [(g, N, k, args.solver, config, use_cache, args.no_plots) 
                 for g in groups for N in Ns for k in ks]
        
        print(f"Processando {len(cases)} casos em paralelo ({parallel_workers} workers)...\n")
        
        with Pool(parallel_workers) as pool:
            results = pool.map(_process_case_wrapper, cases)
        
        # Organizar resultados por grupo (manter ordem original)
        results_dict = {(g, N, k): (res, u_num, u_ex) for g, N, k, res, u_num, u_ex in results}
        
        for g in groups:
            thetas = group_thetas[g]
            print(f"\n--- Grupo {g} | thetas={thetas} ---")
            group_rows: List[CaseResult] = []
            
            for N in Ns:
                for k in ks:
                    res, u_num, u_ex = results_dict[(g, N, k)]
                    
                    group_rows.append(res)
                    all_rows.append(res)
                    
                    print(
                        f"  h={res.h:.6e} | err_rel_L2={res.err_rel_l2:.6e} "
                        f"| build={res.t_build_s:.3f}s | solve={res.t_solve_s:.3f}s | solver={res.solver}"
                    )
                    
                    if not args.no_plots and u_num is not None:
                        plot_surface_3d(u_num, group=g, N=N, k=k, what=args.plot_what)
                        plot_cuts(u_num, u_ex, group=g, N=N, k=k, what=args.plot_what)
            
            # salva CSV do grupo
            save_metrics(group_rows, path=f"outputs/metrics_group{g}.csv")
    else:
        # Execução sequencial (original, otimizada com cache)
        for g in groups:
            thetas = group_thetas[g]
            print(f"\n--- Grupo {g} | thetas={thetas} ---")
            group_rows: List[CaseResult] = []

            for N in Ns:
                for k in ks:
                    print(f"Rodando group={g}, N={N}, k={k} ...")
                    u_num, u_ex, res = solve_case(group=g, N=N, k=float(k), solver=args.solver, config=config, use_cache=use_cache)
                    group_rows.append(res)
                    all_rows.append(res)

                    print(
                        f"  h={res.h:.6e} | err_rel_L2={res.err_rel_l2:.6e} "
                        f"| build={res.t_build_s:.3f}s | solve={res.t_solve_s:.3f}s | solver={res.solver}"
                    )

                    if not args.no_plots:
                        plot_surface_3d(u_num, group=g, N=N, k=k, what=args.plot_what)
                        plot_cuts(u_num, u_ex, group=g, N=N, k=k, what=args.plot_what)

            # salva CSV do grupo
            save_metrics(group_rows, path=f"outputs/metrics_group{g}.csv")

    # salva CSV consolidado
    save_metrics(all_rows, path="outputs/metrics_all_groups.csv")

    # Gerar gráficos de convergência e análise de taxa
    if len(all_rows) > 0 and not args.no_plots:
        print("\n=== Gerando gráficos de convergência ===")
        convergence_plots = []
        for g in groups:
            for k in ks:
                plot_path = plot_convergence(all_rows, group=g, k=float(k))
                if plot_path:
                    convergence_plots.append(plot_path)
        
        if convergence_plots:
            print(f"Grafos de convergência gerados: {len(convergence_plots)}")
        
        # Calcular e imprimir taxas de convergência
        print("\n=== Taxas de Convergência ===")
        print(f"{'Grupo':>6} {'k':>6} {'Taxa':>10} {'Esperado':>10}")
        for g in groups:
            for k in ks:
                rate = calculate_convergence_rate(all_rows, group=g, k=float(k))
                if rate is not None:
                    print(f"{g:6d} {int(k):6d} {rate:10.2f} {2.0:10.2f}")

    # resumo
    print("\n=== Resumo (consolidado) ===")
    print(f"{'grp':>4} {'N':>6} {'k':>6} {'h':>12} {'err_rel_L2':>16} {'build(s)':>10} {'solve(s)':>10} {'solver':>10}")
    for r in all_rows:
        print(f"{r.group:4d} {r.N:6d} {int(r.k):6d} {r.h:12.4e} {r.err_rel_l2:16.6e} {r.t_build_s:10.3f} {r.t_solve_s:10.3f} {r.solver:>10}")

    print("\nArquivos gerados:")
    print(" - outputs/metrics_all_groups.csv")
    print(" - outputs/metrics_group1.csv ... outputs/metrics_group6.csv")
    if not args.no_plots:
        print(" - figures/group_1/*.pdf ... figures/group_6/*.pdf")
        print("   (incluindo gráficos de convergência)")


if __name__ == "__main__":
    main()
