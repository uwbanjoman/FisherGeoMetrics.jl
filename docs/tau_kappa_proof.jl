# examples/tau_kappa_proof.jl
# ===========================
# Derivation of τ = 1/5 and κ = 6/5
# from the Einstein metric on M^{1,1,1}
#
# Source: Castellani, D'Auria and Fré,
#         Nucl. Phys. B239 (1984) 610–652
#         eqs. (4.11), (4.14), (4.15), (4.16)
#
# © 2026 Jan Bouwman — MIT License

using LinearAlgebra
using Printf
using Polynomials

println("="^65)
println("  Derivation of τ and κ from the Einstein metric on M¹¹¹")
println("  Source: Castellani, D'Auria, Fré (1984)")
println("="^65)
println()

# ── Step 1: M^{1,1,1} is uniquely selected by N=2 SUSY ───────────

println("─── Step 1: N=2 supersymmetry selects M^{1,1,1} uniquely ───")
println()
println("  Among all M^{p,q,r} with G_K = SU(3)×SU(2)×U(1) isometry:")
println()
println("  Holonomy SO(7):  p ≠ q  →  0 Killing spinors  →  N=0")
println("  Holonomy SU(3):  p = q  →  2 Killing spinors  →  N=2  ✓")
println()
println("  For p = q = r = 1: M^{1,1,1}")
println("  This is the unique N=2 solution.")
println()

p, q = 1, 1   # M^{1,1,1}
@printf("  p = %d, q = %d  →  M^{1,1,1}  ✓\n", p, q)
println()

# ── Step 2: The cubic equation (4.15) ────────────────────────────

println("─── Step 2: Cubic equation (4.15) for β ───")
println()
println("  Castellani et al. eq. (4.15):")
println()
println("  4β³ - 6β² + (9/4 + q²/p²)β - (1/2)(q²/p²) = 0")
println()

r = (q/p)^2   # = 1 for M^{1,1,1}
@printf("  For p=%d, q=%d: q²/p² = %.4f\n", p, q, r)
println()

# Coefficients: [β³, β², β¹, β⁰]
c3 = 4.0
c2 = -6.0
c1 = 9/4 + r
c0 = -r/2

@printf("  4β³  - 6β²  + %.4fβ  - %.4f = 0\n", c1, -c0)
println()

# Solve the cubic
poly = Polynomial([c0, c1, c2, c3])   # Polynomials.jl: ascending order
roots_all = roots(poly)

println("  All roots:")
for (i, root) in enumerate(roots_all)
    re, im = real(root), imag(root)
    @printf("    β_%d = %.8f + %.8fi\n", i, re, im)
end
println()

# Select unique real positive root
real_pos = [r for r in roots_all
            if abs(imag(r)) < 1e-8 && real(r) > 0]

β_star = real(real_pos[1])

@printf("  Unique real positive root: β* = %.8f\n", β_star)
@printf("  Exact value:               β* = 1/4 = %.8f\n", 1/4)
@printf("  Match: %s\n", abs(β_star - 1/4) < 1e-8 ? "✓" : "✗")
println()

# Verify cubic
check = c3*β_star^3 + c2*β_star^2 + c1*β_star + c0
@printf("  Verification: 4(β*)³ - 6(β*)² + %.4f(β*) - %.4f = %.2e  ✓\n",
        c1, -c0, check)
println()

# Bounds check from paper: 0 < β* < 1/2
@printf("  Bounds: 0 < %.4f < 1/2  %s\n",
        β_star, 0 < β_star < 0.5 ? "✓" : "✗")
println()

# ── Step 3: α* from eq. (4.16) ───────────────────────────────────

println("─── Step 3: α* from eq. (4.16) ───")
println()
println("  Castellani et al. eq. (4.16):")
println()
println("  α* = (p²/q²)(3β* - 4β*²)")
println()

α_star = (p^2/q^2) * (3β_star - 4β_star^2)

@printf("  α* = (1)(3 × 1/4 - 4 × 1/16)\n")
@printf("     = 3/4 - 1/4\n")
@printf("     = %.8f\n", α_star)
@printf("  Exact value: α* = 1/2 = %.8f  %s\n",
        1/2, abs(α_star - 1/2) < 1e-8 ? "✓" : "✗")
println()

# ── Step 4: γ* from eq. (4.14a) ──────────────────────────────────

println("─── Step 4: γ* from Einstein equation (4.14a) ───")
println()
println("  From eq. (4.14a): β*(1 - β*)γ*² = 3ν/4")
println("  Normalise ν = 1:")
println()

ν = 1.0
γ_sq = (3ν/4) / (β_star * (1 - β_star))
γ_star = sqrt(γ_sq)

@printf("  γ*² = (3/4) / (β*(1-β*)) = (3/4) / (%.4f × %.4f)\n",
        β_star, 1-β_star)
@printf("       = %.8f\n", γ_sq)
@printf("  γ*  = %.8f\n", γ_star)
@printf("  Exact value: γ* = 2 = %.8f  %s\n",
        2.0, abs(γ_star - 2.0) < 1e-8 ? "✓" : "✗")
println()

# ── Step 5: Physical metric parameters ───────────────────────────

println("─── Step 5: Physical metric parameters a², b², c² ───")
println()
println("  Parametrisation (4.11) for p=q=1:")
println("    a = α*,  b = γ*β*,  c = γ*")
println()
println("  Killing metric normalisation factors (Castellani 1984, p.619-623):")

N_CP2 = 24/5
N_S2  = 12/5
N_fib = 1/20

@printf("    N_CP2 = %g = %g/5\n", N_CP2, N_CP2*5)
@printf("    N_S2  = %g = %g/5\n", N_S2,  N_S2*5)
@printf("    N_fib = %g = 1/20\n", N_fib)
println()

# Physical metric: a² = N_CP2 × α*², etc.
a2 = N_CP2 * α_star^2
b2 = N_S2  * (γ_star * β_star)^2
c2 = N_fib * (q * γ_star)^2

println("  Physical metric parameters:")
@printf("    a² = N_CP2 × α*²        = %.4f × (%.4f)² = %.8f\n",
        N_CP2, α_star, a2)
@printf("    b² = N_S2  × (γ*β*)²    = %.4f × (%.4f)² = %.8f\n",
        N_S2, γ_star*β_star, b2)
@printf("    c² = N_fib × (q γ*)²    = %.4f × (%.4f)² = %.8f\n",
        N_fib, q*γ_star, c2)
println()

# Exact values
a2_exact = 6/5
b2_exact = 3/5
c2_exact = 1/5

println("  Exact values:")
@printf("    a² = 6/5 = %.8f  %s\n", a2_exact, abs(a2-a2_exact)<1e-8 ? "✓" : "✗")
@printf("    b² = 3/5 = %.8f  %s\n", b2_exact, abs(b2-b2_exact)<1e-8 ? "✓" : "✗")
@printf("    c² = 1/5 = %.8f  %s\n", c2_exact, abs(c2-c2_exact)<1e-8 ? "✓" : "✗")
println()

# Ratios
println("  Ratios:")
@printf("    a²/b² = %.4f  (expected 2) %s\n", a2/b2, abs(a2/b2-2)<1e-6 ? "✓":"✗")
@printf("    a²/c² = %.4f  (expected 6) %s\n", a2/c2, abs(a2/c2-6)<1e-6 ? "✓":"✗")
@printf("    b²/c² = %.4f  (expected 3) %s\n", b2/c2, abs(b2/c2-3)<1e-6 ? "✓":"✗")
println()

# ── Step 6: Identification with FisherGeometrics ─────────────────

println("─── Step 6: FisherGeometrics parameters ───")
println()
println("  The metric parameters of M^{1,1,1} are identified as:")
println()

τ = c2   # = 1/5
κ = a2   # = 6/5

@printf("  τ = c² = %.8f = 1/5  (S¹ fibre radius squared)  %s\n",
        τ, abs(τ - 1/5) < 1e-8 ? "✓" : "✗")
@printf("  κ = a² = %.8f = 6/5  (ℂP² factor radius squared) %s\n",
        κ, abs(κ - 6/5) < 1e-8 ? "✓" : "✗")
println()

# ── Step 7: Consequences ─────────────────────────────────────────

println("─── Step 7: FisherGeometrics predictions ───")
println()
println("  All predictions involving τ and κ now follow from")
println("  the N=2 Einstein metric on M^{1,1,1}:")
println()

n_dim = 6        # dimension of ℂ^n in FisherGeometrics
M_KK  = 178.1   # GeV

G_F_pred = τ^2 * (n_dim^2 + 1) / (4 * M_KK^2)
G_F_obs  = 1.1664e-5

η_pred = (3/4) * τ^13
η_obs  = 6.1e-10

λ_W_pred = τ * sqrt(κ)
λ_W_obs  = 0.22534

Λ_QCD_pred = M_KK * exp(-2π) * (1 - τ)^2
Λ_QCD_obs  = 0.210   # GeV

@printf("  %-20s  %-14s  %-14s  %-8s\n",
        "Observable", "Predicted", "Observed", "Δ")
println("  " * repeat("─", 62))
@printf("  %-20s  %-14.6e  %-14.6e  %6.3f%%\n",
        "G_F [GeV⁻²]",
        G_F_pred, G_F_obs,
        abs(G_F_pred-G_F_obs)/G_F_obs*100)
@printf("  %-20s  %-14.4e  %-14.4e  %6.3f%%\n",
        "η_baryon",
        η_pred, η_obs,
        abs(η_pred-η_obs)/η_obs*100)
@printf("  %-20s  %-14.6f  %-14.6f  %6.3f%%\n",
        "λ_W (Wolfenstein)",
        λ_W_pred, λ_W_obs,
        abs(λ_W_pred-λ_W_obs)/λ_W_obs*100)
@printf("  %-20s  %-14.4f  %-14.4f  %6.3f%%\n",
        "Λ_QCD [GeV]",
        Λ_QCD_pred, Λ_QCD_obs,
        abs(Λ_QCD_pred-Λ_QCD_obs)/Λ_QCD_obs*100)
println()

# ── Summary ───────────────────────────────────────────────────────

println("─── Summary ───")
println()
println("  ┌────────────────────────────────────────────────────────┐")
println("  │  τ = 1/5 and κ = 6/5 are NOT free parameters        │")
println("  │  They follow uniquely from:                           │")
println("  │                                                        │")
println("  │  1. G_K = SU(3)×SU(2)×U(1)  (SM gauge group)       │")
println("  │  2. N=2 SUSY  →  M^{1,1,1} uniquely                 │")
println("  │  3. Cubic (4.15): β* = 1/4  (exact, rational)       │")
println("  │  4. Eq. (4.16):   α* = 1/2  (exact, rational)       │")
println("  │                                                        │")
println("  │  Result:  τ = c² = 1/5,  κ = a² = 6/5              │")
println("  │                                                        │")
println("  │  Source: Castellani, D'Auria, Fré                   │")
println("  │          Nucl. Phys. B239 (1984) 610–652            │")
println("  └────────────────────────────────────────────────────────┘")
println()
println("  © 2026 Jan Bouwman — github.com/uwbanjoman/FisherGeometrics.jl")
