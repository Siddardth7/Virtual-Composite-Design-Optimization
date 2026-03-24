import numpy as np
import copy
import random

# ----------------------------
# Material and Ply Definitions
# ----------------------------

# Material properties (Pa)
E1 = 89e9       # Longitudinal modulus
E2 = 8e9        # Transverse modulus
G12 = 4.5e9     # Shear modulus
nu12 = 0.3      # Major Poisson's ratio
nu21 = nu12 * E2 / E1  # Minor Poisson's ratio

denom = 1 - nu12 * nu21
Q11 = E1 / denom
Q22 = E2 / denom
Q12 = nu12 * E2 / denom
Q66 = G12

# Ply thickness (m) and allowed range for full laminate
t = 0.001       # 1 mm per ply, adjust as needed

# Allowed ply orientations (degrees)
allowed_angles = [0, 45, -45, 90]

# Minimum percentage (of full laminate) for each orientation type
minPercent = 0.1  # at least 10% plies in each required orientation

# Minimum and Maximum number of plies for the full (symmetric) laminate.
# (Because we optimize only half, min_half and max_half refer to the half-laminate.)
min_full_plies = 4     # full laminate minimum (=> min_half = 2)
max_full_plies = 40    # full laminate maximum (=> max_half = 20)
min_half = min_full_plies // 2
max_half = max_full_plies // 2

# ----------------------------
# Classical Laminate Theory Functions
# ----------------------------

def compute_Qbar(angle_deg):
    """Compute the transformed reduced stiffness matrix Qbar for a given ply orientation."""
    theta = np.radians(angle_deg)
    m = np.cos(theta)
    n = np.sin(theta)
    m2 = m**2
    n2 = n**2
    m4 = m**4
    n4 = n**4
    m3n = m**3 * n
    mn3 = m * n**3
    m2n2 = m2 * n2

    Qbar11 = Q11 * m4 + 2*(Q12 + 2*Q66) * m2n2 + Q22 * n4
    Qbar12 = (Q11 + Q22 - 4*Q66) * m2n2 + Q12 * (m4 + n4)
    Qbar22 = Q11 * n4 + 2*(Q12 + 2*Q66) * m2n2 + Q22 * m4
    Qbar16 = (Q11 - Q12 - 2*Q66) * m3n - (Q22 - Q12 - 2*Q66) * mn3
    Qbar26 = (Q11 - Q12 - 2*Q66) * mn3 - (Q22 - Q12 - 2*Q66) * m3n
    Qbar66 = (Q11 + Q22 - 2*Q12 - 2*Q66) * m2n2 + Q66 * (m4 + n4)

    return np.array([[Qbar11, Qbar12, Qbar16],
                     [Qbar12, Qbar22, Qbar26],
                     [Qbar16, Qbar26, Qbar66]])

def compute_D_matrix(full_seq):
    """
    Compute the bending stiffness matrix D (3x3) for a symmetric laminate using CLT.
    full_seq: array (or list) of ply orientations for the full laminate.
    The total thickness is len(full_seq)*t.
    """
    N = len(full_seq)
    total_thickness = N * t
    # z-coordinates from bottom (-h/2) to top (h/2)
    z = np.linspace(-total_thickness/2, total_thickness/2, N+1)
    D = np.zeros((3,3))
    for i in range(N):
        Qbar = compute_Qbar(full_seq[i])
        D += Qbar * (z[i+1]**3 - z[i]**3) / 3.0
    return D

def laminate_deflection_metric(full_seq):
    """
    Compute an objective metric based on bending stiffness.
    For bending about the 0° axis, we use D11 from the D matrix.
    A higher D11 corresponds to a stiffer laminate.
    """
    D = compute_D_matrix(full_seq)
    D11 = D[0, 0]
    if D11 <= 0:
        return 1e12  # Penalty for nonphysical stiffness
    return 1.0 / D11

def laminate_weight(full_seq):
    """Compute laminate weight as number of plies times a per-ply weight."""
    ply_weight = 0.005  # kg per ply (example value)
    return len(full_seq) * ply_weight

def laminate_objective(full_seq):
    """
    Overall objective: combine a deflection metric and weight.
    Lower values are better. Penalty terms are added if constraints are violated.
    """
    # Trade-off: lower deflection metric and lower weight are desired.
    return laminate_deflection_metric(full_seq) + laminate_weight(full_seq) + penalty(full_seq)

# ----------------------------
# Constraint Penalty Functions
# ----------------------------

def penalty(full_seq):
    """
    Apply penalty terms if the laminate does not meet symmetry, balance, or minimum ply requirements.
    For symmetry, the full laminate must be symmetric.
    For balance, the count of +45 must equal that of -45.
    Also, require a minimum percentage of plies for 0°, 90°, and (45 or -45).
    """
    pen = 0
    N = len(full_seq)
    # Symmetry: for a symmetric laminate, first half equals reversed second half.
    half = N // 2
    if not np.allclose(full_seq[:half], full_seq[half:][::-1]):
        pen += 1e6

    # Balance: number of 45° vs. -45° must be equal.
    count_p45 = np.sum(np.array(full_seq) == 45)
    count_m45 = np.sum(np.array(full_seq) == -45)
    pen += 1e6 * abs(count_p45 - count_m45)

    # Minimum ply percentages for each orientation type.
    required = int(N * minPercent)
    count_0 = np.sum(np.array(full_seq) == 0)
    count_90 = np.sum(np.array(full_seq) == 90)
    count_45_total = count_p45 + count_m45

    if count_0 < required:
        pen += 1e6 * (required - count_0)
    if count_90 < required:
        pen += 1e6 * (required - count_90)
    if count_45_total < required:
        pen += 1e6 * (required - count_45_total)
    return pen

def objective(half_seq):
    """
    Build the full symmetric stacking sequence from the half sequence,
    then compute the overall objective value.
    """
    full_seq = half_seq + half_seq[::-1]
    return laminate_objective(full_seq)

# ----------------------------
# Simulated Annealing Optimizer (Variable-Length)
# ----------------------------

def simulated_annealing(n_iterations=10000, initial_temp=1.0, cooling_rate=0.999):
    """
    Optimize both the ply sequence (orientations) and the total number of plies.
    The design variable is the half-laminate (which is then mirrored).
    Moves include:
      - Changing an orientation.
      - Inserting a ply (if below maximum half length).
      - Deleting a ply (if above minimum half length).
    """
    # Initialize with a random half-sequence of random length between min_half and max_half.
    current_length = random.randint(min_half, max_half)
    current_half = [random.choice(allowed_angles) for _ in range(current_length)]
    current_obj = objective(current_half)
    best_half = copy.deepcopy(current_half)
    best_obj = current_obj
    T = initial_temp

    for it in range(n_iterations):
        # Determine allowed move types based on current half length.
        moves = ['change']
        if len(current_half) < max_half:
            moves.append('insert')
        if len(current_half) > min_half:
            moves.append('delete')
        move = random.choice(moves)
        
        candidate_half = current_half.copy()
        
        if move == 'change':
            # Change orientation of one random ply.
            idx = random.randint(0, len(candidate_half) - 1)
            # Choose a new angle different from the current one.
            candidate_half[idx] = random.choice([a for a in allowed_angles if a != candidate_half[idx]])
        
        elif move == 'insert':
            # Insert a new ply at a random position.
            idx = random.randint(0, len(candidate_half))
            candidate_half.insert(idx, random.choice(allowed_angles))
        
        elif move == 'delete':
            # Delete a ply at a random position.
            idx = random.randint(0, len(candidate_half) - 1)
            candidate_half.pop(idx)
        
        candidate_obj = objective(candidate_half)
        
        # Accept candidate move if improved or probabilistically.
        delta = candidate_obj - current_obj
        if delta < 0 or random.random() < np.exp(-delta / T):
            current_half = candidate_half
            current_obj = candidate_obj
            if current_obj < best_obj:
                best_half = copy.deepcopy(current_half)
                best_obj = current_obj

        # Cool down
        T *= cooling_rate

        # Optional: print progress every 1000 iterations
        if (it + 1) % 1000 == 0:
            print(f"Iteration {it+1}: Best Obj = {best_obj:.3e}, Half sequence length = {len(best_half)}")

    best_full_seq = best_half + best_half[::-1]
    return best_full_seq, best_obj

# ----------------------------
# Run the Optimization
# ----------------------------

if __name__ == "__main__":
    best_seq, best_value = simulated_annealing(n_iterations=10000, initial_temp=1.0, cooling_rate=0.999)
    print("\nOptimized Layup Sequence (degrees):")
    print(best_seq)
    print("\nNumber of plies (full laminate):", len(best_seq))
    print("\nObjective Value (deflection metric + weight + penalties):", best_value)
