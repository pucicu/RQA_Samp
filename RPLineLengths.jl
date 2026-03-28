using Random
using Base.Threads


"""
    get_hist_diagonal_RP(x, ε)

Compute the histogram of diagonal line lengths in a recurrence plot.

# Arguments
- `x::AbstractMatrix`: Input data matrix where rows represent observations and columns represent variables.
- `ε::Real`: Recurrence threshold for determining if two state vectors are recurrent.

# Returns
- `L::Vector{Int}`: Histogram where `L[n]` counts the number of diagonal lines of length `n` in the recurrence plot.

# Description
The function:
1. Constructs the recurrence plot `R`, marking pairs of points as recurrent if their Euclidean distance is ≤ `e`.
2. Counts diagonal line segments in lower triangle of `R`, storing their length frequencies in `L`.
"""
function get_hist_diagonal_RP(x, ε)
    N, dim = size(x)                  # Number of observations and variables
    L = zeros(Int, N)                 # Histogram for line lengths
    R = falses(N, N)                  # Boolean recurrence plot
    
    # Step 1: Construct recurrence plot
    for i in 1:N
        for j in 1:i
            D2 = 0.0                  # Squared distance accumulator
            for k in 1:dim  # Loop over dimensions
                D2 += (x[i, k] - x[j, k])^2
            end
            R[i, j] = sqrt(D2) <= ε   # Mark recurrence if within threshold
        end
    end

    # Step 2: Count diagonal line lengths
    for i in 2:N # start with 2 to remove main diagonal
        cnt = 0
        for j in 1:(N - i + 1)
            if R[i + j - 1, j]
                cnt += 1              # Extend current diagonal
            else
                if cnt > 0
                    L[cnt] += 1       # Store completed line length
                end
                cnt = 0
            end
        end
        if cnt > 0
            L[cnt] += 1
        end
    end

    return L
end


"""
    get_hist_diagonal_woRP(x, ε)

Compute the histogram of diagonal line lengths directly from a time series without 
explicitly building the recurrence plot.

# Arguments
- `x::AbstractMatrix`: Input data matrix where rows are observations and columns are variables.
- `ε::Real`: Recurrence threshold (maximum allowed Euclidean distance for recurrence).

# Returns
- `L::Vector{Int}`: Histogram where `L[n]` counts the number of diagonal lines of length `n`.

# Description
This version avoids constructing the full recurrence plot matrix, computing line 
lengths directly by comparing shifted segments of the input. It iterates through 
offsets and counts consecutive recurrent points to build the histogram.
"""
function get_hist_diagonal_woRP(x, ε)
    N, dim = size(x)                  # Number of observations and variables
    L = zeros(Int, N)                 # Histogram for line lengths
    ε2 = ε^2                          # Precompute squared threshold (avoids sqrt)

    @inbounds for i in 2:N # start with 2 to remove main diagonal
        cnt = 0
        @inbounds for j in 1:(N - i + 1)
            # Direct squared Euclidean distance computation

            # variant 1 (needs to replace test D2 < ε2 with D < ε)
#             @views diff .= x[i + j - 1, :] .- x[j, :]  # Avoid creating a new array
#             #diff .= x[i + j - 1, :] .- x[j, :]  ## Causes many allocations
#             D = norm(diff)  # Reuse `diff` instead of allocating new memory
            
            # variant 2 (needs to replace test D2 < ε2 with D < ε)
            # D = norm(x[i + j - 1, :] - x[j, :])
            
            # variant 3
            D2 = 0.0
            @inbounds for k in 1:dim
                D2 += (x[i + j - 1, k] - x[j, k])^2
            end

            if D2 <= ε2
                cnt += 1              # Extend current diagonal
            else
                if cnt > 0
                    L[cnt] += 1       # Store completed line length
                    cnt = 0
                end
            end
        end
        if cnt > 0
            L[cnt] += 1
        end
    end

    return L
end


"""
    get_hist_diagonal_sampled(x, ε, M)

Compute an approximate histogram of diagonal line lengths by random sampling,  
reducing computation cost compared to a full scan.

# Arguments
- `x::AbstractMatrix`: Input data matrix where rows are observations and columns are variables.
- `ε::Real`: Recurrence threshold (maximum Euclidean distance for recurrence).
- `M::Int`: Number of valid diagonal lines to sample.

# Returns
- `L::Vector{Int}`: Histogram where `L[n]` counts the number of sampled diagonal lines of length `n`.
- `N::Int`: Total number of searches.

# Description
This version speeds up diagonal line length computation by:
1. Randomly selecting starting indices `i` instead of checking all.
2. Stopping each search as soon as a line is found (`break`), avoiding unnecessary calculations.
3. Repeating until `M` lines have been successfully found.

The result is an **approximate** histogram that is much faster for large datasets,  
at the cost of reduced completeness compared to full enumeration.
"""
function get_hist_diagonal_sampled(x::AbstractMatrix{T}, ε::T, M::Int) where {T<:AbstractFloat}
    N, dim = size(x)                  # Number of observations and variables
    L_local = zeros(Int, N)           # Histogram for line lengths
    L_local_old = zeros(Int, N)       # Histogram for line lengths
    ε2 = ε^2                          # Squared threshold (avoids sqrt)
    count = 0                         # Number of valid lines found
    countAll = 0                      # Number of searches
    total_pairs = N * (N - 1) ÷ 2     # Number of (i,j) paires with i > j (excl. LOI)

    det_history = Float64[]
    l_history   = Float64[]

    # Convergenz parameter
    K = 50              # newly calculate all K lines
    W = 5               # window size for sliding average
    tol = 1e-3          # relative toleranz
    min_samples = 500   # minimum number of samples before break

    sum_n_L  = 0.0   # sum(n * L[n]) for all n >= 1
    sum_n2_L = 0.0   # sum(n * L[n]) for n >= 2
    sum_L2   = 0.0   # sum(L[n]) for n >= 2
    
    while count < M
        countAll += 1                 # Count number of searches
        idx = rand(1:total_pairs)     # Random start pair (i,j) in linear notation
        i_start = ceil(Int, (1 + sqrt(1 + 8*idx)) / 2)     # Translate linear index to i
        j_start = idx - (i_start - 1) * (i_start - 2) ÷ 2  # Translate linear index to j

        #println("i: ", i_start, "  j: ", j_start)

        # Check if R(i_start,j_start) = 1 (start point)
        D2 = zero(T)
        @inbounds for k in 1:dim
            D2 += (x[i_start,k] - x[j_start,k])^2
        end
        if D2 > ε2
            continue  # kein Linienstart, nächster Versuch
        end

        # Check if preceeding point R(i_start-1, j_start-1) = 0 (beginning of a line)
        if i_start == 1 || j_start == 1
            # We are at the lower border – no previous points
        else
            D2_prev = 0.0
            @inbounds for k in 1:dim
                D2_prev += (x[i_start-1,k] - x[j_start-1,k])^2
            end
            if D2_prev <= ε2
                continue  # no begin of a line, try again
            end
        end

        # Line found, now count line length
        cnt = 0
        @inbounds for offset in 0:(N - i_start)
            D2_line = 0.0
            for k in 1:dim
                D2_line += (x[i_start + offset,k] - x[j_start + offset,k])^2
            end
            if D2_line <= ε2          # Count points belonging to diagonal
                cnt += 1              # Extend diagonal
            else
                break                 # Line ends
            end
        end

        # Store line length to histogram variable
        if cnt > 0 #&& !((i_start, j_start) in seen)
            L_local[cnt] += 1         # Store completed line length
            count += 1                # Count found lines (only counted if a line was found)

            sum_n_L  += cnt
            sum_n2_L += cnt >= 2 ? cnt : 0
            sum_L2   += cnt >= 2 ? 1  : 0



            # im Loop, nach count += 1:
            if count % K == 0 && count >= min_samples

                # DET und L direkt aus den Summen:
                det_val = sum_n2_L / sum_n_L
                l_val   = sum_n2_L / sum_L2

                push!(det_history, det_val)
                push!(l_history,   l_val)


                if length(det_history) >= W
                    # relative change over last window
                    det_window = det_history[end-W+1:end]
                    l_window   = l_history[end-W+1:end]

                    rel_det = (maximum(det_window) - minimum(det_window)) / mean(det_window)
                    rel_l   = (maximum(l_window)   - minimum(l_window))   / mean(l_window)

                    # Break if difference is not changing for last K iterations
                    if rel_det < tol && rel_l < tol
                        println("Converged after $count samples (rel_det=$rel_det, rel_l=$rel_l)")
                        break
                    end
                end
            end
        end
    
    end

    return L_local, countAll
end


"""
    rqa(P, N)
Compute Recurrence Quantification Analysis (RQA) measures from a diagonal line length histogram.

# Arguments
- `P::Vector{Int}`: Histogram where `P[n]` counts the number of diagonal lines of length `n`.
- `N::Int`: Total number of points in the upper (or lower) triangle of the recurrence matrix (number of searches using the sampling approach).

# Returns
- `RR::Float64`: Recurrence Rate – fraction of recurrence points in the recurrence matrix.
- `DET::Float64`: Determinism – fraction of recurrence points forming diagonal lines of length ≥ 2.
- `L::Float64`: Average diagonal line length (lines of length ≥ 2).
- `ENTR::Float64`: Shannon entropy of the diagonal line length distribution.

# Description
RQA extracts dynamical invariants from the recurrence plot:
- **RR** measures overall recurrence density.
- **DET** measures predictability/determinism of the system.
- **L** relates to the mean prediction horizon.
- **ENTR** measures the complexity of the diagonal line structure.
"""
function rqa(P, N)
    p = P / sum(P)                    # Normalized probability distribution of line lengths

    idx = findall(!iszero, P)         # Indices of non-zero histogram entries

    # Recurrence Rate: ratio of recurrent points to total points
    # Numerator:   total number of points in all lines (length * count)
    # Denominator: total matrix points = N (isolated) + points in lines of length > 1
    RR = sum((1:length(P)) .* P) / (N + sum(idx .* (P[idx] .- 1)))

    # Determinism: fraction of recurrent points in lines of length >= 2
    DET = sum((2:length(P)) .* P[2:end]) / sum((1:length(P)) .* P)

    # Average line length: mean length of diagonal lines with length >= 2
    L = sum((2:length(P)) .* P[2:end]) / sum(P[2:end])

    # Shannon entropy of the line length distribution (only non-zero probabilities)
    ENTR = -sum(pi > 0 ? pi * log(pi) : 0 for pi in p)

    return RR, DET, L, ENTR
end
