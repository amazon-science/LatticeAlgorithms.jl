# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
  
# Licensed under the Apache License, Version 2.0 (the "License").
# You may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

################ Exact MLD by BSV for unrotated surface codes

# First define some aux functions

# An aux function for bsv(ws)
function SimulateHorizontal(j, M, logΓ, ws)
    d = size(M, 1) ÷ 2
    s = zeros(d)
    t = zeros(d)
    for i in 1 : d
        e = (j-1) * (2d-1) + i
        # println("SimulateHorizontal, $j, $e")
        logΓ += log10((1+ws[e]^2)/2)
        t[i] = (1-ws[e]^2)/(1+ws[e]^2)
        s[i] = 2ws[e] / (1+ws[e]^2)
    end

    big_t = transpose(hcat(t, zeros(length(t))))[:][1:end-1]
    big_s = transpose(hcat(s, s))[:]

    A = Matrix(Tridiagonal(-big_t, zeros(2d), big_t))
    B = diagm(big_s)

    logΓ += 1/2 * log10(det(M + A))
    M = A - B * inv(M + A) * B

    return M, logΓ

end

# An aux function for bsv(ws)
function SimulateVertical(j, M, logΓ, ws)
    d = size(M, 1) ÷ 2
    s = zeros(d-1)
    t = zeros(d-1)
    for i in 1 : d-1
        e = (j-1) * (2d-1) + d + i
        # println("SimulateVertical, $j, $e")
        logΓ += log10(1+ws[e]^2)
        s[i] = (1-ws[e]^2)/(1+ws[e]^2)
        t[i] = 2ws[e] / (1+ws[e]^2)
    end

    big_t = transpose(hcat(zeros(length(t)), t))[:]
    push!(big_t, 0)
    big_s = transpose(hcat(s, s))[:]
    pushfirst!(big_s, 1)
    push!(big_s, 1)

    A = Matrix(Tridiagonal(-big_t, zeros(2d), big_t))
    B = diagm(big_s)

    logΓ += 1/2 * log10(det(M + A))
    M = A - B * inv(M + A) * B

    return M, logΓ

end


"""
    bsv(ws::Vector{Float64})

return the expectation value of a matching gate quantum circuit using the BSV method

Notes:
    The expectation value is defined as 
        Z(ws) = log(Σ_{g∈G}∏_{i=1}^N (ws[i])^{g[i]}
    where N is the number of qubits, G is the full stabilizer group of the unrotated surface 
    code, and ws is the weight of the data qubits. We note that g should be treated as an 
    N-component vector of 0 and 1.
    
    This expectation value is shown to be equivalent to a matching gate
    quantum circuit with parameters ws, and can be evaluated efficiently with the BSV method.
    See Sec V C in https://arxiv.org/pdf/1405.4883.pdf.

    The log in Z(ws) is included such that the numbers are not too small. Also to get the MLD
    of the unrotated surface, we need to include 
        logπf ≡ sum([fe==1 ? -log10(wse) : 0 for (fe, wse) in zip(f, ws)])
    where f is the errant data qubits
"""
function bsv(ws::Vector)
    N = length(ws) # The number of qubits in the unrotated surface code N = d^2 + (d-1)^2
    d = (√(2N-1) + 1)/2 # The distance of the unrotated surface code

    if !isinteger(d)
        error("The length of `ws` has to be d^2+(d-1)^2 for an integer d")
    end    
    d = round.(Int, d)

    M0 = Matrix(0I, 2d, 2d)
    M0[2:2d-1, 2:2d-1] = BlockDiagonal([[0 1; -1 0] for _ in 1 : d-1])
    M0[1, 2d] = 1
    M0[2d, 1] = -1

    M = copy(M0)

    logΓ = (d-1) * log10(2)

    for j in 1 : d
        M, logΓ = SimulateHorizontal(j, M, logΓ, ws)

        # QR decomposition for numerical stability
        Q, R = qr(M)
        R = diagm(sign.(diag(R)))
        M = (Q * R - transpose(Q * R))/2
        if j < d
            M, logΓ = SimulateVertical(j, M, logΓ, ws)

            # QR decomposition for numerical stability
            Q, R = qr(M)
            R = diagm(sign.(diag(R)))
            M = (Q * R - transpose(Q * R))/2
        end
    end

    Zw = 1/2 * (logΓ - log10(2)) + 1/4 * log10(det(M + M0))
    return Zw
end

